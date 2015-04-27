{-# LANGUAGE OverloadedStrings #-}

module AmbKF (cu) where

import Plover.Types
import Plover.Macros
import Plover.Compile

import qualified Utils
{-
typedef struct {
  u32 state_dim;
  u32 obs_dim;
  double amb_drift_var;
  double decor_mtx[MAX_OBS_DIM * MAX_OBS_DIM]; //the decorrelation matrix. takes raw measurements and decorrelates them
  double decor_obs_mtx[MAX_STATE_DIM * MAX_OBS_DIM]; //the observation matrix for decorrelated measurements
  double decor_obs_cov[MAX_OBS_DIM]; //the diagonal of the decorrelated observation covariance (for cholesky is ones)
  double null_basis_Q[(MAX_STATE_DIM - 3) * MAX_OBS_DIM];
  double state_mean[MAX_STATE_DIM];
  double state_cov_U[MAX_STATE_DIM * MAX_STATE_DIM];
  double state_cov_D[MAX_STATE_DIM];
} nkf_t;
-}

nkf_t = StructDecl "nkf_t" (
            ST External [
              ("state_dim", IntType),
              ("obs_dim", IntType),
              ("amb_drift_var", NumType),
              ("decor_mtx", VecType ["kf" :-> "obs_dim", "kf" :-> "obs_dim"] NumType),
              ("decor_obs_mtx", VecType ["kf" :-> "state_dim", "kf" :-> "obs_dim"] NumType),
              ("decor_obs_cov", VecType ["kf" :-> "obs_dim"] NumType),
              ("null_basis_Q", VecType ["kf" :-> "state_dim" - 3, "kf" :-> "obs_dim"] NumType),
              ("state_mean", VecType ["kf" :-> "state_dim"] NumType),
              ("state_cov_U", VecType ["kf" :-> "state_dim", "kf" :-> "state_dim"] NumType),
              ("state_cov_D", VecType ["kf" :-> "state_dim"] NumType)
              ]
            )

simpleAmbMeasurement :: FunctionDefinition
simpleAmbMeasurement = ("simple_amb_measurement", FnT []
                  [ ("carrier", NumType)
                  , ("code", NumType)
                  ] NumType, body)
  where
    body :: CExpr
    body = Utils.constants :> (Return $ Utils.simpleAmbMeas "carrier" "code")

incorporateScalarMeasurementT = FnT [("state_dim", IntType)]
                                [ ("h", VecType ["state_dim"] NumType)
                                , ("R", NumType)
                                , ("U", VecType ["state_dim", "state_dim"] NumType)
                                , ("D", VecType ["state_dim"] NumType)
                                , ("k", VecType ["state_dim"] NumType)
                                ] Void

-- TODO: Make use of triangular nature of U
{-
incorporateScalarMeasurement :: FunctionDefinition
incorporateScalarMeasurement = ("incorporate_scalar_measurement2",
                                FnT [("state_dim", IntType)]
                                [ ("h", VecType ["state_dim"] NumType)
                                , ("R", NumType)
                                , ("U", VecType ["state_dim", "state_dim"] NumType)
                                , ("D", VecType ["state_dim"] NumType)
                                , ("k", VecType ["state_dim"] NumType)
                                ] Void, body)
  where
    body :: CExpr
    body = seqList $ [
      "f" := transpose "U" * "h",
      "g" := "D" * "f",
      "alpha" := "f" `dot` "g" + "R",

      Return 0
      ]
-}

{-
  for (u32 i=0; i<kf->obs_dim; i++) {
    double *h = &kf->decor_obs_mtx[kf->state_dim * i]; /* vector of length kf->state_dim. */
    double R = kf->decor_obs_cov[i]; /* scalar. */
    double k[kf->state_dim]; /*  vector of length kf->state_dim. */

    /* updates cov and sets k. */
    incorporate_scalar_measurement(kf->state_dim, h, R, kf->state_cov_U, kf->state_cov_D, &k[0]);

    double predicted_obs = 0;
    /* TODO take advantage of sparsity of h. */
    for (u32 j=0; j<kf->state_dim; j++) {
      predicted_obs += h[j] * kf->state_mean[j];
    }
    double obs_minus_predicted_obs = decor_obs[i] - predicted_obs;

    for (u32 j=0; j<kf->state_dim; j++) {
      kf->state_mean[j] += k[j] * obs_minus_predicted_obs; /* uses k to update mean. */
    }
  }
-}

incorporateObs :: FunctionDefinition
incorporateObs = ("incorporate_obs2", FnT []
                  [ ("kf", PtrType $ TypedefType "nkf_t")
                  , ("decor_obs", VecType ["kf" :-> "obs_dim"] NumType)
                  ] Void, body)
  where
    body :: CExpr
    body = seqList $ [
      Extern "incorporate_scalar_measurement" $ FnType incorporateScalarMeasurementT,
      Vec "i" ("kf" :-> "obs_dim") $ seqList [
        Declare (VecType ["kf" :-> "state_dim"] NumType) "k",
        App "incorporate_scalar_measurement"
          [ (transpose $ "kf" :-> "decor_obs_mtx") :! "i"
          , ("kf" :-> "decor_obs_cov") :! "i"
          , "kf" :-> "state_cov_U"
          , "kf" :-> "state_cov_D"
          , "k"
          ]
        ],
      Return 0
      ]


cu :: CompilationUnit
cu = CU
  { unitName = "amb_kf_plover"
  , sourceDefs = [simpleAmbMeasurement, incorporateObs]
  , sourceIncs = ["stdio.h", "constants.h", "math.h", "amb_kf.h", "plover/amb_kf_plover.h"]
  , headerDefs = [nkf_t]
  , headerIncs = ["common.h", "amb_kf.h"]
  }

