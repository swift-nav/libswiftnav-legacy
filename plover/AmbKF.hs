{-# LANGUAGE OverloadedStrings #-}

module AmbKF (cu) where

import Plover.Types
import Plover.Macros
import Plover.Compile

import qualified Utils

simpleAmbMeasurement :: FunctionDefinition
simpleAmbMeasurement = ("simple_amb_measurement", FnT []
                  [ ("carrier", NumType)
                  , ("code", NumType)
                  ] NumType, body)
  where
    body :: CExpr
    body = Utils.constants :> (Return $ Utils.simpleAmbMeas "carrier" "code")

-- TODO: Make use of triangular nature of U
incorporateScalarMeasurement :: FunctionDefinition
incorporateScalarMeasurement = ("incorporate_scalar_measurement",
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

cu :: CompilationUnit
cu = CU
  { unitName = "amb_kf"
  , sourceDefs = [simpleAmbMeasurement]
  , sourceIncs = ["stdio.h", "constants.h", "math.h", "plover/amb_kf.h"]
  , headerDefs = []
  , headerIncs = ["common.h"]
  }

