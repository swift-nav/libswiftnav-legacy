{-# LANGUAGE OverloadedStrings #-}

module AmbiguityTest (cu) where

import Plover.Types
import Plover.Macros
import Plover.Compile

testStruct = StructDecl "test_struct_t" (
               ST Generated [
                 ("n", IntType),
                 ("xs", VecType [22] NumType)
                 ]
               )

helloPlover :: FunctionDefinition
helloPlover = ("hello_world", FnT [] [] IntType, body)
  where
    body :: CExpr
    body = seqList $ [
      Ref "printf" :$ StrLit "hello world\n",
      Return 0
      ]
{-
/*  Presumes that the first alm entry is the reference sat. */
s8 assign_de_mtx(u8 num_sats, const sdiff_t *sats_with_ref_first,
                 const double ref_ecef[3], double *DE)
{
  if (num_sats <= 1) {
    log_debug("assign_de_mtx: not enough sats\n");
    return -1;
  }

  assert(num_sats > 1);

  /* Vector to reference satellite. */
  double e_0[3];
  vector_subtract(3, sats_with_ref_first[0].sat_pos, ref_ecef, e_0);
  vector_normalize(3, e_0);

  for (u8 i=1; i<num_sats; i++) {
    /* Vector to satellite i */
    double e_i[3];
    vector_subtract(3, sats_with_ref_first[i].sat_pos, ref_ecef, e_i);
    vector_normalize(3, e_i);
    /* DE row = e_i - e_0 */
    vector_subtract(3, e_i, e_0, &DE[3*(i-1)]);
  }

  return 0;
}

typedef struct {
  double pseudorange;
  double carrier_phase;
  double doppler;
  double sat_pos[3];
  double sat_vel[3];
  double snr;
  u8 prn;
} sdiff_t;

-}
sdiff_t = StructDecl "sdiff_t" (
            ST External [
              ("pseudorange", NumType),
              ("carrier_phase", NumType),
              ("doppler", NumType),
              ("sat_pos", VecType [3] NumType),
              ("sat_vel", VecType [3] NumType),
              ("snr", NumType),
              ("prn", IntType)
              ]
            )

assignDE :: FunctionDefinition
assignDE = ("assign_de_mtxP", FnT [("num_sats", IntType)]
           [ ("sds", VecType ["num_sats"] (TypedefType "sdiff_t"))
           , ("ref_ecef", VecType [3] NumType)
           , ("DE", VecType ["num_sats"-1, 3] NumType)
           ] IntType, body)
  where
    body :: CExpr
    body = seqList [
      sdiff_t,
      "e_0" := (("sds" :! 0) :. "sat_pos" - "ref_ecef"),
      "e_0_" := normalize "e_0",
      "DE" :<= Vec "i" ("num_sats"-1) (e_i - "e_0"),
      Return 0
      ]
    e_i = normalize $ ("sds" :! "i") :. "sat_pos" - "ref_ecef"

cu :: CompilationUnit
cu = CU
  { unitName = "ambiguity_test"
  , sourceDefs = [helloPlover, assignDE]
  , sourceIncs = ["stdio.h", "single_diff.h", "math.h", "plover/ambiguity_test.h"]
  , headerDefs = [testStruct]
  , headerIncs = ["common.h"]
  }

