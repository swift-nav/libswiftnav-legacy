{-# LANGUAGE OverloadedStrings #-}

module PVT where

import Control.Monad.Free

import Plover

-- The PVT Example --
decls :: CExpr
decls = seqList [
  Ext "GPS_OMEGAE_DOT" numType,
  Ext "GPS_C" numType 
 ]

losLoop :: CExpr
losLoop = Lam "j" (R "n_used") $ seqList [
  "tau" := norm ("rx_state" - "sat_pos" :! "j") / "GPS_C",
  "we_tau" := "GPS_OMEGAE_DOT" * "tau",
  -- TODO rewrite issue forces this onto its own line
  "rot" := rot_small "we_tau",
  "xk_new" := "rot" * ("sat_pos" :! "j"),
  --"xk_new" := rot_small "we_tau" * ("sat_pos" :! "j"),
  "xk_new" - "rx_state"
 ]

pvtSig = FnT
  { ft_imp = [("n_used", IntType)]
  , ft_exp =
      [("sat_pos", ExprType [R "n_used", 3])
      ,("pseudo", ExprType [R "n_used"])
      ,("rx_state", ExprType [3])
      ,("correction", ExprType [4])
      ,("G", ExprType ["n_used", 4])
      ,("X", ExprType [4, "n_used"])
      ]
  , ft_out = Void
  }

pvtBody = seqList [
    decls,
    "los" :=  losLoop,
    "G" :< Lam "j" (R "n_used") (normalize ((- "los") :! "j") :# (Lam "i" 1 1)),
    "omp" := "pseudo" - Lam "j" (R "n_used") (norm ("los" :! "j")),
    "X" :< inverse (transpose "G" * "G") * transpose "G",
    "correction" :< "X" * "omp"
 ]

pvtDef :: (Variable, FunctionType CExpr, CExpr)
pvtDef = ("pvt", pvtSig, pvtBody)

pvt :: CExpr
pvt = FnDef "pvt" pvtSig $ pvtBody

testPVT = do
  -- Generate random arguments, call "pvt" defined above
  test1 <- generateTestArguments "pvt" pvtSig
  -- Print n_used
  let pnused = ("printInt" :$ "n_used")
  -- Call the wrapped libswiftnav version
  let test2 = Free (App (R "pvt2") (map (R . fst) (ft_exp pvtSig)))
  n <- freshName
  let printer = Lam n 4 ("printDouble" :$ ("correction" :! R n))
  -- Definition of pvt, then main that calls test code
  return
    $ Ext "pvt2" (FnType $ pvtSig )
    :> pvt
    :> (wrapMain $ seqList
         [ test1
         , pnused
         , newline "generated output:"
         , printer
         , test2
         , newline "reference output:"
         , printer
         , newline ""
         ])

