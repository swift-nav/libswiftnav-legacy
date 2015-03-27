{-# LANGUAGE OverloadedStrings #-}

import Plover.Types
import Plover.Macros
import Plover.Compile (compileLib)

--import qualified Plover as P (main)

import PVT

-- TODO remove
dothfile = "void plover_test();"
plover_test_sig = FnT [] [] Void
plover_test = FnDef "plover_test" plover_test_sig $ seqList [
  "printf" :$ (Str $ "plover test success! " ++ show 22 ++ "\n")
 ]

main = do
  compileLib "generated" [pvtDef]
