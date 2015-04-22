{-# LANGUAGE OverloadedStrings #-}
import System.Environment

import Plover.Types
import Plover.Compile

import qualified AmbiguityTest (defs, includes)

files :: [([FunctionDefinition], String, [String])]
files = [
  (AmbiguityTest.defs, "ambiguity_test", AmbiguityTest.includes)
  ]

gen :: String -> String -> ([FunctionDefinition], String, [String]) -> IO ()
gen c_dir h_dir (defs, name, includes) = generate name c_dir h_dir includes defs

main = do
  args <- getArgs
  case args of
    [c_dir, h_dir] -> do
        mapM_ (gen c_dir h_dir) files
    _ -> error "Requires exactly two arguments: C and H file output directories"

