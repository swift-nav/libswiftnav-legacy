{-# LANGUAGE OverloadedStrings #-}
import System.Environment

import Plover.Types
import Plover.Compile

import qualified AmbiguityTest

units :: [CompilationUnit]
units = [AmbiguityTest.cu]

main = do
  args <- getArgs
  case args of
    [c_dir, h_dir] -> do
        mapM_ (generateMain c_dir h_dir) units
    _ -> error "Requires exactly two arguments: C and H file output directories"

