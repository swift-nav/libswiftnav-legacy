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
    [cdir, hdir] -> do
        mapM_ (generateMain hdir cdir) units
    _ -> error "Requires exactly two arguments: C and H file output directories"

