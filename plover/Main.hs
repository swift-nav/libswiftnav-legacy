{-# LANGUAGE OverloadedStrings #-}
import System.Environment

import Plover.Types
import Plover.Compile

import qualified AmbiguityTest
import qualified AmbKF
import qualified Utils

units :: [CompilationUnit]
units = [ AmbiguityTest.cu
        , AmbKF.cu
        ]

main = do
  args <- getArgs
  case args of
    [cdir, hdir] -> do
        mapM_ (generateMain hdir cdir) units
    _ -> error "Requires exactly two arguments: C and H file output directories"

