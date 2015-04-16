{-# LANGUAGE OverloadedStrings #-}
import System.Environment

import Plover.Types
import Plover.Compile (defaultMain)

mainBody :: CExpr
mainBody =
  Ref "printf" :$ StrLit "hello world\n" :> Return 0

defs :: [FunctionDefinition]
defs = [("helloWorld", [], FnT [] [] IntType, mainBody)]

main = do
  args <- getArgs
  case args of
    [fn] -> do
      putStrLn $ "generating file: " ++ fn
      defaultMain fn [] defs
    _ -> error "Requires exactly one argument: output file name"
