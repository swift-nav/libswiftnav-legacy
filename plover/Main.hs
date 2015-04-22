{-# LANGUAGE OverloadedStrings #-}
import System.Environment

import Plover.Types
import Plover.Compile

mainBody :: CExpr
mainBody =
  Ref "printf" :$ StrLit "hello world\n" :> Return 0

defs :: [FunctionDefinition]
defs = [("helloWorld", [VoidExpr], FnT [] [] IntType, mainBody)]

logGenerate :: String -> IO ()
logGenerate fn = putStrLn $ "generating file: " ++ fn

main = do
  args <- getArgs
  case args of
    [c_dir, h_dir] -> do
        files <- generate "generated" c_dir h_dir ["stdio.h", "plover/generated.h"] defs
        mapM_ logGenerate files
    _ -> error "Requires exactly one argument: output file name"

