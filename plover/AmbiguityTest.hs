{-# LANGUAGE OverloadedStrings #-}

module AmbiguityTest (cu) where

import Plover.Types
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
    body = Ref "printf" :$ StrLit "hello world\n" :> Return 0

cu :: CompilationUnit
cu = CU
  { unitName = "ambiguity_test"
  , definitions = [helloPlover]
  , includes = ["stdio.h", "plover/ambiguity_test.h"]
  , headerDefs = [testStruct]
  }

