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

cu :: CompilationUnit
cu = CU
  { unitName = "ambiguity_test"
  , sourceDefs = [helloPlover]
  , sourceIncs = ["stdio.h", "plover/ambiguity_test.h"]
  , headerDefs = [testStruct]
  , headerIncs = ["common.h"]
  }

