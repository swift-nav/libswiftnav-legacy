{-# LANGUAGE OverloadedStrings #-}

module AmbiguityTest (defs, includes) where

import Plover.Types

mainBody :: CExpr
mainBody =
  Ref "printf" :$ StrLit "hello world\n" :> Return 0

defs :: [FunctionDefinition]
defs = [("helloWorld", [VoidExpr], FnT [] [] IntType, mainBody)]

includes :: [String]
includes = ["stdio.h", "plover/ambiguity_test.h"]

