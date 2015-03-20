{-# LANGUAGE OverloadedStrings #-}

import Plover.Types
import Plover.Macros
import Plover.Compile (writeProgram)

--dotcfile = "#include <stdio.h>\n#include \"plover.h\"\nvoid plover_test() { printf(\"plover test success! " ++ show fundamentalConstant ++ "\\n\"); }\n"


dothfile = "void plover_test();"

plover_test_sig = FD { fd_params = [], fd_output = Void}
plover_test = FnDeclare "plover_test" plover_test_sig $ seqList [
  "printf" :$ (Str $ "plover test success! " ++ show 22 ++ "\n")
 ]

main = do
  writeProgram "plover.c" (return plover_test)
  writeFile "plover.h" dothfile
