/* -*- C -*-  (not really, but good for syntax highlighting) */

%typemap(in) s8, s16, s32, s64 {
    $1 = ($1_type)PyInt_AsLong($input);
}

%typemap(out) s8, s16, s32, s64 {
    $result = PyInt_FromLong($1);
}

%typemap(in) u8, u16, u32, u64 {
    $1 = ($1_type)PyLong_AsUnsignedLong($input);
}

%typemap(out) u8, u16, u32, u64  {
    $result = PyLong_FromUnsignedLong($1);
}

%include "carrays.i"
%include "cpointer.i"

%array_class(s8, carray_s8)
%array_class(s16, carray_s16)
%array_class(s32, carray_s32)
%array_class(s64, carray_s64)

%array_class(u8, carray_u8)
%array_class(u16, carray_u16)
%array_class(u32, carray_u32)
%array_class(u64, carray_u64)

%array_class(float, carray_float)
%array_class(double, carray_double)
