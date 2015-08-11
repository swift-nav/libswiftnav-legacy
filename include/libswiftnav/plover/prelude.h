#ifndef PLOVER_GENERATED_prelude
#define PLOVER_GENERATED_prelude



#include <stdint.h>
#include <stdbool.h>
#include <inttypes.h>
#ifndef COMMON_INT_TYPES
#define COMMON_INT_TYPES
/** \defgroup common_inttypes Integer types
 * Specified-width integer type definitions for shorter and nicer code.
 *
 * These should be used in preference to unspecified width types such as
 * `int` which can lead to portability issues between different platforms.
 * \{ */

/** Signed 8-bit integer. */
typedef int8_t s8;
/** Signed 16-bit integer. */
typedef int16_t s16;
/** Signed 32-bit integer. */
typedef int32_t s32;
/** Signed 64-bit integer. */
typedef int64_t s64;
/** Unsigned 8-bit integer. */
typedef uint8_t u8;
/** Unsigned 16-bit integer. */
typedef uint16_t u16;
/** Unsigned 32-bit integer. */
typedef uint32_t u32;
/** Unsigned 64-bit integer. */
typedef uint64_t u64;

#endif

/** \} */

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int ipow(int base, int exp);
double dipow(double base, int exp);
double rand_uniform(void);
double rand_normal (void);
double norm (const s32 n, const double * v);
void normalize (const s32 n, const double * v, double * result);
void print_vec (const s32 n, const double * v);
void print_mat (const s32 n, const s32 m, const double * A);
s32 pl_matrix_inverse (const s32 n, const double * A, double * B);
double det (const s32 n, const double * A);


#endif /* PLOVER_GENERATED_prelude */
