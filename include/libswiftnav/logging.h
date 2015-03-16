/*
 * Copyright (C) 2015 Swift Navigation Inc.
 * Contact: Fergus Noble <fergus@swift-nav.com>
 *
 * This source is subject to the license found in the file 'LICENSE' which must
 * be be distributed together with this source. All other rights reserved.
 *
 * THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
 * EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
 */

#ifndef LIBSWIFTNAV_LOGGING_H
#define LIBSWIFTNAV_LOGGING_H

#include <stdio.h>

/* DEBUG off by default, enable it on a per-file basis. */
#ifndef DEBUG
#define DEBUG 0
#endif

/** \defgroup logging Logging
 * Logging
 *
 * Logging at the `DEBUG` level is turned off by default and should be enabled
 * on a per-file basis by adding the following line to the source file *before*
 * including `logging.h`:
 *
 *    #define DEBUG 1
 *
 * \{ */

/** Log an error.
 * \param args `printf` style format and arguments.
 */
#define log_error(args...) printf("ERROR: " args)

/** Log a warning.
 * \param args `printf` style format and arguments.
 */
#define log_warn(args...)  printf("WARNING: " args)

/** Log an information message.
 * \param args `printf` style format and arguments.
 */
#define log_info(args...)  printf("INFO: " args)

/** Log a debugging message.
 * \param args `printf` style format and arguments.
 */
#define log_debug(args...)  \
do {                        \
  if (DEBUG) {              \
    printf("DEBUG: " args); \
  }                         \
} while (0)

/** Log a debug message indicating entry to a function.
 * Logs a debug message of the form `<function_name>` to indicate entry to a
 * function. `function_name` is automatically filled in with the name of the
 * current function by GCC magic.
 */
#define DEBUG_ENTRY()              \
do {                               \
  if (DEBUG) {                     \
    log_debug("<%s>\n", __func__); \
  }                                \
} while (0)

/** Log a debug message indicating exit to a function.
 * Logs a debug message of the form `</function_name>` to indicate exit from a
 * function. `function_name` is automatically filled in with the name of the
 * current function by GCC magic.
 */
#define DEBUG_EXIT()                \
do {                                \
  if (DEBUG) {                      \
    log_debug("</%s>\n", __func__); \
  }                                 \
} while (0)

/** \} */

#endif /* LIBSWIFTNAV_LOGGING_H */

