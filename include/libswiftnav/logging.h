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
#include <string.h>

#include <libswiftnav/common.h>

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
 *    \#define DEBUG 1
 *
 * \{ */

#define __FNAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

void log_(u8 level, const char *file, const int line, const char *msg, ...) __attribute__ ((weak));

#define LOG_EMERG       0       /* system is unusable */
#define LOG_ALERT       1       /* action must be taken immediately */
#define LOG_CRIT        2       /* critical conditions */
#define LOG_ERROR       3       /* error conditions */
#define LOG_WARN        4       /* warning conditions */
#define LOG_NOTICE      5       /* normal but significant condition */
#define LOG_INFO        6       /* informational */
#define LOG_DEBUG       7       /* debug-level messages */

#define log_emerg(args...) log_(LOG_EMERG, __FNAME__, __LINE__, args)

/** Log an alert.
 * \param args `printf` style format and arguments.
 */
#define log_alert(args...) log_(LOG_ALERT, __FNAME__, __LINE__, args)

/** Log a critical event.
 * \param args `printf` style format and arguments.
 */
#define log_crit(args...) log_(LOG_CRIT, __FNAME__, __LINE__, args)

/** Log an error.
 * \param args `printf` style format and arguments.
 */
#define log_error(args...) log_(LOG_ERROR, __FNAME__, __LINE__, args)

/** Log a warning.
 * \param args `printf` style format and arguments.
 */
#define log_warn(args...)  log_(LOG_WARN, __FNAME__, __LINE__, args)

/** Log a notice.
 * \param args `printf` style format and arguments.
 */
#define log_notice(args...)  log_(LOG_NOTICE, __FNAME__, __LINE__, args)

/** Log an information message.
 * \param args `printf` style format and arguments.
 */
#define log_info(args...)  log_(LOG_INFO, __FNAME__, __LINE__, args)

/** Log a debugging message.
 * \param args `printf` style format and arguments.
 */
#define log_debug(args...)  \
do {                        \
  if (DEBUG) {              \
    log_(LOG_DEBUG, __FNAME__, __LINE__, args);  \
  }                         \
} while (0)

/** Log a debug message indicating entry to a function.
 * Logs a debug message of the form `\<function_name\>` to indicate entry to a
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
 * Logs a debug message of the form `\</function_name\>` to indicate exit from a
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
