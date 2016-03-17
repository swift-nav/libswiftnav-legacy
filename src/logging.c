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

#include <stdarg.h>

#include <libswiftnav/logging.h>

/** \defgroup logging Logging
 * Logging functions.
 * \{ */

const char *level_string[] = {
  "EMERGENCY",
  "ALERT",
  "CRITICAL",
  "ERROR",
  "WARNING",
  "NOTICE",
  "INFO",
  "DEBUG",
};

/** Log message by level.
 *
 * \param level Log level
 * \param file file name
 * \param line line number
 * \param msg Log contents
 */
void log_(u8 level, const char *file, const int line, const char *msg, ...)
{
  (void)file;
  (void)line;

  va_list ap;

  printf("%s: ", level_string[level]);
  va_start(ap, msg);
  vprintf(msg, ap);
  va_end(ap);
}

/* \} */
