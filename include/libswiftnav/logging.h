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

/** \defgroup logging Logging
 * Logging
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
#define log_debug(args...) printf("DEBUG: " args)

/** \} */

#endif /* LIBSWIFTNAV_LOGGING_H */
