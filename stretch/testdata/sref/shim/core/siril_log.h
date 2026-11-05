/* Stand-in for Siril core/siril_log.h: logging is a no-op in the oracle. */
#ifndef SREF_SHIM_SIRIL_LOG_H
#define SREF_SHIM_SIRIL_LOG_H
#define siril_log_message(...) ((void)0)
#define siril_log_color_message(...) ((void)0)
#define siril_debug_print(...) ((void)0)
#endif
