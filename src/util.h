#ifndef SRC_UTIL_H_
#define SRC_UTIL_H_
#include <stdint.h>
void utl_error(const char *src, int err, const char *format, ...);
void check_read(bool b);
void check_write(bool b);
void check_alloc(void* p);
uint64_t get_vsize();
#endif
