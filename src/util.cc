#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
void utl_error(const char *src, int err, const char *format, ...) {
	fprintf(stderr, "ERR-%s-%03d: ", src, err);
	va_list args;
	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	exit(-1);
}
void check_write(bool b) {
        if (!b) utl_error("utl", 1, "cannot write to file, is disk full?\n");
}
void check_read(bool b) {
        if (!b) utl_error("utl", 2, "cannot read from file, is the file corrupted?\n");
}
void check_alloc(void* p) {
	if (p == NULL) utl_error("utl", 3, "cannot allocate memory");
}
uint64_t get_vsize() {
	uint64_t rtn = 0;
	FILE *f = fopen("/proc/self/status", "r");
	if (f) {
		size_t len = 0;
		char *line = NULL;
		while (getline(&line, &len, f) != -1) {
			if (strncmp(line, "VmSize:\t", 8) == 0) {
				rtn = strtol(line+8, NULL, 10) * 1024;
				break;
			}
		}
		if (line) free(line);
		fclose(f);
	}
	return rtn;
}
