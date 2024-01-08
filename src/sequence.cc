#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "sequence.h"
char *sequence::buf_4_less_than = 0;
bool sequence::bisulfite = false;
sequence::sequence(uint64_t loc, short len, unsigned int idx) {
/*
	_M_storage = (uint64_t) len;
	_M_storage = _M_storage << 54;
	_M_storage = _M_storage | loc;
*/
	_M_storage = (((uint64_t)len) << 54) | loc;
	_M_idx = idx;
}
uint64_t sequence::get_loc() const {
	return _M_storage & 0x003FFFFFFFFFFFFFLLU;
}
short sequence::get_len() const {
	return (short)(_M_storage >> 54);
}
unsigned int sequence::get_idx() {
	return _M_idx;
}
char sequence::getchar(char *seq, short idx) {
	if (idx >= get_len())
		return 0;
	else
		return seq[get_loc() + idx];
}
bool sequence::operator< (const sequence& ano) const {
	short me_len = get_len();
	short ano_len = ano.get_len();
	short cmp_len = me_len > ano_len ? ano_len : me_len;
	int diff = 0;
	if (bisulfite) {
		for (short i=0; i<cmp_len; i++) {
			char me_c = buf_4_less_than[get_loc()];
			char ano_c = buf_4_less_than[ano.get_loc()];
			diff = me_c - ano_c;
			if (diff == 17 || diff == -17)
				diff = 0;
			else if (diff)
				break;
		}
	} else
		diff = memcmp(buf_4_less_than + get_loc(), buf_4_less_than + ano.get_loc(), cmp_len);
	if (diff)
		return diff < 0;
	else
		return me_len < ano_len;
}
void sequence::DEBUG(char *seq) {
	char *me = seq + get_loc();
	for (short i=0; i<get_len(); i++) fputc(me[i], stdout);
//	printf(":%lu:%d", get_loc(), get_len());
//	fputc('\n', stdout);
}
