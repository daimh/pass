#ifndef SRC_SORTER_H_
#define SRC_SORTER_H_
#include <algorithm>
#include <iostream>
#include "util.h"
#include "sorternode.h"
#include "sequence.h"
class sorter {
public:
	sorter(char* seq, FILE *out, const char *tmpfile, uint64_t mem_size, bool bisulfite,
		uint64_t io_size = 0);
	~sorter();
	void append(const sequence& val);
private:
	const static uint64_t MIN_IO_SIZE = 32768;
	FILE *_M_out, *_M_f_tmp;
	const char *_M_s_tmp;
	uint64_t _M_sort_size, _M_sort_tail, _M_io_size, _M_total_num;
	sequence *_M_buf, *_M_cpy;
	char* _M_seq;
	int _M_bin_len;
	bool _M_bisulfite;

	inline sorter_node<sequence>* _M_add_node(sorter_node<sequence>* root,
		uint64_t* bin_idx, sequence** bin_buf, FILE** bin_arr, int from);
	void _M_qsort_normal(sequence *src, sequence *dst, short level, uint64_t head,
		uint64_t tail);
	void _M_qsort_bisulfite(sequence *src, sequence *dst, short level, uint64_t head,
		uint64_t tail);
};
#endif
