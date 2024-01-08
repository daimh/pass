#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include "sorter.h"
sorter::sorter(char *seq, FILE *out, const char *tmpfile, uint64_t mem_size,
	bool bisulfite, uint64_t io_size
) {
	_M_seq = seq;
	if (!out) utl_error("sorter", 1, "Output file is null\n");
	_M_out = out;
	if (!tmpfile) utl_error("sorter", 2, "tmpfile is null\n");
	if (!*tmpfile) utl_error("sorter", 3, "tmpfile is empty\n");
	_M_s_tmp = tmpfile;
	_M_f_tmp = NULL;
	_M_sort_size = mem_size / sizeof(sequence) / 2;
	_M_buf = (sequence*) malloc(_M_sort_size * sizeof(sequence));
	if (!_M_buf) utl_error("sorter", 4, "memory allocation failed\n");
	_M_cpy = (sequence*) malloc(_M_sort_size * sizeof(sequence));
	if (!_M_cpy) utl_error("sorter", 5, "memory allocation failed\n");
	_M_bisulfite = bisulfite;
	_M_io_size = io_size;
	_M_total_num = _M_sort_tail = 0;
	_M_bin_len = 0;
}
sorter::~sorter() {
	FILE **bin_arr = NULL;
	if (_M_bin_len) {
		fsync(fileno(_M_f_tmp));
		fclose(_M_f_tmp);
		bin_arr = (FILE**) malloc(sizeof(FILE*) * _M_bin_len);
		for (int i=0; i<_M_bin_len; i++) {
			bin_arr[i] = fopen(_M_s_tmp, "r");
			if (!bin_arr[i])
				utl_error("sorter", 6, "cannot open tmpfile for read, check 'ulimit -n'\n");
			if (fseek(bin_arr[i], sizeof(sequence) * _M_sort_size * i, SEEK_SET))
				utl_error("sorter", 7, "cannot fseek tmpfile for read\n");
		}
	}
	_M_buf = (sequence*) realloc(_M_buf, _M_sort_tail * sizeof(sequence));
	_M_cpy = (sequence*) realloc(_M_cpy, _M_sort_tail * sizeof(sequence));
	_M_total_num += _M_sort_tail;
	if (fwrite(&_M_total_num, sizeof(uint64_t), 1, _M_out) != 1) 
		utl_error("sorter", 8, "cannot write to index file!\n");
	if (_M_bisulfite)
		_M_qsort_bisulfite(_M_buf, _M_cpy, 0, 0, _M_sort_tail);
	else
		_M_qsort_normal(_M_buf, _M_cpy, 0, 0, _M_sort_tail);
	if (_M_bin_len) {
		//merge them into one big file
		uint64_t *bin_idx = (uint64_t*) calloc(sizeof(uint64_t), _M_bin_len+1);
		_M_io_size /= _M_bin_len + 1;
		if (_M_io_size < MIN_IO_SIZE) _M_io_size = MIN_IO_SIZE;
		sequence **bin_buf = (sequence**) malloc(sizeof(sequence*) * _M_bin_len);
		if (!bin_buf) utl_error("sorter", 9, "memory allocation failed\n");
		for (int i=0; i<_M_bin_len; i++) {
			bin_buf[i] = (sequence*) malloc(sizeof(sequence) * _M_io_size);
			if (!bin_buf[i]) utl_error("sorter", 10, "memory allocation failed\n");
		}
		sequence::buf_4_less_than = _M_seq;
		sequence::bisulfite = _M_bisulfite;
		sorter_node<sequence>* root = NULL;
		for (int i=0; i<=_M_bin_len; i++) root = _M_add_node(root, bin_idx, bin_buf, bin_arr, i);
		sequence* writer_buf = (sequence*) malloc(sizeof(sequence) * _M_io_size);
		uint64_t writer_idx = 0;
		while (root) {
			int from;
			root = root->get_min(writer_buf + writer_idx++, from);
			if (writer_idx == _M_io_size) {
				if (fwrite(writer_buf, sizeof(sequence) * _M_io_size, 1, _M_out) != 1)
					utl_error("sorter", 11, "cannot write to output file\n");
//for (uint64_t i=0; i<writer_idx; i++) writer_buf[i].DEBUG();
				writer_idx = 0;
			}
			root = _M_add_node(root, bin_idx, bin_buf, bin_arr, from);
		}
		if (writer_idx)
			if (fwrite(writer_buf, sizeof(sequence) * writer_idx, 1, _M_out) != 1)
				utl_error("sorter", 12, "cannot write to output file\n");
//for (uint64_t i=0; i<writer_idx; i++) writer_buf[i].DEBUG();
		free(writer_buf);
		free(bin_idx);
		for (int i=0; i<_M_bin_len; i++) free(bin_buf[i]);
		free(bin_buf);
		//merge is done.
		for (int i=0; i<_M_bin_len; i++) fclose(bin_arr[i]);
		free(bin_arr);
		unlink(_M_s_tmp);
	} else {
		if (fwrite(_M_buf, sizeof(sequence), _M_sort_tail, _M_out) != _M_sort_tail) 
			utl_error("sorter", 13, "cannot write to index file!\n");
	}
	free(_M_buf);
	free(_M_cpy);
}
void sorter::append(const sequence& val) {
	if (_M_sort_tail == _M_sort_size) {
		if (!_M_f_tmp) {
			_M_f_tmp = fopen(_M_s_tmp, "w");
			if (!_M_f_tmp) utl_error("sorter", 14, "fopen %s failed\n", tmpfile);
		}
		if (_M_bisulfite)
			_M_qsort_bisulfite(_M_buf, _M_cpy, 0, 0, _M_sort_tail);
		else
			_M_qsort_normal(_M_buf, _M_cpy, 0, 0, _M_sort_tail);
		_M_bin_len++;
		if (fwrite(_M_buf, sizeof(sequence) * _M_sort_tail, 1, _M_f_tmp) != 1)
			utl_error("sorter", 15, "cannot write to tmporary file\n");
		_M_total_num += _M_sort_tail;
		_M_sort_tail = 0;
	}
	_M_buf[_M_sort_tail++] = val;
}
void sorter::_M_qsort_normal(sequence *src, sequence *dst, short level, uint64_t head,
	uint64_t tail
) {
	uint64_t dst_head = head;
	uint64_t dst_tail = tail;
	uint64_t startC = 0;
	for (uint64_t i=head; i<tail; i++)
		switch (src[i].getchar(_M_seq, level)) {
			case 0:
				memcpy(dst + dst_head++, src + i, sizeof(sequence));
				break;
			case 'A':
				startC++;
				break;
			case 'T':
				memcpy(dst + --dst_tail, src + i, sizeof(sequence));
				break;
		}
	uint64_t startA = dst_head;
	uint64_t startT = dst_tail;
	startC += startA;
	for (uint64_t i=head; i<tail; i++)
		switch (src[i].getchar(_M_seq, level)) {
			case 'A':
				memcpy(dst + dst_head++, src + i, sizeof(sequence));
				break;
			case 'C':
				memcpy(dst + startC++, src + i, sizeof(sequence));
				break;
			case 'G':
				memcpy(dst + --dst_tail, src + i, sizeof(sequence));
				break;
		}
	startC = dst_head;
	uint64_t startG = dst_tail;
	level++;
	if (startA + 1 < startC) _M_qsort_normal(dst, src, level, startA, startC);
	if (startC + 1 < startG) _M_qsort_normal(dst, src, level, startC, startG);
	if (startG + 1 < startT) _M_qsort_normal(dst, src, level, startG, startT);
	if (startT + 1 < tail) _M_qsort_normal(dst, src, level, startT, tail);
	if (src == _M_buf) {
		if (startA > head) memcpy(_M_buf + head, _M_cpy + head, sizeof(sequence) * (startA - head));
		if (startA + 1 == startC) memcpy(_M_buf + startA, _M_cpy + startA, sizeof(sequence));
		if (startC + 1 == startG) memcpy(_M_buf + startC, _M_cpy + startC, sizeof(sequence));
		if (startG + 1 == startT) memcpy(_M_buf + startG, _M_cpy + startG, sizeof(sequence));
		if (startT + 1 == tail) memcpy(_M_buf + startT, _M_cpy + startT, sizeof(sequence));
	}
}
void sorter::_M_qsort_bisulfite(sequence *src, sequence *dst, short level,
	uint64_t head, uint64_t tail
) {
	uint64_t dst_head = head;
	uint64_t dst_tail = tail;
	for (uint64_t i=head; i<tail; i++)
		switch (src[i].getchar(_M_seq, level)) {
			case 0:
				memcpy(dst + dst_head++, src + i, sizeof(sequence));
				break;
			case 'C':
			case 'T':
				memcpy(dst + --dst_tail, src + i, sizeof(sequence));
				break;
		}
	uint64_t startA = dst_head;
	uint64_t startCT = dst_tail;
	for (uint64_t i=head; i<tail; i++)
		switch (src[i].getchar(_M_seq, level)) {
			case 'A':
				memcpy(dst + dst_head++, src + i, sizeof(sequence));
				break;
			case 'G':
				memcpy(dst + --dst_tail, src + i, sizeof(sequence));
				break;
		}
	uint64_t startG = dst_tail;
	level++;
	if (startA + 1 < startG) _M_qsort_bisulfite(dst, src, level, startA, startG);
	if (startG + 1 < startCT) _M_qsort_bisulfite(dst, src, level, startG, startCT);
	if (startCT + 1 < tail) _M_qsort_bisulfite(dst, src, level, startCT, tail);
	if (src == _M_buf) {
		if (startA > head) memcpy(_M_buf + head, _M_cpy + head, sizeof(sequence) * (startA - head));
		if (startA + 1 == startG) memcpy(_M_buf + startA, _M_cpy + startA, sizeof(sequence));
		if (startG + 1 == startCT) memcpy(_M_buf + startG, _M_cpy + startG, sizeof(sequence));
		if (startCT + 1 == tail) memcpy(_M_buf + startCT, _M_cpy + startCT, sizeof(sequence));
	}
}
sorter_node<sequence>* sorter::_M_add_node(sorter_node<sequence>* root,
		uint64_t* bin_idx, sequence** bin_buf, FILE** bin_arr, int from
) {
	sequence* val;
	if (from < _M_bin_len)
		if (bin_idx[from] >= _M_sort_size)
			return root;
		else {
			uint64_t loc = bin_idx[from] % _M_io_size;
			if (!loc) {
				uint64_t len = _M_sort_size - bin_idx[from];
				if (len > _M_io_size) len = _M_io_size;
				if (fread(bin_buf[from], sizeof(sequence) * len, 1, bin_arr[from]) != 1)
					utl_error("sorter", 16, "Cannot read from tmpfile\n");
					
			}
			bin_idx[from]++;
			val = bin_buf[from] + loc;
		}
	else if (bin_idx[from] >= _M_sort_tail)
		return root;
	else
		val = _M_buf + bin_idx[from]++;
	return sorter_node<sequence>::insert(root, *val, from);
}
