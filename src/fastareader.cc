#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "util.h"
#include "fastareader.h"
#define BUFSIZE 1048576
fasta_reader* fasta_reader::_G_fa_ref;
fasta_reader::fasta_reader(char *fa_name, bool is_ref, char *id_from, char *id_to, FILE *fout) {
	check_alloc(_M_counter_other = (uint64_t*) calloc(256, sizeof(uint64_t)));
	_M_is_ref = is_ref;
	_M_id_num = 0;
	_M_id_loc = NULL;;
	_M_id_buf = NULL;
	_M_seq_loc = NULL;
	uint64_t seq_buf_size = BUFSIZE;
	check_alloc(_M_seq_buf = (char*)malloc(seq_buf_size));
	uint64_t seq_buf_tail = 0;
	uint64_t id_buf_size = 0, id_size = 0, id_buf_tail = 0;

	FILE *file_fa = fopen(fa_name, "r");
	if (!file_fa) utl_error("fastareader", 1, "cannot open file '%s'\n", fa_name);
	int file_type = fgetc(file_fa);
	if (file_type == EOF) utl_error("fastareader", 2, "cannot read from file '%s'\n", fa_name);
	ungetc(file_type, file_fa);
	size_t line_len;
	int64_t line_read;
	char *line = NULL;
	uint64_t row = 0;
	_M_has_score = file_type == '@';
	if (_M_has_score) {
		check_alloc(_M_score_buf = (unsigned char*)malloc(seq_buf_size));
		int64_t seq_len = 0;
		while ((line_read = getline(&line, &line_len, file_fa)) > 0) {
			row += 1;
			if (line_read > BUFSIZE)
				utl_error("fastareader", 3,
					"line of fasta file '%s' has more than %d character\n",
					fa_name, BUFSIZE);
			while (line_read > 0 && (line[line_read-1] == '\r' || line[line_read-1] == '\n')) line_read--;
			if (line_read == 0)
				utl_error("fastareader", 4,
					"Line %lu of fastq file '%s' is empty.\n",
					row, fa_name);
			line[line_read] = 0;
			int mod = row % 4;
			if (mod == 1) 
				_M_put_id(line, id_from, id_to, seq_buf_tail, &id_buf_tail, &id_buf_size, &id_size);
			else if (mod == 2) {
				seq_len = line_read;
				_M_put_seq(line, line_read, &seq_buf_tail, &seq_buf_size);
				_M_seq_buf[seq_buf_tail++] = 0;
			} else if (mod == 3) {
				if (strncmp(line, "+", 1))
					utl_error("fastareader", 5,
						"Line %lu of fastq file '%s' is not '+'\n", row, fa_name);
			} else if (mod == 0) {
				if (line_read != seq_len)
					utl_error("fastareader", 6,
						"Length of line %lu is not equal to its preceding sequence'\n",
						row);
				memcpy(_M_score_buf + seq_buf_tail - line_read - 1, line, line_read);
//				for (int i=0; i<line_read; i++) _M_score_buf[seq_buf_tail - line_read - 1 + i] = line[i];
				_M_score_buf[seq_buf_tail - 1] = 0;
			}
		}
		if (row % 4 != 0)
			utl_error("fastareader", 7,
				"Line number of the file is not multiple of 4.\n");
		_M_seq_buf[seq_buf_tail++] = -1;
	} else {
		_M_score_buf = NULL;
		bool wait_seq = false;
		while ((line_read = getline(&line, &line_len, file_fa)) > 0) {
			row += 1;
			if (line_read > BUFSIZE)
				utl_error("fastareader", 8,
					"line of fasta file '%s' has more than %d character\n",
					fa_name, BUFSIZE);
			while (line_read > 0 && (line[line_read-1] == '\r' || line[line_read-1] == '\n')) line_read--;
			if (line_read == 0) continue;
			line[line_read] = 0;
			if (*line == '>') {
				_M_finish_seq_read(&seq_buf_tail, fa_name, wait_seq, line);
				wait_seq = true;
				_M_put_id(line, id_from, id_to, seq_buf_tail, &id_buf_tail, &id_buf_size, &id_size);
			} else {
				wait_seq = false;
				_M_put_seq(line, line_read, &seq_buf_tail, &seq_buf_size);
			}
		}
		_M_finish_seq_read(&seq_buf_tail, fa_name, wait_seq, line);
		_M_seq_buf[seq_buf_tail++] = -1;
	}
	fclose(file_fa);
	if (line) free(line);
	for (int i=0; i<256; i++) 
		if (_M_counter_other[i])
			fprintf(stderr, "%ju '%c' is converted to 'N'\n",
				_M_counter_other[i], (char)i);
	free(_M_counter_other);
	check_write(fwrite(&_M_id_num, sizeof(uint64_t), 1, fout) == 1);;
	check_write(fwrite(_M_id_loc, sizeof(uint64_t), _M_id_num, fout) == _M_id_num);;
	check_write(fwrite(_M_seq_loc, sizeof(uint64_t), _M_id_num, fout) == _M_id_num);;
	check_write(fwrite(&id_buf_tail, sizeof(uint64_t), 1, fout) == 1);;
	check_write(fwrite(_M_id_buf, sizeof(char), id_buf_tail, fout) == id_buf_tail);;
	check_write(fwrite(&seq_buf_tail, sizeof(uint64_t), 1, fout) == 1);;
	check_write(fwrite(_M_seq_buf, sizeof(char), seq_buf_tail, fout) == seq_buf_tail);;
	check_write(fwrite(&_M_has_score, sizeof(bool), 1, fout) == 1);;
	if (_M_has_score) 
		check_write(fwrite(_M_score_buf, sizeof(unsigned char), seq_buf_tail, fout) == seq_buf_tail);;
	_M_size = sizeof(uint64_t)* id_size * 2 + id_buf_size + seq_buf_size;
	if (_M_is_ref && (id_buf_tail >= _M_get_limit() || seq_buf_tail >= _M_get_limit())) 
			utl_error("fastareader", 9,
				"cannot generate reversed sequence, try to reduce size of fasta file '%s'\n",
				fa_name);
	_M_finish_new(id_buf_tail, seq_buf_tail);
}
fasta_reader::fasta_reader(bool is_ref, FILE *fin, int limit, unsigned char mis_thr,
	int quality_score_base, bool print_seq, bool print_all
) {
	_M_counter_other = NULL;
	_M_print_seq = print_seq;
	_M_print_all = print_all;
	_M_is_ref = is_ref;
	uint64_t id_buf_tail, seq_buf_tail;
	check_read(fread(&_M_id_num, sizeof(uint64_t), 1, fin) == 1);
	check_alloc(_M_id_loc = (uint64_t*) malloc(sizeof(uint64_t) * _M_id_num));
	check_alloc(_M_seq_loc = (uint64_t*) malloc(sizeof(uint64_t) * _M_id_num));
	check_read(fread(_M_id_loc, sizeof(uint64_t), _M_id_num, fin) == _M_id_num);;
	check_read(fread(_M_seq_loc, sizeof(uint64_t), _M_id_num, fin) == _M_id_num);
	check_read(fread(&id_buf_tail, sizeof(uint64_t), 1, fin) == 1);
	check_alloc(_M_id_buf = (char*) malloc(sizeof(char) * id_buf_tail));
	check_read(fread(_M_id_buf, sizeof(char), id_buf_tail, fin) == id_buf_tail);
	check_read(fread(&seq_buf_tail, sizeof(uint64_t), 1, fin) == 1);;
	check_alloc(_M_seq_buf = (char*) malloc(sizeof(char) * seq_buf_tail));
	check_read(fread(_M_seq_buf, sizeof(char), seq_buf_tail, fin) == seq_buf_tail);
	check_read(fread(&_M_has_score, sizeof(bool), 1, fin) == 1);;
	if (_M_has_score) {
		if (quality_score_base != -1) {
			check_alloc(_M_score_buf = (unsigned char*) malloc(sizeof(char) * seq_buf_tail));
			check_read(fread(_M_score_buf, sizeof(unsigned char), seq_buf_tail, fin) == seq_buf_tail);
			for (uint64_t i=0; i<seq_buf_tail; i++)
				if (_M_score_buf[i]) {
					if ((int) _M_score_buf[i] < quality_score_base)
						utl_error("fastareader", 10,
							"parameter --quality-score-base larger than the quality score in original fastq file.\n");
					_M_score_buf[i] -= quality_score_base;
				}
		} else {
			_M_score_buf = NULL;
			fseek(fin, seq_buf_tail, SEEK_CUR);
			_M_has_score = false;
		}
	} else {
		_M_score_buf = NULL;
		if (quality_score_base != -1)
			fprintf(stderr, "Original query file is not FASTQ file, --quality-score-base is ignored\n");
	}
	_M_size = sizeof(uint64_t)* _M_id_num * 2 + id_buf_tail + seq_buf_tail;
	_M_finish_new(id_buf_tail, seq_buf_tail);
	_M_limit = limit;
	_M_mis_thr = mis_thr;
	_M_match_cnt = (unsigned char*)calloc(sizeof(unsigned char), _M_id_num);
	if (limit > 0)
		_M_match_buf = (match_info**)calloc(sizeof(match_info*), _M_id_num);
	else 
		_M_match_buf = NULL;
}
void fasta_reader::_M_finish_new(uint64_t id_buf_tail, uint64_t seq_buf_tail) {
	if (_M_is_ref) {
		_M_size *= 2;
		uint64_t size = sizeof(uint64_t)* _M_id_num * 2;
		check_alloc(_M_id_loc = (uint64_t*)realloc(_M_id_loc, size));
		check_alloc(_M_seq_loc = (uint64_t*)realloc(_M_seq_loc, size + 1));
		check_alloc(_M_id_buf = (char*)realloc(_M_id_buf, id_buf_tail*2));
		seq_buf_tail--;
		check_alloc(_M_seq_buf = (char*)realloc(_M_seq_buf, seq_buf_tail * 2));
		memcpy(_M_id_buf+id_buf_tail, _M_id_buf, id_buf_tail);
		for (uint64_t i=0; i<_M_id_num; i++) {
			_M_id_loc[_M_id_num+i] = id_buf_tail + _M_id_loc[i];
			_M_seq_loc[_M_id_num+i] = seq_buf_tail + _M_seq_loc[i] - 1;
			uint64_t len = strlen(_M_seq_buf + _M_seq_loc[i]);
			for (uint64_t j=0; j<len; j++) {
				switch (_M_seq_buf[_M_seq_loc[i] + len - j - 1]) {
					case 'A':
						_M_seq_buf[_M_seq_loc[i] + seq_buf_tail + j - 1] = 'T';
						break;
					case 'C':
						_M_seq_buf[_M_seq_loc[i] + seq_buf_tail + j - 1] = 'G';
						break;
					case 'G':
						_M_seq_buf[_M_seq_loc[i] + seq_buf_tail + j - 1] = 'C';
						break;
					case 'T':
						_M_seq_buf[_M_seq_loc[i] + seq_buf_tail + j - 1] = 'A';
						break;
					case 'N':
						_M_seq_buf[_M_seq_loc[i] + seq_buf_tail + j - 1] = 'N';
						break;
					default:
						utl_error("fastareader", 11,
							"Contact developer please, errorcode('%d').\n",
							(int)_M_seq_buf[_M_seq_loc[i] + len - j - 1]);
				}
			}
		}
		seq_buf_tail = seq_buf_tail * 2 - 1;
		_M_seq_buf[seq_buf_tail] = -1;
		_M_id_num *= 2;
		_M_seq_loc[_M_id_num] = seq_buf_tail;
	} else {
		check_alloc(_M_seq_loc = (uint64_t*)realloc(_M_seq_loc, sizeof(uint64_t)*(_M_id_num+1)));
		_M_seq_loc[_M_id_num] = seq_buf_tail - 1;
	}
}
fasta_reader::~fasta_reader() {
	if (_M_seq_buf) free(_M_seq_buf);
	if (_M_id_buf) free(_M_id_buf);
	if (_M_seq_loc) free(_M_seq_loc);
	if (_M_id_loc) free(_M_id_loc);
	if (_M_limit > 0) {
		for (uint64_t i=0; i<_M_id_num; i++)
			if (_M_match_buf[i]) free(_M_match_buf[i]);
		free(_M_match_buf);
	}
	free(_M_match_cnt);
}
uint64_t fasta_reader::get_size() {
	uint64_t rtn = get_vsize();
	if (rtn == 0)
		return _M_size;
	else
		return rtn;
}
char* fasta_reader::get_seq_buffer() {
	return _M_seq_buf+1;
}
unsigned char* fasta_reader::get_score_buffer() {
	if (_M_score_buf)
		return _M_score_buf+1;
	else
		return NULL;
}
void fasta_reader::_M_finish_seq_read(uint64_t *seq_buf_tail, char *fa_name,
	bool wait_seq, char *line
) {
	if (wait_seq)
		utl_error("fastareader", 12,
			"Fasta file '%s' has no sequence following line '%s'\n",
			fa_name, line);
	_M_seq_buf[(*seq_buf_tail)++] = 0;
}
bool fasta_reader::ignorable(unsigned int qry_idx) {
//	return _M_limit > 0 && _M_match_cnt[qry_idx] > _M_limit;
	return _M_match_cnt[qry_idx] == UCHAR_MAX;
}
void fasta_reader::add_match(unsigned int qry_idx, unsigned int ref_idx,
	uint64_t ref_loc, int score, char* mismatch, unsigned char mis_idx
) {
	if (_M_limit < 0) {
		fputs(_M_id_buf + _M_id_loc[qry_idx], stdout);
		if (_M_print_seq) {
			fputc('\t', stdout);
			fputs(_M_seq_buf + _M_seq_loc[qry_idx], stdout);
		}
		fputc('\t', stdout);
		_G_fa_ref->print_ref_fasta_id(ref_idx, ref_loc);
		if (_M_has_score) {
			fputc('\t', stdout);
			printf("%d", score);
		}
		fputs(mismatch, stdout);
		fputc('\n', stdout);
		_M_match_cnt[qry_idx] = 1;
	} else if (_M_match_cnt[qry_idx] < _M_limit) {
		_M_match_buf[qry_idx] = (match_info*)realloc(_M_match_buf[qry_idx],
			sizeof(match_info) * (_M_match_cnt[qry_idx] + 1));
		_M_match_buf[qry_idx][_M_match_cnt[qry_idx]].ref_idx = ref_idx;
		_M_match_buf[qry_idx][_M_match_cnt[qry_idx]].ref_loc = ref_loc;
		_M_match_buf[qry_idx][_M_match_cnt[qry_idx]].score = score;
		_M_match_buf[qry_idx][_M_match_cnt[qry_idx]].mismatch = strdup(mismatch);
		_M_match_cnt[qry_idx]++;
	} else if (_M_match_cnt[qry_idx] == _M_limit) {
		if (_M_limit > 0) free(_M_match_buf[qry_idx]);
		_M_match_cnt[qry_idx] = UCHAR_MAX;
		if (_M_print_all) {
			fputs(_M_id_buf + _M_id_loc[qry_idx], stdout);
			if (_M_print_seq) {
				fputc('\t', stdout);
				fputs(_M_seq_buf + _M_seq_loc[qry_idx], stdout);
			}
			fputs("\tFT\n", stdout);
		}
	} else {
		fputs(_M_id_buf + _M_id_loc[qry_idx], stderr);
		fputc('\n', stderr);
		fprintf(stderr, "daimh100:interesting bug1,%u,%u\n", _M_match_cnt[qry_idx], _M_limit);
		exit(1);
	}
}
void fasta_reader::print_match() {
	for (uint64_t qry_idx=0; qry_idx<_M_id_num; qry_idx++) {
		if (_M_match_cnt[qry_idx] == UCHAR_MAX || _M_match_cnt[qry_idx] == 0) continue;
		if (_M_limit > 0) {
			for (int j=0; j<_M_match_cnt[qry_idx]; j++) {
				fputs(_M_id_buf + _M_id_loc[qry_idx], stdout);
				if (_M_print_seq) {
					fputc('\t', stdout);
					fputs(_M_seq_buf + _M_seq_loc[qry_idx], stdout);
				}
				fputc('\t', stdout);
				_G_fa_ref->print_ref_fasta_id(_M_match_buf[qry_idx][j].ref_idx,
					_M_match_buf[qry_idx][j].ref_loc);
				if (_M_has_score) {
					fputc('\t', stdout);
					printf("%d", _M_match_buf[qry_idx][j].score);
				}
				fputs(_M_match_buf[qry_idx][j].mismatch, stdout);
				fputc('\n', stdout);
			}
			free(_M_match_buf[qry_idx]);
		}
		_M_match_cnt[qry_idx] = UCHAR_MAX;
	}
}
void fasta_reader::print_nomatch() {
	if (!_M_print_all) return;
	for (uint64_t qry_idx=0; qry_idx<_M_id_num; qry_idx++) {
		if (_M_match_cnt[qry_idx] == 0) {
			fputs(_M_id_buf + _M_id_loc[qry_idx], stdout);
			if (_M_print_seq) {
				fputc('\t', stdout);
				fputs(_M_seq_buf + _M_seq_loc[qry_idx], stdout);
			}
			fputs("\tNM\n", stdout);
		}
	}
}
int64_t fasta_reader::get_ref_loc(unsigned int ref_idx, uint64_t ref_start, uint64_t ref_end) {
	if (ref_idx >= _M_id_num >> 1)
		return ref_end + 2 - _M_seq_loc[ref_idx+1];
	else
		return ref_start + 2 - _M_seq_loc[ref_idx];
}
void fasta_reader::print_ref_fasta_id(unsigned int ref_idx, int64_t ref_loc) {
	fputs(_M_id_buf + _M_id_loc[ref_idx], stdout);
	if (ref_loc < 0) {
		fputs("\t-\t", stdout);
		printf("%jd", -ref_loc);
	} else {
		fputs("\t+\t", stdout);
		printf("%jd", ref_loc);
	}
}
uint64_t fasta_reader::_M_get_limit() {
//10 bit reserved for short sequence length
//54 bit reserved for merged sequences length
	return _M_is_ref ? 0x001FFFFFFFFFFFFFLLU : 0x003FFFFFFFFFFFFFLLU;
}
void fasta_reader::_M_put_id(char *line, char *id_from, char *id_to,
	uint64_t seq_buf_tail, uint64_t *id_buf_tail, uint64_t *id_buf_size,
	uint64_t *id_size
) {
	char *id = line+1;
	if (id_from && *id_from) {
		id = strstr(id, id_from);
		if (id == NULL)
			fprintf(stderr, "Warning: Cannot find start tag in '%s' in file\n", line);
		else
			id += strlen(id_from);
	}
	if (id_to && *id_to) {
		char *tail = strstr(id, id_to);
		if (tail == NULL)
			fprintf(stderr, "Warning: Cannot find end tag in '%s' in file\n", line);
		else
			*tail = 0;
	}
	uint64_t id_len = strlen(id) + 1;
	if (*id_buf_tail + id_len > *id_buf_size) {
		*id_buf_size += BUFSIZE;
		if (*id_buf_size >= _M_get_limit())
			utl_error("fastareader", 13,
				"ID string buffer is full, try to reduce size of the file\n");
		check_alloc(_M_id_buf = (char*)realloc(_M_id_buf, *id_buf_size));
	}
	memcpy(_M_id_buf + *id_buf_tail, id, id_len);
	if (_M_id_num >= *id_size) {
		*id_size += BUFSIZE;
		if (*id_size >= _M_get_limit())
			utl_error("fastareader", 14,
				"ID array buffer is full, try to reduce size of the file\n");
		check_alloc(_M_id_loc = (uint64_t*)realloc(_M_id_loc, sizeof(uint64_t)*(*id_size)));
		check_alloc(_M_seq_loc = (uint64_t*)realloc(_M_seq_loc, sizeof(uint64_t)*(*id_size)));
	}
	_M_id_loc[_M_id_num] = *id_buf_tail;
	_M_seq_loc[_M_id_num++] = seq_buf_tail;
	*id_buf_tail += id_len;
}
void fasta_reader::_M_put_seq(char *line, int64_t line_read, uint64_t *seq_buf_tail,
	uint64_t *seq_buf_size
) {
	if (*seq_buf_tail + line_read + 1 >= *seq_buf_size) {
		*seq_buf_size += BUFSIZE;
		if (*seq_buf_size >= _M_get_limit())
			utl_error("fastareader", 15,
				"Sequence buffer is full, try to reduce size of the file\n");
		check_alloc(_M_seq_buf = (char*)realloc(_M_seq_buf, *seq_buf_size));
		if (_M_has_score)
			check_alloc(_M_score_buf = (unsigned char*)realloc(_M_score_buf, *seq_buf_size));
	}
	for (int64_t i=0; i<line_read; i++) {
		switch (line[i]) {
			case 'a':
			case 'c':
			case 'g':
			case 'n':
			case 't':
				line[i] -= 32;
				break;
			case 'A':
			case 'C':
			case 'G':
			case 'N':
			case 'T':
				break;
			case 'u':
			case 'U':
				line[i] = 'T';
				break;
			default:
//				_M_counter_other[(unsigned char)line[i]]++;
				line[i] = 'N';
		}
	}
	memcpy(_M_seq_buf + *seq_buf_tail, line, line_read);
	*seq_buf_tail += line_read;
}
