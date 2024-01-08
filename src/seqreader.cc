#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "util.h"
#include "seqreader.h"
#include "sequence.h"
#define UNIT 33554432
int seq_reader::_M_mismatch_idx = 0;
char* seq_reader::_M_mismatch_buffer = (char*) malloc(1048576);
seq_reader::seq_reader(fasta_reader *fasta, FILE *fin, const char *filename) {
	_M_fin = fin;
	_M_fasta = fasta;
	_M_filename = filename;
	_M_ftell = ftell(fin);
	check_alloc(_M_ss_buf = (sequence*)malloc(sizeof(sequence) * UNIT));
}
seq_reader::~seq_reader() {
	free(_M_ss_buf);
}
void seq_reader::init_ref() {
	fseek(_M_fin, _M_ftell, SEEK_SET);
	_M_seed_idx = 0;
	check_read(fread(&_M_seed_cnt, sizeof(uint64_t), 1, _M_fin) == 1);
	_M_ss_size = fread(_M_ss_buf, sizeof(sequence), _M_read_size(UNIT), _M_fin);
	_M_seed_idx += _M_ss_size;
	_M_ss_head = 0;
}
void seq_reader::jump_to(int step) {
	uint64_t seed_max_len;
	for (int i=0; i<step; i++) {
		check_read(fread(&seed_max_len, sizeof(uint64_t), 1, _M_fin) == 1);
		check_read(fread(&_M_seed_cnt, sizeof(uint64_t), 1, _M_fin) == 1);
		if (fseek(_M_fin, sizeof(sequence) * _M_seed_cnt, SEEK_CUR)) 
			utl_error("seqreader", 1, "query sequence is corrupted\n");
	}
}
uint64_t seq_reader::get_max_seed_len(uint64_t ref_max_len) {
	uint64_t seed_max_len;
	_M_seed_idx = 0;
	check_read(fread(&seed_max_len, sizeof(uint64_t), 1, _M_fin) == 1);
	if (seed_max_len > ref_max_len)
		utl_error("seqreader", 2,
			"reference sequence is built with '--max-len=%lu', but this alignement need %lu\n",
			ref_max_len, seed_max_len);
	check_read(fread(&_M_seed_cnt, sizeof(uint64_t), 1, _M_fin) == 1);
	_M_ss_size = fread(_M_ss_buf, sizeof(sequence), _M_read_size(UNIT), _M_fin);
	_M_seed_idx += _M_ss_size;
	_M_ss_head = 0;
	return seed_max_len;
}
void seq_reader::align(seq_reader *rd_qry, seq_reader *rd_ref, unsigned char mis_thr,
	uint64_t qry_max_len, bool bisulfite, unsigned char gap
) {
	mismatch_info *stk_mismatch, *tmp_mismatch;
	check_alloc(stk_mismatch = (mismatch_info*) malloc(sizeof(mismatch_info) *
		qry_max_len * (mis_thr + 1)));
	check_alloc(tmp_mismatch = (mismatch_info*) malloc(sizeof(mismatch_info) *
		qry_max_len * (mis_thr + 1)));
	while (rd_qry->_M_has_next() && rd_ref->_M_has_next()) {
		sequence *s_qry = rd_qry->_M_get_ss(0);
		if (rd_qry->_M_fasta->ignorable(s_qry->get_idx())) {
			rd_qry->_M_slide(1);
			continue;
		}
		int diff = _M_qry_vs_ref(s_qry, rd_ref->_M_get_ss(0),
			rd_qry->_M_fasta->get_seq_buffer(),
			rd_ref->_M_fasta->get_seq_buffer(),
			bisulfite);
		if (diff < 0) 
			rd_qry->_M_slide(1);
		else if (diff > 0) 
			rd_ref->_M_slide(1);
		else {
			rd_qry->_M_extend(-1, bisulfite);
			rd_ref->_M_extend(s_qry->get_len(), bisulfite);
			for (uint64_t q = 0; q<rd_qry->_M_get_window(); q++) {
				sequence* ss_qry = rd_qry->_M_get_ss(q);
				if (rd_qry->_M_fasta->ignorable(ss_qry->get_idx())) continue;
				uint64_t unit_len = 0;
				for (uint64_t i=ss_qry->get_loc(); rd_qry->_M_fasta->get_seq_buffer()[i]; i++) unit_len++;
				for (uint64_t i=ss_qry->get_loc() - 1; rd_qry->_M_fasta->get_seq_buffer()[i]; i--) unit_len++;
				unit_len /= mis_thr + 1;
				for (uint64_t r = 0; r<rd_ref->_M_get_window(); r++) {
					sequence* ss_ref = rd_ref->_M_get_ss(r);
					_M_try_match(ss_qry, ss_ref, rd_qry->_M_fasta, rd_ref->_M_fasta, 
						unit_len, mis_thr, stk_mismatch, tmp_mismatch, bisulfite, gap);
					if (rd_qry->_M_fasta->ignorable(ss_qry->get_idx())) break;
				}
			}
			if (rd_qry->_M_ss_tail < rd_qry->_M_ss_size 
					&& s_qry->get_len() == rd_qry->_M_get_ss(rd_qry->_M_get_window())->get_len())
				rd_ref->_M_slide(rd_ref->_M_ss_tail - rd_ref->_M_ss_head);
			rd_qry->_M_slide(rd_qry->_M_ss_tail - rd_qry->_M_ss_head);
		}
	}
	free(stk_mismatch);
	free(tmp_mismatch);
}
bool seq_reader::_M_has_next() {
	return _M_ss_head < _M_ss_size;
}
sequence* seq_reader::_M_get_ss(uint64_t idx) {
	return _M_ss_buf + _M_ss_head + idx;
}
int seq_reader::_M_qry_vs_ref (sequence *s_qry, sequence *s_ref, char *c_qry,
	char *c_ref, bool bisulfite
) {
	short i_qry = s_qry->get_len();
	short i_ref = s_ref->get_len();
	short i_cmp = i_qry <= i_ref ? i_qry : i_ref;
	c_qry += s_qry->get_loc();
	c_ref += s_ref->get_loc();
	if (bisulfite) {
		for (short i=0; i<i_cmp; i++) {
			char cq = c_qry[i];
			char cr = c_ref[i];
			if (cq == 'C') cq = 'T';
			if (cr == 'C') cr = 'T';
			char diff = cq - cr;
			if (diff) return diff;
		}
	} else {
		int diff = memcmp(c_qry, c_ref, i_cmp);
		if (diff) return diff;
/*
		for (short i=0; i<i_cmp; i++) {
			char diff = c_qry[i] - c_ref[i];
			if (diff) return diff;
		}
*/
	}
	return i_qry - i_cmp;
}
void seq_reader::_M_slide(uint64_t window) {
	_M_ss_head += window;
	if (_M_ss_head == _M_ss_size) {
		_M_ss_size = fread(_M_ss_buf, sizeof(sequence), _M_read_size(UNIT), _M_fin);
		_M_seed_idx += _M_ss_size;
		_M_ss_head = 0;
	}
}
void seq_reader::_M_extend(short len, bool bisulfite) {
//if len == -1, extend query sequence
//other, extend reference sequence
	_M_ss_tail = _M_ss_head + 1;
	for (;;) {
		if (_M_ss_tail == _M_ss_size) {
			if (_M_ss_head) {
				memmove(_M_ss_buf, _M_ss_buf + _M_ss_head, (_M_ss_tail - _M_ss_head) * sizeof(sequence));
				_M_ss_tail -= _M_ss_head;
				_M_ss_size -= _M_ss_head;
				_M_ss_head = 0;
				uint64_t tmp_read = fread(_M_ss_buf + _M_ss_tail, sizeof(sequence),
					_M_read_size(UNIT-_M_ss_tail), _M_fin);
				_M_ss_size += tmp_read;
				_M_seed_idx += tmp_read;
				if (_M_ss_size == _M_ss_tail) return;
			} else if (_M_ss_size == UNIT) {
				fprintf(stderr, "Warning: because sequence '");
				for (short i=0; i<_M_ss_buf[_M_ss_head].get_len(); i++)
					fputc(_M_fasta->get_seq_buffer()[_M_ss_buf[_M_ss_head].get_loc() + i], stderr);
				fprintf(stderr, "' repeated too many times, its alignment could be incomplete\n");
				return;
			}
		}
		if (len == -1) {
			if (_M_ss_buf[_M_ss_head].get_len() != _M_ss_buf[_M_ss_tail].get_len()) 
				break;
		} else if (_M_ss_buf[_M_ss_tail].get_len() < len)
			break;
		for (short i=0; i< (len == -1 ? _M_ss_buf[_M_ss_head].get_len() : len ); i++) {
			char diff = _M_fasta->get_seq_buffer()[_M_ss_buf[_M_ss_head].get_loc() + i] -
					_M_fasta->get_seq_buffer()[_M_ss_buf[_M_ss_tail].get_loc() + i];
			if (diff && (!bisulfite || (diff != 17 && diff != -17)))
				return;
		}	
		_M_ss_tail++;
	}
}
/*
void seq_reader::DEBUG(char *seq) {
	for (uint64_t i=_M_ss_head; i<_M_ss_tail; i++) {
		fprintf(stderr, "%lu:", i);
		_M_ss_buf[i].DEBUG(seq);
	}
}
*/
uint64_t seq_reader::_M_get_window() {
	return _M_ss_tail - _M_ss_head;
}
/*
step(1/-1)		forward or backward
ref_ending		store the ending location of reference sequence
qry_ending		store the ending location of qury sequence
seq_qry			query sequence 
seq_ref			reference sequence 
l_qry			query sequence location
l_ref			reference sequence location
mis_thr			mismatch threhshold
bisulfite		bisulfite conversion
gap			maximum allowed gap width
shift_qry		shift one sequence only, <0 means shift reference sequence, >0 means shift query sequence
stk_mismatch_start	mismatch information start location
stk_mis_num		final mismatch number
stk_mismatch_idx	final mismatch location
stk_mismatch		final mismatch information array
tmp_mis_num		temp mismatch number
tmp_mismatch_idx	temp mismatch location
tmp_mismatch		temp mismatch information array
unit_match		if a match is found in this seed, useful only when going backwards. to eliminate duplicated output
unit_idx		index in this seed
unit_len		seed length
*/
bool seq_reader::_M_gap_match(
	char step, uint64_t *ref_ending, uint64_t *qry_ending, char *seq_qry, char *seq_ref,
	uint64_t l_qry, uint64_t l_ref, unsigned char mis_thr, bool bisulfite, unsigned char gap,
	int shift_qry, uint64_t stk_mismatch_start, unsigned char *stk_mis_num,
	uint64_t *stk_mismatch_idx, mismatch_info *stk_mismatch, unsigned char tmp_mis_num,
	uint64_t tmp_mismatch_idx, mismatch_info *tmp_mismatch, bool unit_match,
	uint64_t unit_idx, uint64_t unit_len
) {
	if (step == -1 && unit_idx == unit_len) {
		if (unit_match) return false;
		unit_idx = 0;
		unit_match = true;
	}
	if (!seq_qry[l_qry]) {
		if (*stk_mis_num > tmp_mis_num || *stk_mis_num == UCHAR_MAX) {
			*ref_ending = l_ref-step;
			*qry_ending = l_qry-step;
			*stk_mis_num = tmp_mis_num;
			*stk_mismatch_idx = tmp_mismatch_idx;
			memcpy(stk_mismatch+stk_mismatch_start, tmp_mismatch+stk_mismatch_start,
				sizeof(mismatch_info)*(tmp_mismatch_idx-stk_mismatch_start));
		}
		return true;
	}
	if (bisulfite && seq_qry[l_qry] == seq_ref[l_ref] + 17) {
		tmp_mismatch[tmp_mismatch_idx].type = 'b';
		tmp_mismatch[tmp_mismatch_idx].loc_qry = l_qry;
		bool rtn = _M_gap_match(step, ref_ending, qry_ending, seq_qry, seq_ref,
			l_qry+step, l_ref+step, mis_thr, bisulfite, gap, 0, stk_mismatch_start,
			stk_mis_num, stk_mismatch_idx, stk_mismatch, tmp_mis_num, tmp_mismatch_idx+1,
			tmp_mismatch, unit_match, unit_idx+1, unit_len);
		if (!rtn) return false;
	} else if (seq_qry[l_qry] == seq_ref[l_ref] && seq_qry[l_qry] != 'N') {
		bool rtn = _M_gap_match(step, ref_ending, qry_ending, seq_qry, seq_ref,
			l_qry+step, l_ref+step, mis_thr, bisulfite, gap, 0, stk_mismatch_start,
			stk_mis_num, stk_mismatch_idx, stk_mismatch, tmp_mis_num, tmp_mismatch_idx,
			tmp_mismatch, unit_match, unit_idx+1, unit_len);
		if (!rtn) return false;
	} else if (tmp_mis_num < mis_thr && seq_ref[l_ref]) {
		tmp_mismatch[tmp_mismatch_idx].type = 'm';
		tmp_mismatch[tmp_mismatch_idx].loc_qry = l_qry;
		tmp_mismatch[tmp_mismatch_idx].loc_supp1 = l_ref;
		_M_gap_match(step, ref_ending, qry_ending, seq_qry, seq_ref, l_qry+step,
			l_ref+step, mis_thr, bisulfite, gap, 0, stk_mismatch_start, stk_mis_num,
			stk_mismatch_idx, stk_mismatch, tmp_mis_num+1, tmp_mismatch_idx+1,
			tmp_mismatch, false, unit_idx+1, unit_len);
	}
	if (shift_qry <= 0 && seq_ref[l_ref]) {
		if (shift_qry == 0 || shift_qry+gap == 0) {
			if (tmp_mis_num < mis_thr) {
				tmp_mismatch[tmp_mismatch_idx].type = 'i';
				tmp_mismatch[tmp_mismatch_idx].loc_qry = l_qry;
				tmp_mismatch[tmp_mismatch_idx].loc_supp1 = l_ref;
				tmp_mismatch[tmp_mismatch_idx].loc_supp2 = l_ref+step;
				_M_gap_match(step, ref_ending, qry_ending, seq_qry, seq_ref, l_qry,
					l_ref+step, mis_thr, bisulfite, gap, -1, stk_mismatch_start,
					stk_mis_num, stk_mismatch_idx, stk_mismatch, tmp_mis_num+1,
					tmp_mismatch_idx+1, tmp_mismatch, unit_idx == 0, unit_idx, unit_len);
			}
		} else {
			tmp_mismatch[tmp_mismatch_idx-1].loc_supp2 = l_ref+step;
			_M_gap_match(step, ref_ending, qry_ending, seq_qry, seq_ref, l_qry,
				l_ref+step, mis_thr, bisulfite, gap, shift_qry-1, stk_mismatch_start,
				stk_mis_num, stk_mismatch_idx, stk_mismatch, tmp_mis_num,
				tmp_mismatch_idx, tmp_mismatch, unit_match, unit_idx, unit_len);
		}
	}
	if (shift_qry >= 0) {
		if (shift_qry == 0 || shift_qry == gap) {
			if (tmp_mis_num < mis_thr) {
				tmp_mismatch[tmp_mismatch_idx].type = 'd';
				tmp_mismatch[tmp_mismatch_idx].loc_qry = l_qry;
				tmp_mismatch[tmp_mismatch_idx].loc_supp1 = l_qry+step;
				_M_gap_match(step, ref_ending, qry_ending, seq_qry, seq_ref, l_qry+step,
					l_ref, mis_thr, bisulfite, gap, 1, stk_mismatch_start, stk_mis_num,
					stk_mismatch_idx, stk_mismatch, tmp_mis_num+1, tmp_mismatch_idx+1,
					tmp_mismatch, false, unit_idx+1, unit_len);
			}
		} else {
			tmp_mismatch[tmp_mismatch_idx-1].loc_supp1 = l_qry+step;
			_M_gap_match(step, ref_ending, qry_ending, seq_qry, seq_ref, l_qry+step,
				l_ref, mis_thr, bisulfite, gap, shift_qry+1, stk_mismatch_start,
				stk_mis_num, stk_mismatch_idx, stk_mismatch, tmp_mis_num,
				tmp_mismatch_idx, tmp_mismatch, false, unit_idx+1, unit_len);
		}
	}
	return true;
}
void seq_reader::_M_try_match(sequence *ss_qry, sequence *ss_ref, fasta_reader *fa_qry,
	fasta_reader *fa_ref, uint64_t unit_len, unsigned char mis_thr,
	mismatch_info *stk_mismatch, mismatch_info *tmp_mismatch, bool bisulfite,
	unsigned char gap
) {
	unsigned char mis_num = 0;
	char *seq_qry = fa_qry->get_seq_buffer();
	char *seq_ref = fa_ref->get_seq_buffer();
	uint64_t shift_ref = ss_ref->get_loc() - ss_qry->get_loc();
	uint64_t stk_mismatch_mid = 0;
	if (bisulfite) {
		for (uint64_t idx_qry=ss_qry->get_loc();
			idx_qry < ss_qry->get_loc() + ss_qry->get_len();
			idx_qry++
		) {
			char c_ref = seq_ref[idx_qry + shift_ref];
			char c_qry = seq_qry[idx_qry];
			if (c_qry == c_ref + 17) {
				stk_mismatch[stk_mismatch_mid].type = 'b';
				stk_mismatch[stk_mismatch_mid++].loc_qry = idx_qry;
			} else if (c_ref != c_qry)
//				if (++mis_num > mis_thr) return;
				return;
		}
	}
	uint64_t ref_start, ref_end, qry_start, stk_mismatch_prev, stk_mismatch_next;
	if (!gap) {
		bool unit_match = true;
		uint64_t unit_idx = 0;
		stk_mismatch_prev = stk_mismatch_mid;
		uint64_t idx_qry = ss_qry->get_loc() - 1;
		uint64_t idx_ref = ss_ref->get_loc() - 1;
		while (seq_qry[idx_qry]) {
			if (bisulfite && seq_qry[idx_qry] == seq_ref[idx_ref] + 17) {
				stk_mismatch[stk_mismatch_prev].type = 'b';
				stk_mismatch[stk_mismatch_prev++].loc_qry = idx_qry;
				idx_ref--;
			} else if (seq_qry[idx_qry] != seq_ref[idx_ref] || seq_qry[idx_qry] == 'N') {
				if (++mis_num > mis_thr) return;
				unit_match = false;
				stk_mismatch[stk_mismatch_prev].type = 'm';
				stk_mismatch[stk_mismatch_prev].loc_qry = idx_qry;
				stk_mismatch[stk_mismatch_prev++].loc_supp1 = idx_ref;
				if (seq_ref[idx_ref]) idx_ref--;
			} else
				idx_ref--;
			if (++unit_idx == unit_len) {
				if (unit_match) return;
				unit_match = true;
				unit_idx = 0;
			}
			idx_qry --;
		}
		qry_start = idx_qry + 1;
		ref_start = idx_ref + 1;
		idx_qry = ss_qry->get_loc() + ss_qry->get_len(); 
		idx_ref = ss_ref->get_loc() + ss_qry->get_len(); 
		stk_mismatch_next = stk_mismatch_prev;
		while (seq_qry[idx_qry]) {
			if (bisulfite && seq_qry[idx_qry] == seq_ref[idx_ref] + 17) {
				stk_mismatch[stk_mismatch_next].type = 'b';
				stk_mismatch[stk_mismatch_next++].loc_qry = idx_qry;
				idx_ref ++;
			} else if (seq_qry[idx_qry] != seq_ref[idx_ref] || seq_qry[idx_qry] == 'N') {
				if (++mis_num > mis_thr) return;
				stk_mismatch[stk_mismatch_next].type = 'm';
				stk_mismatch[stk_mismatch_next].loc_qry = idx_qry;
				stk_mismatch[stk_mismatch_next++].loc_supp1 = idx_ref;
				if (seq_ref[idx_ref]) idx_ref++;
			} else
				idx_ref++;
			idx_qry ++;
		}
		ref_end = idx_ref - 1;	
	} else {
		stk_mismatch_prev = stk_mismatch_mid;
		unsigned char stk_mis_num = UCHAR_MAX;
		uint64_t qry_end;
		_M_gap_match(-1, &ref_start, &qry_start, seq_qry, seq_ref, ss_qry->get_loc()-1,
			ss_ref->get_loc()-1, mis_thr, bisulfite, gap, 0, stk_mismatch_mid,
			&stk_mis_num, &stk_mismatch_prev, stk_mismatch, mis_num, stk_mismatch_mid,
			tmp_mismatch, true, 0, unit_len);
		if (stk_mis_num == UCHAR_MAX) return;
		mis_num = stk_mis_num;
		stk_mismatch_next = stk_mismatch_prev;
		stk_mis_num = UCHAR_MAX;
		_M_gap_match(1, &ref_end, &qry_end, seq_qry, seq_ref,
			ss_qry->get_loc()+ss_qry->get_len(),
			ss_ref->get_loc()+ss_qry->get_len(), mis_thr, bisulfite, gap, 0,
			stk_mismatch_prev, &stk_mis_num, &stk_mismatch_next, stk_mismatch,
			mis_num, stk_mismatch_prev, tmp_mismatch, true, 0, unit_len);
		if (stk_mis_num == UCHAR_MAX) return;
	}
//	if (_G_limit && mis_num != mis_thr) return;
	unsigned int ref_idx = ss_ref->get_idx();
	int64_t ref_loc = fa_ref->get_ref_loc(ref_idx, ref_start, ref_end);
	_M_mismatch_idx = *_M_mismatch_buffer = 0;
	char delim = '\t';
	int score = 0;
	unsigned char *score_qry = fa_qry->get_score_buffer();
	for (uint64_t i=stk_mismatch_prev; i>stk_mismatch_mid;) 
		delim = _M_mismatch_print(true, --i, delim, seq_qry, seq_ref, shift_ref,
			stk_mismatch, qry_start, score_qry, &score);
	for (uint64_t i=0; i<stk_mismatch_mid; i++)
		delim = _M_mismatch_print(false, i, delim, seq_qry, seq_ref, shift_ref,
			stk_mismatch, qry_start, score_qry, &score);
	for (uint64_t i=stk_mismatch_prev; i<stk_mismatch_next; i++) 
		delim = _M_mismatch_print(false, i, delim, seq_qry, seq_ref, shift_ref,
			stk_mismatch, qry_start, score_qry, &score);
	fa_qry->add_match(ss_qry->get_idx(), ref_idx, ref_loc, score, _M_mismatch_buffer, mis_thr);
}
char seq_reader::_M_mismatch_print(bool backwards, uint64_t idx, char delim,
	char *seq_qry, char *seq_ref, uint64_t shift_ref,
	mismatch_info *stk_mismatch, uint64_t qry_start, unsigned char *score_qry,
	int *score
) {
	_M_mismatch_buffer[_M_mismatch_idx++] = delim;
	char type = stk_mismatch[idx].type;
	if (type == 'm') {
		uint64_t loc_qry = stk_mismatch[idx].loc_qry;
		_M_mismatch_idx += sprintf(_M_mismatch_buffer + _M_mismatch_idx, "%ju", loc_qry - qry_start + 1);
		_M_mismatch_buffer[_M_mismatch_idx++] = ':';
		_M_mismatch_buffer[_M_mismatch_idx++] = type;
		uint64_t loc_ref = stk_mismatch[idx].loc_supp1;
		_M_mismatch_buffer[_M_mismatch_idx++] = seq_qry[loc_qry];
		if (seq_ref[loc_ref])
			_M_mismatch_buffer[_M_mismatch_idx++] = seq_ref[loc_ref];
		else
			_M_mismatch_buffer[_M_mismatch_idx++] = '0';
		if (score_qry) {
			*score += score_qry[loc_qry];
			_M_mismatch_buffer[_M_mismatch_idx++] = ':';
			_M_mismatch_idx += sprintf(_M_mismatch_buffer + _M_mismatch_idx, "%u", score_qry[loc_qry]);
		}
	} else if (type == 'i') {
		uint64_t loc_qry = stk_mismatch[idx].loc_qry;
		if (backwards) {
			_M_mismatch_idx += sprintf(_M_mismatch_buffer + _M_mismatch_idx, "%ju",
				loc_qry - qry_start + 1);
			_M_mismatch_buffer[_M_mismatch_idx++] = ':';
			_M_mismatch_buffer[_M_mismatch_idx++] = type;
			for (uint64_t i=stk_mismatch[idx].loc_supp2+1; i<=stk_mismatch[idx].loc_supp1; i++)
				_M_mismatch_buffer[_M_mismatch_idx++] = seq_ref[i];
		} else {
			_M_mismatch_idx += sprintf(_M_mismatch_buffer + _M_mismatch_idx, "%ju", loc_qry - qry_start);
			_M_mismatch_buffer[_M_mismatch_idx++] = ':';
			_M_mismatch_buffer[_M_mismatch_idx++] = type;
			for (uint64_t i=stk_mismatch[idx].loc_supp1; i<stk_mismatch[idx].loc_supp2; i++)
				_M_mismatch_buffer[_M_mismatch_idx++] = seq_ref[i];
		}
	} else if (type == 'd') {
		if (backwards) {
			_M_mismatch_idx += sprintf(_M_mismatch_buffer + _M_mismatch_idx, "%ju",
				stk_mismatch[idx].loc_supp1 - qry_start + 2);
			_M_mismatch_buffer[_M_mismatch_idx++] = ':';
			_M_mismatch_buffer[_M_mismatch_idx++] = type;
			for (uint64_t i=stk_mismatch[idx].loc_supp1+1; i<=stk_mismatch[idx].loc_qry; i++)
				_M_mismatch_buffer[_M_mismatch_idx++] = seq_qry[i];
		} else {
			uint64_t loc_qry = stk_mismatch[idx].loc_qry;
			_M_mismatch_idx += sprintf(_M_mismatch_buffer + _M_mismatch_idx,
				"%ju", loc_qry - qry_start + 1);
			_M_mismatch_buffer[_M_mismatch_idx++] = ':';
			_M_mismatch_buffer[_M_mismatch_idx++] = type;
			for (uint64_t i=stk_mismatch[idx].loc_qry; i<stk_mismatch[idx].loc_supp1; i++)
				_M_mismatch_buffer[_M_mismatch_idx++] = seq_qry[i];
		}
	} else if (type == 'b') {
		uint64_t loc_qry = stk_mismatch[idx].loc_qry;
		_M_mismatch_idx += sprintf(_M_mismatch_buffer + _M_mismatch_idx, "%ju", loc_qry - qry_start + 1);
		_M_mismatch_buffer[_M_mismatch_idx++] = ':';
		_M_mismatch_buffer[_M_mismatch_idx++] = type;
	} else
		utl_error("seqreader", 3, "wrong mismatch type, contact developer\n");
	_M_mismatch_buffer[_M_mismatch_idx] = 0;
	return ',';
}	
uint64_t seq_reader::_M_read_size(uint64_t len) {
	uint64_t rtn = _M_seed_cnt - _M_seed_idx;
	if (len < rtn) rtn = len;
	return rtn;
}
