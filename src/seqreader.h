#ifndef SRC_SEQREADER_H_
#define SRC_SEQREADER_H_
#include <stdio.h>
#include <stdint.h>
#include "fastareader.h"
#include "sequence.h"
struct mismatch_info {
	char type;
	uint64_t loc_qry, loc_supp1, loc_supp2;
};
class seq_reader {
public:
	seq_reader(fasta_reader *fasta, FILE *fin, const char *filename);
	~seq_reader();
	static void align(seq_reader *rd_qry, seq_reader *rd_ref, unsigned char mis_thr,
		uint64_t qry_max_len, bool bisulfite, unsigned char gap);
	void init_ref();
	void jump_to(int step);
	uint64_t get_max_seed_len(uint64_t ref_max_len);
private:
	sequence *_M_ss_buf;
	uint64_t _M_ss_size, _M_ss_head, _M_ss_tail, _M_seed_cnt, _M_seed_idx;
	int64_t _M_ftell;
	FILE *_M_fin;
	fasta_reader *_M_fasta;
	const char * _M_filename;
	static char * _M_mismatch_buffer;
	static int _M_mismatch_idx;

	bool _M_has_next();
	void _M_slide(uint64_t window);
	void _M_extend(short len, bool bisulfite);
	sequence* _M_get_ss(uint64_t idx);
	uint64_t _M_get_window();
	uint64_t _M_read_size(uint64_t len);
	static int _M_qry_vs_ref (sequence *s_qry, sequence *s_ref, char *c_qry,
		char *c_ref, bool bisulfite);
	static void _M_try_match(sequence *ss_qry, sequence *ss_ref, fasta_reader *fa_qry,
		fasta_reader *fa_ref, uint64_t unit_len, unsigned char mis_thr,
		mismatch_info *stk_mismatch, mismatch_info *tmp_mismatch, bool bisulfite,
		unsigned char gap);
	static char _M_mismatch_print(bool backwards, uint64_t idx, char delim,
		char *seq_qry, char *seq_ref, uint64_t shift_ref,
		mismatch_info *stk_mismatch, uint64_t qry_start, unsigned char *score_qry,
		int *score);
	static bool _M_gap_match(
		char step, uint64_t *ref_ending, uint64_t *qry_ending, char *seq_qry, char *seq_ref,
		uint64_t l_qry, uint64_t l_ref, unsigned char mis_thr, bool bisulfite, unsigned char gap,
		int shift_qry, uint64_t stk_mismatch_start, unsigned char *stk_mis_num,
		uint64_t *stk_mismatch_idx, mismatch_info *stk_mismatch, unsigned char tmp_mis_num,
		uint64_t tmp_mismatch_idx, mismatch_info *tmp_mismatch, bool unit_match,
		uint64_t unit_idx, uint64_t unit_len
	);
//	void DEBUG(char *seq);
};
#endif
