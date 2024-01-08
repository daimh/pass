#ifndef SRC_FASTAREADER_H_
#define SRC_FASTAREADER_H_
#include <stdint.h>
struct match_info {
	unsigned int ref_idx;
	int64_t ref_loc;
	int score;
	char* mismatch;
} __attribute__((packed));
class fasta_reader {
public:
	fasta_reader(char *fa_name, bool is_ref, char *id_from, char *id_to, FILE *fout);
	fasta_reader(bool is_ref, FILE *fin, int limit, unsigned char mis_thr,
		int quality_score_base, bool print_seq, bool print_all);
	~fasta_reader();
	char *get_seq_buffer();
	unsigned char *get_score_buffer();
	void add_match(unsigned int qry_idx, unsigned int ref_idx, uint64_t ref_loc,
		int score, char* mismatch, unsigned char mis_idx);
	void print_match();
	void print_nomatch();
	void print_ref_fasta_id(unsigned int ref_idx, int64_t ref_loc);
	int64_t get_ref_loc(unsigned int ref_idx, uint64_t ref_start, uint64_t ref_end);
	uint64_t get_size();
	bool ignorable(unsigned int qry_idx);
//	void DEBUG();
	static fasta_reader* _G_fa_ref;
private:
	bool _M_is_ref;
	uint64_t _M_id_num;
	bool _M_has_score;
	char *_M_id_buf;
	uint64_t *_M_id_loc;
	char *_M_seq_buf;
	unsigned char *_M_score_buf;
	bool _M_print_seq, _M_print_all;
	uint64_t *_M_seq_loc;
	uint64_t _M_size;
	int _M_limit;
	match_info** _M_match_buf;
	unsigned char *_M_match_cnt;
	unsigned char _M_mis_thr;
	uint64_t *_M_counter_other;
	
	void _M_finish_seq_read(uint64_t *seq_buf_tail, char *fa_name, bool miss_seq,
		char *line);
	uint64_t _M_get_limit();
	void _M_finish_new(uint64_t id_buf_tail, uint64_t seq_buf_tail);
	void _M_put_id(char *line, char *id_from, char *id_to,
		uint64_t seq_buf_tail, uint64_t *id_buf_tail, uint64_t *id_buf_size,
		uint64_t *id_size);
	void _M_put_seq(char *line, int64_t line_read, uint64_t *seq_buf_tail,
		uint64_t *seq_buf_size);
};
#endif
