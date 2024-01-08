#ifndef SRC_SEQUENCE_H_
#define SRC_SEQUENCE_H_
#include <stdint.h>
class sequence {
public:
	sequence(uint64_t loc, short len, unsigned int idx);
	uint64_t get_loc() const;
	short get_len() const;
	unsigned int get_idx();
	char getchar(char *seq, short idx);
	bool operator< (const sequence& ano) const;
	static char *buf_4_less_than;
	static bool bisulfite;

	void DEBUG(char *seq);
private:
	uint64_t _M_storage;
	unsigned int _M_idx;
} __attribute__((packed));
#endif
