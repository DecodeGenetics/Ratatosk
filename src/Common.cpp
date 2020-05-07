#include "Common.hpp"

// returns maximal entropy score for random sequences which is log2|A| where A is the alphabet
// returns something close to 0 for highly repetitive sequences
double getEntropy(const char* s, const size_t len){

	double counts[4] = {0, 0, 0, 0};
	double logs[4];

	for (size_t i = 0; i != len; ++i){

        const char c = s[i] & 0xDF;

		if (isDNA(c)){

			const size_t x = (c & 4) >> 1;

			++counts[x + ((x ^ (c & 2)) >> 1)];
		}
	}

	counts[0] /= len;
	counts[1] /= len;
	counts[2] /= len;
	counts[3] /= len;

	logs[0] = (counts[0] == 0) ? 0 : log2(counts[0]);
	logs[1] = (counts[1] == 0) ? 0 : log2(counts[1]);
	logs[2] = (counts[2] == 0) ? 0 : log2(counts[2]);
	logs[3] = (counts[3] == 0) ? 0 : log2(counts[3]);

	return 0 - ((counts[0] * logs[0]) + (counts[1] * logs[1]) + (counts[2] * logs[2]) + (counts[3] * logs[3]));
}

size_t getMaxPaths(const double seq_entropy, const size_t max_len_path, const size_t k){

	const double entropy_factor = 3. - seq_entropy;
	const size_t v = rndup(max_len_path);

	return static_cast<double>(v * (__builtin_ctzll(v) + k)) * entropy_factor;
}

size_t getMaxBranch(const double seq_entropy, const size_t max_len_path, const size_t k){

	const double entropy_factor = 3. - seq_entropy;
	const double base_nb_nodes = static_cast<double>(__builtin_ctzll(rndup(max_len_path - k + 1)));

	return (base_nb_nodes + 1) * entropy_factor;
}

bool hasEnoughSharedPairID(const PairID& a, const PairID& b, const size_t min_shared_ids) {

	size_t nb_shared = 0;

	PairID::const_iterator a_it_s = a.begin();
	PairID::const_iterator b_it_s = b.begin();

	PairID::const_iterator a_it_e = a.end();
	PairID::const_iterator b_it_e = b.end();

	const size_t a_card = a.cardinality();
	const size_t b_card = b.cardinality();

	if (a_card >= 10 * b_card){

		while ((b_it_s != b_it_e) && (nb_shared < min_shared_ids)){

			nb_shared += static_cast<size_t>(a.contains(*b_it_s));
			++b_it_s;
		}
	}
	else if (b_card >= 10 * a_card){

		while ((a_it_s != a_it_e) && (nb_shared < min_shared_ids)){

			nb_shared += static_cast<size_t>(b.contains(*a_it_s));
			++a_it_s;
		}
	}
	else {

		while ((a_it_s != a_it_e) && (b_it_s != b_it_e) && (nb_shared < min_shared_ids)){

			const uint32_t val_a = *a_it_s;
			const uint32_t val_b = *b_it_s;

			if (val_a == val_b){

				++nb_shared;
				++a_it_s;
				++b_it_s;
			}
			else if (val_a < val_b) ++a_it_s;
			else ++b_it_s;
		}
	}

	return (nb_shared == min_shared_ids);
}

bool hasEnoughSharedPairID(const TinyBloomFilter<uint32_t>& tbf_a, const PairID& a, const PairID& b, const size_t min_shared_ids) {

	size_t nb_shared = 0;

	PairID::const_iterator a_it_s = a.begin();
	PairID::const_iterator b_it_s = b.begin();

	PairID::const_iterator a_it_e = a.end();
	PairID::const_iterator b_it_e = b.end();

	const size_t a_card = a.cardinality();
	const size_t b_card = b.cardinality();

	if (a_card >= 10 * b_card){

		while ((b_it_s != b_it_e) && (nb_shared < min_shared_ids)){

			const uint32_t val_b = *b_it_s;

			nb_shared += static_cast<size_t>(tbf_a.query(val_b) && a.contains(val_b));
			++b_it_s;
		}
	}
	else if (b_card >= 10 * a_card){

		while ((a_it_s != a_it_e) && (nb_shared < min_shared_ids)){

			nb_shared += static_cast<size_t>(b.contains(*a_it_s));
			++a_it_s;
		}
	}
	else {

		while ((a_it_s != a_it_e) && (b_it_s != b_it_e) && (nb_shared < min_shared_ids)){

			const uint32_t val_a = *a_it_s;
			const uint32_t val_b = *b_it_s;

			if (val_a == val_b){

				++nb_shared;
				++a_it_s;
				++b_it_s;
			}
			else if (val_a < val_b) ++a_it_s;
			else ++b_it_s;
		}
	}

	return (nb_shared == min_shared_ids);
}

size_t getNumberSharedPairID(const PairID& a, const PairID& b) {

	if (a.isEmpty() || b.isEmpty()) return 0;

	size_t nb_shared = 0;

	PairID::const_iterator a_it_s = a.begin();
	PairID::const_iterator b_it_s = b.begin();

	PairID::const_iterator a_it_e = a.end();
	PairID::const_iterator b_it_e = b.end();

	while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

		if (*a_it_s == *b_it_s){

			++nb_shared;
			++a_it_s;
			++b_it_s;
		}
		else if (*a_it_s < *b_it_s) ++a_it_s;
		else ++b_it_s;
	}

	return nb_shared;
}