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

	const size_t a_card = a.cardinality();
	const size_t b_card = b.cardinality();

	if ((a_card >= min_shared_ids) && (b_card >= min_shared_ids)) {

		const size_t log2_a = b_card * approximate_log2(a_card);
		const size_t log2_b = a_card * approximate_log2(b_card);

		const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

		if (min_a_b == (a_card + b_card)) {

			PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();
			PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

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
		else if (min_a_b == log2_a){

			PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

			while ((b_it_s != b_it_e) && (nb_shared < min_shared_ids)){

				nb_shared += static_cast<size_t>(a.contains(*b_it_s));
				++b_it_s;
			}
		}
		else {

			PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();

			while ((a_it_s != a_it_e) && (nb_shared < min_shared_ids)){

				nb_shared += static_cast<size_t>(b.contains(*a_it_s));
				++a_it_s;
			}
		}
	}

	return (nb_shared == min_shared_ids);
}

bool hasEnoughSharedPairID(const TinyBloomFilter<uint32_t>& tbf_a, const PairID& a, const PairID& b, const size_t min_shared_ids) {

	size_t nb_shared = 0;

	const size_t a_card = a.cardinality();
	const size_t b_card = b.cardinality();

	if ((a_card >= min_shared_ids) && (b_card >= min_shared_ids)) {

		const size_t log2_a = b_card * (approximate_log2(a_card) + tbf_a.getNumberHashFunctions());
		const size_t log2_b = a_card * approximate_log2(b_card);

		const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

		if (min_a_b == log2_a) {

			PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

			while ((b_it_s != b_it_e) && (nb_shared < min_shared_ids)){

				const uint32_t val_b = *b_it_s;

				nb_shared += static_cast<size_t>(tbf_a.query(val_b) && a.contains(val_b));
				++b_it_s;
			}
		}
		else return hasEnoughSharedPairID(a, b, min_shared_ids);
	}

	return (nb_shared == min_shared_ids);
}

size_t getNumberSharedPairID(const PairID& a, const PairID& b) {

	if (a.isEmpty() || b.isEmpty()) return 0;

	size_t nb_shared = 0;

	const size_t a_card = a.cardinality();
	const size_t b_card = b.cardinality();

	const size_t log2_a = b_card * approximate_log2(a_card);
	const size_t log2_b = a_card * approximate_log2(b_card);

	const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

	if (min_a_b == (a_card + b_card)) {

		PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();
		PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

		while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

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
	else if (min_a_b == log2_a){

		PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

		while (b_it_s != b_it_e){

			nb_shared += static_cast<size_t>(a.contains(*b_it_s));
			++b_it_s;
		}
	}
	else {

		PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();

		while (a_it_s != a_it_e){

			nb_shared += static_cast<size_t>(b.contains(*a_it_s));
			++a_it_s;
		}
	}

	return nb_shared;
}

size_t getNumberSharedPairID(const TinyBloomFilter<uint32_t>& tbf_a, const PairID& a, const PairID& b) {

	if (a.isEmpty() || b.isEmpty()) return 0;

	size_t nb_shared = 0;

	const size_t a_card = a.cardinality();
	const size_t b_card = b.cardinality();

	const size_t log2_a = b_card * (approximate_log2(a_card) + tbf_a.getNumberHashFunctions());
	const size_t log2_b = a_card * approximate_log2(b_card);

	const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

	if (min_a_b == log2_a) {

		PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

		while (b_it_s != b_it_e){

			const uint32_t val_b = *b_it_s;

			nb_shared += static_cast<size_t>(tbf_a.query(val_b) && a.contains(val_b));
			++b_it_s;
		}
	}
	else return getNumberSharedPairID(a, b);

	return nb_shared;
}

size_t approximate_log2(size_t v) {

	static const uint8_t tab64[64] = {
	    63,  0, 58,  1, 59, 47, 53,  2,
	    60, 39, 48, 27, 54, 33, 42,  3,
	    61, 51, 37, 40, 49, 18, 28, 20,
	    55, 30, 34, 11, 43, 14, 22,  4,
	    62, 57, 46, 52, 38, 26, 32, 41,
	    50, 36, 17, 19, 29, 10, 13, 21,
	    56, 45, 25, 31, 35, 16,  9, 12,
	    44, 24, 15,  8, 23,  7,  6,  5
	};

    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;

    return tab64[(static_cast<size_t>((v - (v >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}