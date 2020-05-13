#ifndef RATATOSK_COMMON_HPP
#define RATATOSK_COMMON_HPP

#include <iostream>

#include <bifrost/CompactedDBG.hpp>

#include "PairID.hpp"
#include "TinyBloomFilter.hpp"

#define RATATOSK_VERSION "0.1"

struct Correct_Opt : CDBG_Build_opt {

	vector<string> filenames_unmapped_short_in;
	vector<string> filenames_long_in;
	vector<string> filenames_helper_long_in;

	string filename_long_out;
	string filename_unmapped_short_graph_in;

	size_t min_qv;

	int out_qual;
    int trim_qual;

    size_t large_k;

    size_t min_cov_vertices;
    size_t min_cov_edges;
    size_t max_cov_vertices;

    size_t nb_partitions;
    size_t max_time;

    double weak_region_len_factor;

    bool correct;
    bool index;

	Correct_Opt() : out_qual(0), trim_qual(0), min_qv(6), large_k(95), nb_partitions(1000), max_time(20),
					min_cov_vertices(2), min_cov_edges(2), max_cov_vertices(512), weak_region_len_factor(1.25),
					correct(false), index(false) {

        k = 63;

        clipTips = true;
        deleteIsolated = true;
        useMercyKmers = false;
    }
};

// returns maximal entropy score for random sequences which is log2|A| where A is the alphabet
// returns something close to 0 for highly repetitive sequences
double getEntropy(const char* s, const size_t len);

size_t getMaxPaths(const double seq_entropy, const size_t max_len_path, const size_t k);
size_t getMaxBranch(const double seq_entropy, const size_t max_len_path, const size_t k);

bool hasEnoughSharedPairID(const PairID& a, const PairID& b, const size_t min_shared_ids);
bool hasEnoughSharedPairID(const TinyBloomFilter<uint32_t>& tbf_a, const PairID& a, const PairID& b, const size_t min_shared_ids);

size_t getNumberSharedPairID(const PairID& a, const PairID& b);

// This order of nucleotide in this array is important
static const char ambiguity_c[16] = {'.', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

inline void copyFile(const string& dest, const vector<string>& src){

    std::ofstream  out(dest, std::ios::binary | std::ios::app);

    for (const auto& s : src){

	    std::ifstream  in(s, std::ios::binary);

	    out << in.rdbuf();
    }
}

inline size_t countRecordsFASTX(const string& filename){

	const vector<string> v(1, filename);
	
	FileParser fp(v);

	size_t i = 0;
	size_t nb_rec = 0;

    string s;

	while (fp.read(s, i)) ++nb_rec;

	fp.close();

	return nb_rec;
}

inline size_t countRecordsFASTX(const vector<string>& filenames){
	
	FileParser fp(filenames);

	size_t i = 0;
	size_t nb_rec = 0;

    string s;

	while (fp.read(s, i)) ++nb_rec;

	fp.close();

	return nb_rec;
}

inline char getAmbiguity(const bool nuc_a, const bool nuc_c, const bool nuc_g, const bool nuc_t) {

	const uint8_t idx = static_cast<uint8_t>(nuc_a) | (static_cast<uint8_t>(nuc_c) << 1) | (static_cast<uint8_t>(nuc_g) << 2) | (static_cast<uint8_t>(nuc_t) << 3);

	return ambiguity_c[idx];
}

inline uint8_t getAmbiguityIndex(const char nuc_ambiguity){

	const char c = nuc_ambiguity & 0xDF;

	uint8_t idx = 0;

	idx += (static_cast<size_t>(!(c == ambiguity_c[0])) - 1) & 0x0;
	idx += (static_cast<size_t>(!(c == ambiguity_c[1])) - 1) & 0x1;
	idx += (static_cast<size_t>(!(c == ambiguity_c[2])) - 1) & 0x2;
	idx += (static_cast<size_t>(!(c == ambiguity_c[3])) - 1) & 0x3;
	idx += (static_cast<size_t>(!(c == ambiguity_c[4])) - 1) & 0x4;
	idx += (static_cast<size_t>(!(c == ambiguity_c[5])) - 1) & 0x5;
	idx += (static_cast<size_t>(!(c == ambiguity_c[6])) - 1) & 0x6;
	idx += (static_cast<size_t>(!(c == ambiguity_c[7])) - 1) & 0x7;
	idx += (static_cast<size_t>(!(c == ambiguity_c[8])) - 1) & 0x8;
	idx += (static_cast<size_t>(!(c == ambiguity_c[9])) - 1) & 0x9;
	idx += (static_cast<size_t>(!(c == ambiguity_c[10])) - 1) & 0xA;
	idx += (static_cast<size_t>(!(c == ambiguity_c[11])) - 1) & 0xB;
	idx += (static_cast<size_t>(!(c == ambiguity_c[12])) - 1) & 0xC;
	idx += (static_cast<size_t>(!(c == ambiguity_c[13])) - 1) & 0xD;
	idx += (static_cast<size_t>(!(c == ambiguity_c[14])) - 1) & 0xE;
	idx += (static_cast<size_t>(!(c == ambiguity_c[15])) - 1) & 0xF;

	return idx;
}

inline char getAmbiguityIndexRev(const uint8_t idx){

	return ambiguity_c[idx];
}

inline void getAmbiguityRev(const char nuc_ambiguity, bool& nuc_a, bool& nuc_c, bool& nuc_g, bool& nuc_t) {

	const uint8_t idx = getAmbiguityIndex(nuc_ambiguity);

	nuc_a = static_cast<bool>(idx & 0x1);
	nuc_c = static_cast<bool>(idx & 0x2);
	nuc_g = static_cast<bool>(idx & 0x4);
	nuc_t = static_cast<bool>(idx & 0x8);

	return;
}

inline void getStdQual(string& s) {

	for (auto& c : s){

		if (c < static_cast<char>(33)) c = static_cast<char>(33);
		if (c > static_cast<char>(73)) c = static_cast<char>(73);
	}
}

inline char getQual(const double score, const size_t qv_min = 0) {

	const char phred_base_std = static_cast<char>(33);
	const char phred_scale_std = static_cast<char>(40);

	const double qv_score = score * static_cast<double>(phred_scale_std - qv_min);

	return static_cast<char>(qv_score + phred_base_std + qv_min);
}

inline double getScore(const char c, const size_t qv_min = 0) {

	const char phred_base_std =  static_cast<char>(33);
	const char phred_scale_std =  static_cast<char>(40);

	const double qv_score = static_cast<double>(c - phred_base_std - qv_min);

	return (qv_score / static_cast<double>(phred_scale_std - qv_min));
}

#endif