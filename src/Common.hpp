#ifndef RATATOSK_COMMON_HPP
#define RATATOSK_COMMON_HPP

#include <iostream>
#include <random>

#include "CompactedDBG.hpp"

#include "edlib.h"
#include "PairID.hpp"
#include "TinyBloomFilter.hpp"
#include "SharedPairID.hpp"

#define RATATOSK_VERSION "0.9.0"

struct Correct_Opt : CDBG_Build_opt {

	vector<string> filenames_long_in; // Long reads to correct
	vector<string> filenames_long_raw; // Raw long reads to correct. This is the same as "filenames_long_in" during pass 1.

	vector<string> filenames_helper_long_in; // Accurate long reads helping with coloring on the 2nd round
	vector<string> filenames_short_all; // Unmapped short reads

	vector<string> filenames_long_phase; // Phasing files long reads
	vector<string> filenames_short_phase; // Phasing files short reads

	string filename_long_out; // Output filename prefix for long reads

	string fn_graph_in;
	string fn_index_in;

	bool correct;
	bool index;

	bool pass1_only;
	bool pass2_only;

	int out_qual;
    int trim_qual;
    int max_qual;

    size_t small_k;
    size_t insert_sz;

    size_t min_cov_vertices;
    size_t max_cov_vertices;

    size_t max_km_cov;
    double top_km_cov_ratio;

    size_t min_nb_km_unmapped;

    size_t nb_correction_rounds;

    size_t nb_partitions;
    size_t min_bases_partition;

    size_t max_len_weak_region1;
    size_t max_len_weak_region2;

    size_t min_len_2nd_pass;

    size_t buffer_sz_read2disk;
    size_t buffer_sz;
    size_t h_seed;

    double weak_region_len_factor;
    double large_k_factor;
    double min_score;
    double min_color_sharing;
    double sampling_rate;

    double min_confidence_snp_corr;
    double min_confidence_2nd_pass;

    bool force_unres_snp_corr;
    bool force_no_snp_corr;
    bool force_io_order;
    bool force_no_graph_index;
    bool compress_out;

	Correct_Opt() {

		clear();
	}

	void clear() {

		filenames_long_in.clear();
		filenames_long_raw.clear();
		filenames_helper_long_in.clear();
		filenames_short_all.clear();
		filenames_long_phase.clear();
		filenames_short_phase.clear();

		filename_long_out.clear();

		fn_graph_in.clear();
		fn_index_in.clear();

		k = 63;

		clipTips = false;
		deleteIsolated = false;
		useMercyKmers = false;

		correct = false;
		index = false;

		pass1_only = false;
		pass2_only = false;

		max_qual = 40;
		out_qual = 1;
	    trim_qual = 0;

	    small_k = 31;
	    insert_sz = 500;

	    min_cov_vertices = 2;
	    max_cov_vertices = 128;

	    max_km_cov = 128;
	    top_km_cov_ratio = 0.001;

	    nb_correction_rounds = 1;

	    nb_partitions = 1000;
	    min_bases_partition = 100000;

	    max_len_weak_region1 = 1000;
	    max_len_weak_region2 = 5000;

	    min_len_2nd_pass = 3000;

	    buffer_sz_read2disk = 0x00000000ffffffffULL; // Read 4GB of sequence before writing to disk
	    //buffer_sz_read2disk = 0x1000000ULL; // Read 16 MB of sequence before writing to disk -> for testing only
	    buffer_sz = 1048576; // 1 MB reading buffers per thread
	    h_seed = 0;

	    weak_region_len_factor = 0.25;
	    large_k_factor = 1.5;
	    min_score = 0.0;
	    min_color_sharing = 0.5;
	    sampling_rate = 1.0;

	    min_confidence_snp_corr = 0.9;
	    min_confidence_2nd_pass = 0.0;

	    force_unres_snp_corr = false;
	    force_no_snp_corr = false;
	    force_io_order = false;
	    force_no_graph_index = false;
	    compress_out = false;

	    min_nb_km_unmapped = small_k;
	}
};

struct HashSharedPairIDptr {

    size_t operator()(const SharedPairID* ptr) const {

        return wyhash(&ptr, sizeof(const SharedPairID*), 0, _wyp);
    }
};

struct CustomHashUint64_t {

    size_t operator()(const uint64_t v) const {

        return wyhash(&v, sizeof(uint64_t), 0, _wyp);
    }
};

struct CustomHashSize_t {

    size_t operator()(const size_t v) const {

        return wyhash(&v, sizeof(size_t), 0, _wyp);
    }
};

struct CustomHashString {

    size_t operator()(const string& s) const {

        return wyhash(s.c_str(), s.length(), 0, _wyp);
    }
};

struct HapReads {

	unordered_map<uint64_t, uint64_t, CustomHashUint64_t> read2hap;

	unordered_map<string, uint64_t, CustomHashString> hapBlock2id;
	unordered_map<string, uint64_t, CustomHashString> hapType2id;

	vector<PairID> hap2phasedReads;
	PairID hap2unphasedReads;

    size_t block_id;
    size_t type_id;

    HapReads() {

    	clear();
    }

    void clear() {

    	read2hap.clear();

    	hapBlock2id.clear();
    	hapType2id.clear();

    	hap2phasedReads.clear();
    	hap2unphasedReads.clear();

    	block_id = 0;
    	type_id = 0;
    }
};

struct WeightsPairID {

	PairID noWeight_pids;
	PairID weighted_pids;
	PairID all_pids;

	double weight;
	double sum_pids_weights;

	WeightsPairID() {

		clear();
	};

	void clear() {

		noWeight_pids.clear();
		weighted_pids.clear();
		all_pids.clear();

		weight = 2.0;
		sum_pids_weights = 0.0;
	}
};

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

// !!! This order of the nucleotide in this array is important -> DO NOT CHANGE !!!
static const char ambiguity_c[16] = {'.', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

static const EdlibEqualityPair edlib_iupac_alpha[28] = {
															{'M', 'A'}, {'M', 'C'},
															{'R', 'A'}, {'R', 'G'},
															{'S', 'C'}, {'S', 'G'},
															{'V', 'A'}, {'V', 'C'}, {'V', 'G'},
															{'W', 'A'}, {'W', 'T'},
															{'Y', 'C'}, {'Y', 'T'},
															{'H', 'A'}, {'H', 'C'}, {'H', 'T'},
															{'K', 'G'}, {'K', 'T'},
															{'D', 'A'}, {'D', 'G'}, {'D', 'T'},
															{'B', 'C'}, {'B', 'G'}, {'B', 'T'},
															{'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
														};

static const size_t sz_edlib_iupac_alpha = 28;

// returns maximal entropy score for random sequences which is log2|A| where A is the alphabet
// returns something close to 0 for highly repetitive sequences
double getEntropy(const char* s, const size_t len);

size_t getMaxPaths(const double seq_entropy, const size_t max_len_path, const size_t k);
size_t getMaxBranch(const double seq_entropy, const size_t max_len_path, const size_t k);

size_t getNumberSharedPairID(const SharedPairID& a, const PairID& b, const size_t min_shared_ids);
size_t getNumberSharedPairID(const SharedPairID& a, const SharedPairID& b, const size_t min_shared_ids);

size_t getNumberSharedPairID(const SharedPairID& a, const PairID& b);
size_t getNumberSharedPairID(const SharedPairID& a, const SharedPairID& b);

PairID getSharedPairID(const SharedPairID& a, const SharedPairID& b, const size_t min_shared);
PairID getSharedPairID(const SharedPairID& a, const PairID& b, const size_t min_shared);
PairID getSharedPairID(const PairID& a, const PairID& b, const size_t min_shared);

PairID subsample(const SharedPairID& spid, const size_t nb_id_out);
PairID subsample(const PairID& spid, const size_t nb_id_out);

inline void toUpperCase(char* s, const size_t len) {

	for (char* s_tmp = s; s_tmp < (s + len); ++s_tmp) *s_tmp &= 0xDF;
}

inline void copyFile(const string& dest, const vector<string>& src){

    std::ofstream  out(dest, std::ios::binary | std::ios::app);

    for (const auto& s : src){

	    std::ifstream  in(s, std::ios::binary);

	    out << in.rdbuf();
    }
}

inline size_t countRecords(const vector<string>& v_fn, const bool unique, const uint64_t seed){

	size_t i = 0, n = 0;

    string s;

	FileParser fp(v_fn);

	if (v_fn.empty()) return 0;

	if (unique) { // Count approximate number of unique records

		unordered_set<uint64_t> sh;

		while (fp.read(s, i)){

			const uint64_t h = wyhash(fp.getNameString(), strlen(fp.getNameString()), seed, _wyp);

			n += static_cast<size_t>(sh.insert(h).second);
		}
	}
	else {

		while (fp.read(s, i)) ++n;
	}

	return n;
}

inline size_t countRecords(const string& fn, const bool unique, const uint64_t seed){

	const vector<string> v(1, fn);
	
	return countRecords(v, unique, seed);
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

inline void getStdQual(string& s, const size_t qv_max = 40) {

	for (auto& c : s){

		if (c < static_cast<char>(33)) c = static_cast<char>(33);
		if (c > static_cast<char>(33+qv_max)) c = static_cast<char>(33+qv_max);
	}
}

inline char getQual(const double score, const size_t qv_min = 0, const size_t qv_max = 40) {

	const char phred_base_std = static_cast<char>(33);
	const char phred_scale_std = static_cast<char>(qv_max);

	const double qv_score = min(score, 1.0) * static_cast<double>(phred_scale_std - qv_min);

	return static_cast<char>(qv_score + phred_base_std + qv_min);
}

inline double getScore(const char c, const size_t qv_min = 0, const size_t qv_max = 40) {

	const char phred_base_std =  static_cast<char>(33);
	const char phred_scale_std =  static_cast<char>(qv_max);

	const double qv_score = static_cast<double>(c - phred_base_std - qv_min);

	return min(qv_score / static_cast<double>(phred_scale_std - qv_min), 1.0);
}

inline bool isValidHap(const PairID& hap_ids, const uint64_t hap_id) {

	return (hap_ids.isEmpty() || hap_ids.contains(hap_id));
}

inline pair<size_t, size_t> getMinMaxLength(const size_t l, const double len_factor) {

	return {static_cast<size_t>(max(l - (l * len_factor), 1.0)), static_cast<size_t>(max(l + (l * len_factor), 1.0))};
}

inline int count_prints(const string& s){

    return std::count_if(s.begin(), s.end(), [](unsigned char c){ return std::isprint(c); });
}

size_t approximate_log2(size_t v);

bool check_files(vector<string>& v_fn, const bool check_files_format, const bool verbose = false);
bool check_files(const string& fn, const bool check_files_format, const bool verbose = false);

// From https://stackoverflow.com/questions/9345087/choose-m-elements-randomly-from-a-vector-containing-n-elements
// This is a Fisher-Yates shuffle, similar to random_Shuffle, except it stops after n iterations
// Note that this function modifies the input vector
template<class BidiIter>
BidiIter random_unique(BidiIter begin, BidiIter end, size_t n) {

    size_t left = std::distance(begin, end);

    while (n--) {

        BidiIter r = begin;

        std::advance(r, rand()%left);
        std::swap(*begin, *r);

        ++begin;
        --left;
    }

    return begin;
}


#endif
