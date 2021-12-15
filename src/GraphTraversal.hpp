#ifndef RATATOSK_GRAPHTRAVERSAL_HPP
#define RATATOSK_GRAPHTRAVERSAL_HPP

#include <chrono>
#include <iostream>
#include <math.h>
#include <queue>
#include <stack>
#include <set>

#include "CompactedDBG.hpp"

#include "Common.hpp"
#include "Alignment.hpp"
#include "Path.hpp"
#include "TinyBloomFilter.hpp"
#include "UnitigData.hpp"

struct info_traversal {

	Path<UnitigData> p;
	size_t l;

	info_traversal(const Path<UnitigData>& p_, const size_t l_) : p(p_), l(l_) {};
	info_traversal(const Path<UnitigData>& p_) : p(p_), l(0) {};
	info_traversal(const size_t l_) : l(l_) {};
	info_traversal() : l(0) {};
};

struct local_graph_traversal {

	unordered_map<Kmer, const SharedPairID*, KmerHash> m_km;
	queue<const_UnitigMap<UnitigData>> q_um;
};

pair<vector<Path<UnitigData>>, bool> explorePathsBFS(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const WeightsPairID& w_pid, const const_UnitigMap<UnitigData>& um_s,
														const bool long_read_correct, const uint64_t hap_id);

pair<vector<Path<UnitigData>>, bool> explorePathsBFS2(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const WeightsPairID& w_pid, const const_UnitigMap<UnitigData>& um_s, const const_UnitigMap<UnitigData>& um_e,
														const bool long_read_correct, const uint64_t hap_id);


pair<double, double> exploreSubGraph(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len, const size_t max_len_path,
										const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e, const size_t level,
										vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id,
										unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid);

pair<double, double> exploreSubGraphLong(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len, const size_t max_len_path,
											const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e,
											vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id,
											unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid);

pair<double, double> getScorePath(const Correct_Opt& opt, const Path<UnitigData>& path, const char* ref, const size_t ref_len, const bool terminal, const WeightsPairID& w_pid);
pair<double, double> getScorePath(	const Correct_Opt& opt, const Path<UnitigData>& path, const char* ref, const size_t ref_len, const bool terminal, const WeightsPairID& w_pid,
									unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid);

string getScorePath(const Correct_Opt& opt, const Path<UnitigData>& path, const char* ref, const size_t ref_len, const double score_best, const double score_second_best);

vector<Path<UnitigData>> fixRepeats(const Correct_Opt& opt, const vector<Path<UnitigData>>& v_path, const char* ref, const size_t ref_len);

vector<Path<UnitigData>> selectMostContiguous(const vector<Path<UnitigData>>& v_paths, const WeightsPairID& w_pid);
vector<Path<UnitigData>> selectMostContiguous(const vector<Path<UnitigData>>& v_paths, const size_t min_cov, const uint64_t hap_id);

vector<pair<size_t, char>> getAmbiguityVector(const vector<const_UnitigMap<UnitigData>>& v_um, const size_t k);
vector<pair<size_t, char>> getAmbiguityVector(const Path<UnitigData>::PathOut& path, const size_t k, const bool rm_solid_amb = false);

bool isValidSNPcandidate(local_graph_traversal& lgt_fw, local_graph_traversal& lgt_bw,
						const const_UnitigMap<UnitigData>& um_a, const const_UnitigMap<UnitigData>& um_b,
						const size_t min_cov = 2, const size_t limit_sz_stack = 65536);

#endif