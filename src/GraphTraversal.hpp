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

	unordered_map<Kmer, size_t, KmerHash> m_km_fw, m_km_bw;
	queue<pair<const_UnitigMap<UnitigData>, pair<size_t, size_t>>> q_um_fw, q_um_bw;
	PairID p_ids_fw, p_ids_bw;
};

pair<vector<Path<UnitigData>>, bool> explorePathsBFS(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const WeightsPairID& w_pid, const const_UnitigMap<UnitigData>& um_s,
														const bool long_read_correct, const uint64_t hap_id);

pair<vector<Path<UnitigData>>, bool> explorePathsBFS2(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const WeightsPairID& w_pid, const const_UnitigMap<UnitigData>& um_s, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_e,
														const bool long_read_correct, const uint64_t hap_id);


pair<double, double> exploreSubGraph(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
										const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e, const size_t level,
										vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id);

pair<double, double> exploreSubGraph(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
										const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e, const size_t level,
										vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id);

pair<double, double> exploreSubGraphLong(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
											const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e,
											vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id);

pair<double, double> exploreSubGraphLong(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
											const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e,
											vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id);

pair<double, double> getScorePath(const Path<UnitigData>& path, const WeightsPairID& w_pid, const uint64_t hap_id, const char* ref, const size_t ref_len);

vector<Path<UnitigData>> selectMostContiguous(const vector<Path<UnitigData>>& v_paths, const WeightsPairID& w_pid);
vector<Path<UnitigData>> selectMostContiguous(const vector<Path<UnitigData>>& v_paths, const size_t min_cov, const uint64_t hap_id);

vector<pair<size_t, char>> getAmbiguityVector(const vector<const_UnitigMap<UnitigData>>& v_um, const size_t k);
vector<pair<size_t, char>> getAmbiguityVector(const Path<UnitigData>::PathOut& path, const size_t k, const bool rm_solid_amb = false);

bool isValidSNPcandidate(local_graph_traversal& lgt, const const_UnitigMap<UnitigData>& um_a, const const_UnitigMap<UnitigData>& um_b, const size_t limit_length_path = 10000, const size_t limit_sz_stack = 65536);

void annotateConnectedComponents(CompactedDBG<UnitigData>& dbg);

#endif