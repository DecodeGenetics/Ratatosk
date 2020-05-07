#ifndef RATATOSK_GRAPHTRAVERSAL_HPP
#define RATATOSK_GRAPHTRAVERSAL_HPP

#include <chrono>
#include <iostream>
#include <math.h>
#include <queue>
#include <stack>
#include <set>

#include <bifrost/CompactedDBG.hpp>

#include "Common.hpp"
#include "Alignment.hpp"
#include "Path.hpp"
#include "TinyBloomFilter.hpp"
#include "UnitigData.hpp"

pair<vector<Path<UnitigData>>, bool> explorePathsBFS(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const TinyBloomFilter<uint32_t>& bf, const PairID& r,
														const const_UnitigMap<UnitigData>& um_s,
														const size_t min_cov_vertex, const bool long_read_correct);

pair<vector<Path<UnitigData>>, bool> explorePathsBFS2(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const TinyBloomFilter<uint32_t>& bf, const PairID& r,
														const const_UnitigMap<UnitigData>& um_s, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_e,
														const size_t min_cov_vertex, const bool long_read_correct);

void exploreSubGraph(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
					const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e,
					const Path<UnitigData>& p, const PairID& r, const size_t level, const size_t min_cov_vertex,
					double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& over_represented);

void exploreSubGraph(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
					const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e,
					const Path<UnitigData>& p, const PairID& r, const size_t level, const size_t min_cov_vertex,
					double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& over_represented);

void exploreSubGraphLong(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
						const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e,
						const Path<UnitigData>& p, const PairID& r, const size_t min_cov_vertex,
						double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& semi_weak_successors);

void exploreSubGraphLong(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
						const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e,
						const Path<UnitigData>& p, const PairID& r, const size_t min_cov_vertex,
						double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& semi_weak_successors);


#endif