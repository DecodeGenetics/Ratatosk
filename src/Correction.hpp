#ifndef RATATOSK_CORRECTION_HPP
#define RATATOSK_CORRECTION_HPP

#include <chrono> 
#include <iostream>
#include <sstream>
#include <unordered_map>

#include "CompactedDBG.hpp"

#include "Alignment.hpp"
#include "Common.hpp"
#include "GraphTraversal.hpp"
#include "Path.hpp"
#include "ResultCorrection.hpp"
#include "TinyBloomFilter.hpp"
#include "UnitigData.hpp"

pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> extractSemiWeakPaths(	const Correct_Opt& opt, const string& s,
																				const TinyBloomFilter<uint32_t>& bf, const PairID& r, const WeightsPairID& w_pid,
																				const pair<size_t, const_UnitigMap<UnitigData>>& um_solid_start,
																				const pair<size_t, const_UnitigMap<UnitigData>>& um_solid_end, 
																				const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_weak, size_t i_weak,
																				const bool long_read_correct, const uint64_t hap_id);

pair<string, string> correctSequence(	const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const string& seq, const string& qual,
										const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_solid, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_weak,
										const bool long_read_correct, const Roaring* all_partitions, const uint64_t hap_id, const pair<HapReads, HapReads>& hap_reads);

inline bool hasUniquePosition(const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v, const size_t pos){

	return (((pos == 0) || (v[pos].first != v[pos-1].first)) && ((pos == v.size()-1) || (v[pos].first != v[pos+1].first)));
}

inline bool hasMinQual(const string& s, const string& q, const size_t start, const size_t end, const char min_q){

	bool hasQual = true;

	for (size_t i = start; (i < end) && hasQual; ++i) hasQual = ((q[i] >= min_q) || !isDNA(s[i]));

	return hasQual;
}

#endif