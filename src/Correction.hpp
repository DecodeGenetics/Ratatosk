#ifndef RATATOSK_CORRECTION_HPP
#define RATATOSK_CORRECTION_HPP

#include <chrono> 
#include <iostream>
#include <sstream>
#include <unordered_map>

#include <bifrost/CompactedDBG.hpp>

#include "Alignment.hpp"
#include "Common.hpp"
#include "GraphTraversal.hpp"
#include "Path.hpp"
#include "ResultCorrection.hpp"
#include "TinyBloomFilter.hpp"
#include "UnitigData.hpp"

pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> extractSemiWeakPaths(	const Correct_Opt& opt, const string& s,
																				const TinyBloomFilter<uint32_t>& bf, const PairID& r,
																				const pair<size_t, const_UnitigMap<UnitigData>>& um_solid_start,
																				const pair<size_t, const_UnitigMap<UnitigData>>& um_solid_end, 
																				const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_weak, size_t i_weak,
																				const size_t min_cov_vertex, const bool long_read_correct);

pair<string, string> correctSequence(	const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const string& seq, const string& qual,
										const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_solid, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_weak,
										const bool long_read_correct, const Roaring* all_partitions);

inline bool hasUniquePosition(const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v, const size_t pos){

	return (((pos == 0) || (v[pos].first != v[pos-1].first)) && ((pos == v.size()-1) || (v[pos].first != v[pos+1].first)));
}

#endif