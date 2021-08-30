#ifndef RATATOSK_ALIGNMENT_HPP
#define RATATOSK_ALIGNMENT_HPP

#include <iostream>
#include <cfloat>
#include <cmath>

#include "edlib.h"
#include "Path.hpp"
#include "ResultCorrection.hpp"
#include "TinyBloomFilter.hpp"
#include "UnitigData.hpp"

pair<int, int> selectBestAlignment(const vector<Path<UnitigData>>& candidates, const char* ref, const size_t len_ref);

pair<int, int> selectBestPrefixAlignment(const char* ref, const size_t ref_len, const vector<Path<UnitigData>>& candidates, const double cut_threshold_norm_edit = -1.0);
pair<int, int> selectBestPrefixAlignment(const char* ref, const size_t ref_len, const vector<const Path<UnitigData>*>& candidates, const double cut_threshold_norm_edit = -1.0);

pair<int, pair<int, int>> selectBestPrefixAlignment2(const char* ref, const size_t ref_len, const vector<Path<UnitigData>>& candidates, const bool get_end_location);

pair<string, string> generateConsensus(const ResultCorrection* fw_s, const ResultCorrection* bw_s, const string& ref_seq, const double max_norm_edit_distance = -1.0);

void fixAmbiguity(	const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const WeightsPairID& w_pid,
					string& query, string& quality,
					const char* ref_seq, const char* ref_qual, const size_t ref_len,
					const uint64_t hap_id, const vector<pair<size_t, char>>& v_ambiguity);

pair<int, int> selectBestSubstringAlignment(const char* ref, const size_t ref_len, const vector<Path<UnitigData>>& candidates, const double cut_threshold_norm_edit = -1.0);

/*vector<pair<size_t, char>> filterAmbiguity(	const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const WeightsPairID& w_pid,
											const string& query, const vector<pair<size_t, char>>& v_ambiguity);*/

string fixSNPs(const Correct_Opt& opt, const CompactedDBG<UnitigData>& dbg, const string& s);

vector<pair<size_t, const_UnitigMap<UnitigData>>> keep_non_overlap(const Correct_Opt& opt, const char* ref, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v);

bool test(const string& query_seq, const size_t query_len, const string& qual_seq, const size_t qual_len, const string& ref_seq, const size_t ref_len);

#endif