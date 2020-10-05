#ifndef RATATOSK_GRAPH_HPP
#define RATATOSK_GRAPH_HPP

#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <stack>

#include "CompactedDBG.hpp"
#include "roaring.hh"

#include "Correction.hpp"

pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> getSeeds(const Correct_Opt& opt, const CompactedDBG<UnitigData>& dbg,
																													const string& s, const bool long_read_correct, const uint64_t hap_id,
																													const bool filter_non_unique_inexact_kmers = false);

vector<Kmer> addCoverage(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>& read2hap,
						const bool long_read_correct, const bool map_reads, const size_t id_read_start = 0);

void detectSNPs(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const Roaring* part_neighbors = nullptr);
size_t detectSTRs(CompactedDBG<UnitigData>& dbg, const size_t max_len_cycle, const size_t max_cov_km, const size_t nb_threads = 1, const bool verbose = false);

Roaring* createPartitions(CompactedDBG<UnitigData>& dbg, const vector<Kmer>& centroids, const double mean_len_to_neighbor, const size_t nb_threads, const bool verbose);

void addPhasing(CompactedDBG<UnitigData>& dbg, HapReads& hap_r, const Correct_Opt& opt,
				const vector<string>& v_phase_filenames, const vector<string>& v_read_filenames,
				const bool long_read_phasing, const bool mapHapReads);

void writeGraphData(const string& output_filename, const CompactedDBG<UnitigData>& dbg, const bool verbose = false);
bool readGraphData(const string& input_filename, CompactedDBG<UnitigData>& dbg, const bool verbose = false);

unordered_map<uint64_t, uint64_t, CustomHashUint64_t> getPairIDs(const vector<string>& filenames, const uint64_t seed);
void packPairIDs(const CDBG_Build_opt& opt, CompactedDBG<UnitigData>& dbg, unordered_map<uint64_t, uint64_t, CustomHashUint64_t>& name_hmap, const size_t seed, const size_t nb_threads);

void expandConnectedComponent(CompactedDBG<UnitigData>& dbg, const size_t id_mapped, const size_t id_visited);

pair<BlockedBloomFilter, size_t> buildBFF(const vector<string>& v_filenames_in, const Correct_Opt& opt, const size_t k, const size_t g, const bool long_reads);
string retrieveMissingReads(const Correct_Opt& opt);

//vector<pair<size_t, const_UnitigMap<UnitigData>>> searchSeeds(const CompactedDBG<UnitigData>& dbg, const string& seq);

#endif