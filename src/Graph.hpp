#ifndef RATATOSK_GRAPH_HPP
#define RATATOSK_GRAPH_HPP

#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <stack>

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/roaring.hh>

#include "Correction.hpp"

struct CustomHashUint64_t {

    size_t operator()(const uint64_t v) const {

        return XXH64(&v, sizeof(uint64_t), 0);
    }
};

pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> getSeeds(const Correct_Opt& opt, const CompactedDBG<UnitigData>& dbg,
																													const string& s, const bool long_read_correct,
																													const bool filter_non_unique_inexact_kmers = false);

vector<Kmer> addCoverage(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const bool long_read_correct, const bool map_reads, const size_t id_read_start = 0);

void detectSNPs(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const Roaring* part_neighbors = nullptr);
size_t detectSTRs(CompactedDBG<UnitigData>& dbg, const size_t max_len_cycle, const size_t max_cov_km, const size_t nb_threads = 1, const bool verbose = false);

Roaring* createPartitions(CompactedDBG<UnitigData>& dbg, const vector<Kmer>& centroids, const double mean_len_to_neighbor, const size_t nb_threads, const bool verbose);

void mergeGraphUnmapped(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const string& filenameOut);

//void mergeGrapLongReads(CompactedDBG<UnitigData>& dbg_sr, const Correct_Opt& opt, const Roaring* part_neighbors = nullptr);
void mergeGrapLongReads(CompactedDBG<UnitigData>& dbg_sr, const Correct_Opt& opt, const vector<string>& v_ref_filenames);

void writeGraphData(const string& output_filename, const CompactedDBG<UnitigData>& dbg, const bool verbose = false);
bool readGraphData(const string& input_filename, CompactedDBG<UnitigData>& dbg, const bool verbose = false);

unordered_map<uint64_t, uint64_t, CustomHashUint64_t> getPairIDs(const vector<string>& filenames, const uint64_t seed);
void packPairIDs(const CDBG_Build_opt& opt, CompactedDBG<UnitigData>& dbg, unordered_map<uint64_t, uint64_t, CustomHashUint64_t>& name_hmap, const size_t seed, const size_t nb_threads);

void expandConnectedComponent(CompactedDBG<UnitigData>& dbg, const size_t id_mapped, const size_t id_visited);

//BlockedBloomFilter createBloomFilter(const CDBG_Build_opt& opt, const CompactedDBG<UnitigData>& dbg, const size_t nb_threads = 1);

inline void resetUnitigData(CompactedDBG<UnitigData>& dbg, const bool clear_partitions){

	for (auto& um : dbg) um.getData()->clear(clear_partitions);
}


pair<string, string> quickAndDirty(const CompactedDBG<UnitigData>& dbg, const string& s, const string& q);

#endif