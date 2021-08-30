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
																													const string& s, const string& q,
																													const bool long_read_correct, const uint64_t hap_id,
																													unordered_map<Kmer, vector<const_UnitigMap<UnitigData>>, KmerHash>& m_km_um,
																													const bool fill_map);

void addCoverage(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, HapReads& hap_reads, const bool long_read_correct);
void addPhasing(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, HapReads& hap_r, const bool long_read_phasing, const bool mapHapReads);

void detectSNPs(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt);
void detectShortCycles(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt);

void writeGraphData(const string& output_filename, const CompactedDBG<UnitigData>& dbg, const bool verbose = false);
void writeGraphPairID(CompactedDBG<UnitigData>& dbg, const string& out_fn, const Correct_Opt& opt);
bool readGraphData(const string& input_filename, CompactedDBG<UnitigData>& dbg, const bool pid_only, const bool verbose = false);

unordered_map<Kmer, vector<streampos>, KmerHash> mergeDiskPairIDs(const string& fn);

unordered_map<uint64_t, uint64_t, CustomHashUint64_t> getPairIDs(const vector<string>& filenames, const uint64_t seed);
void packPairIDs(const CDBG_Build_opt& opt, CompactedDBG<UnitigData>& dbg, unordered_map<uint64_t, uint64_t, CustomHashUint64_t>& name_hmap, const size_t seed, const size_t nb_threads);

pair<BlockedBloomFilter, unordered_set<uint64_t, CustomHashUint64_t>> buildBBF(const vector<string>& v_filenames_in, const Correct_Opt& opt, const size_t k, const size_t g, const bool long_reads);
BlockedBloomFilter buildBBF(const CompactedDBG<UnitigData>& dbg);

string retrieveMissingReads(const Correct_Opt& opt);

size_t estimateHaplotypeCoverage(const CompactedDBG<UnitigData>& dbg);

inline void annotateBranchingUnitigs(CompactedDBG<UnitigData>& dbg) {

	for (auto& um : dbg) um.getData()->setBranching((um.getPredecessors().cardinality() > 1) || (um.getSuccessors().cardinality() > 1));
}

//Roaring* createPartitions(CompactedDBG<UnitigData>& dbg, const vector<Kmer>& centroids, const double mean_len_to_neighbor, const size_t nb_threads, const bool verbose);
//void expandConnectedComponent(CompactedDBG<UnitigData>& dbg, const size_t id_mapped, const size_t id_visited);

string phasing(const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const string& s, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_solid);

size_t getMaxKmerCoverage(const CompactedDBG<UnitigData>& dbg, const double top_ratio);

#endif