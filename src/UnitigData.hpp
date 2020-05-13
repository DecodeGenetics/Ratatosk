#ifndef RATATOSK_UNITIGDATA_HPP
#define RATATOSK_UNITIGDATA_HPP

#include <iostream>

//#include <bifrost/CompactedDBG.hpp>
#include "CompactedDBG.hpp"

#include "PairID.hpp"
#include "Common.hpp"

class UnitigData : public CDBG_Data_t<UnitigData> {

    public:

        UnitigData() {

            clear();
        }

        UnitigData(const UnitigData& o) : cycle_kmer_cov(o.cycle_kmer_cov), connected_comp_id(o.connected_comp_id), read_ids(o.read_ids), ambiguity_ids(o.ambiguity_ids) {}

        UnitigData(UnitigData&& o) : cycle_kmer_cov(o.cycle_kmer_cov), connected_comp_id(o.connected_comp_id), read_ids(move(o.read_ids)), ambiguity_ids(move(o.ambiguity_ids)) {

            o.clear();
        }

        UnitigData& operator=(const UnitigData& o){

            if (this != &o){

                cycle_kmer_cov = o.cycle_kmer_cov;
                connected_comp_id = o.connected_comp_id;
                read_ids = o.read_ids;
                ambiguity_ids = o.ambiguity_ids;
            }

            return *this;
        }

        UnitigData& operator=(UnitigData&& o){

            if (this != &o){

                cycle_kmer_cov = o.cycle_kmer_cov;
                connected_comp_id = o.connected_comp_id;
                read_ids = move(o.read_ids);
                ambiguity_ids = move(o.ambiguity_ids);

                o.clear();
            }

            return *this;
        }

        void clear(const bool clear_partitions = true){

            if (clear_partitions) connected_comp_id = 0;

            cycle_kmer_cov = 0;

            read_ids.clear();
            ambiguity_ids.clear();
        }

        void clear(const UnitigMap<UnitigData>& um_dest){

            clear();
        }

        void concat(const UnitigMap<UnitigData>& um_dest, const UnitigMap<UnitigData>& um_src){

            const size_t k = um_dest.getGraph()->getK();

            const UnitigData* pli_dest = um_dest.getData();
            const UnitigData* pli_src = um_src.getData();

            vector<pair<size_t, char>> v_dest = pli_dest->get_ambiguity_char(um_dest);
            vector<pair<size_t, char>> v_src = pli_src->get_ambiguity_char(um_src);

            v_dest.reserve(v_dest.size() + v_src.size());

            for (const auto& p : v_src) v_dest.push_back({um_dest.size - k + 1 + p.first, p.second});

            if ((pli_dest->get_readID().cardinality() == 0) || (pli_src->get_readID().cardinality() == 0)) read_ids.clear();
            else {

                read_ids = pli_dest->get_readID() | pli_src->get_readID();

                read_ids.runOptimize();
            }

            cycle_kmer_cov = 0;
            connected_comp_id = max(pli_src->getConnectedComp(), pli_dest->getConnectedComp());

            increaseCoverage(pli_dest->getCoverage());
            increaseCoverage(pli_src->getCoverage());

            string new_um_seq = um_dest.mappedSequenceToString() + um_src.mappedSequenceToString().substr(k-1, um_src.size - k + 1);

            for (const auto& p : v_dest){

                if (isDNA(new_um_seq[p.first])) new_um_seq[p.first] = p.second;
                else {

                    bool is_a_q = false, is_c_q = false, is_g_q = false, is_t_q = false;
                    bool is_a_t = false, is_c_t = false, is_g_t = false, is_t_t = false;

                    getAmbiguityRev(new_um_seq[p.first], is_a_t, is_c_t, is_g_t, is_t_t);
                    getAmbiguityRev(p.second, is_a_q, is_c_q, is_g_q, is_t_q);

                    new_um_seq[p.first] = getAmbiguity(is_a_q || is_a_t, is_c_q || is_c_t, is_g_q || is_g_t, is_t_q || is_t_t);
                }
            }

            ambiguity_ids.clear();

            for (size_t i = 0; i < new_um_seq.length(); ++i){

                if (!isDNA(new_um_seq[i])) add_ambiguity_char(i, new_um_seq[i]);
            }

            ambiguity_ids.runOptimize();
        }

        void merge(const UnitigMap<UnitigData>& um_dest, const const_UnitigMap<UnitigData, void>& um_src){

            const UnitigData* pli_dest = um_dest.getData();
            const UnitigData* pli_src = um_src.getData();

            read_ids |= pli_src->get_readID();

            read_ids.runOptimize();

            connected_comp_id = max(pli_src->getConnectedComp(), pli_dest->getConnectedComp());

            increaseCoverage(pli_src->getKmerCoverage(um_src) * um_src.len);

            if (pli_src->isInCycle()) setCycle();
        }

        void extract(const UnitigMap<UnitigData>& um_src, bool last_extraction) {

            const UnitigData* pli_src = um_src.getData();
            const vector<pair<size_t, char>> v_amb = pli_src->get_ambiguity_char(um_src);

            clear();

            for (const auto p : v_amb) add_ambiguity_char(p.first, p.second);

            connected_comp_id = pli_src->getConnectedComp();
            read_ids = pli_src->get_readID();

            read_ids.runOptimize();
            ambiguity_ids.runOptimize();

            increaseCoverage(pli_src->getKmerCoverage(um_src) * um_src.len);

            if (pli_src->isInCycle()) setCycle();
        }

        string serialize(const const_UnitigMap<UnitigData, void>& um_src) const {

            return string(  "KC:Z:" + std::to_string(um_src.getData()->getKmerCoverage(um_src)) + '\t' + 
                            "RC:Z:" + std::to_string(um_src.getData()->get_readID().cardinality()) + '\t' +
                            "PU:Z:" + std::to_string(um_src.getData()->getConnectedComp())
                            );
        }

        inline void resetCoverage() {

            cycle_kmer_cov &= 0x8000000000000000ULL;
        }

        inline void increaseCoverage(const size_t count){

            if (count <= (0x7fffffffffffffffULL - (cycle_kmer_cov & 0x7fffffffffffffffULL))) cycle_kmer_cov += count;
            else cycle_kmer_cov = (cycle_kmer_cov & 0x8000000000000000ULL) | 0x7fffffffffffffffULL;
        }

        inline void  decreaseKmerCoverage(const size_t count){

            if (count <= (cycle_kmer_cov & 0x7fffffffffffffffULL)) cycle_kmer_cov -= count;
            else cycle_kmer_cov &= 0x8000000000000000ULL;
        }

        inline size_t getCoverage() const {

            return (cycle_kmer_cov & 0x7fffffffffffffffULL);
        }

        inline size_t getCoverage(const size_t nb_km) const {

            return (static_cast<double>(cycle_kmer_cov & 0x7fffffffffffffffULL) / static_cast<double>(nb_km));
        }

        inline double getKmerCoverage(const const_UnitigMap<UnitigData>& um) const {

            return (static_cast<double>(cycle_kmer_cov & 0x7fffffffffffffffULL) / static_cast<double>(um.size - um.getGraph()->getK() + 1));
        }

        inline void setCycle(){

            cycle_kmer_cov |= 0x8000000000000000ULL;
        }

        inline void unsetCycle(){

            cycle_kmer_cov &= 0x7fffffffffffffffULL;
        }

        inline bool isInCycle() const {

            return static_cast<bool>(cycle_kmer_cov >> 63);
        }

        inline const PairID& get_readID() const {

            return read_ids;
        }   

        inline PairID& get_readID() {

            return read_ids;
        }

        inline size_t getConnectedComp() const {

            return connected_comp_id;
        }

        inline void setConnectedComp(const size_t id) {

            connected_comp_id = id;
        }

        inline bool write(ostream& stream_out) const {

            if (stream_out.good()) stream_out.write(reinterpret_cast<const char*>(&cycle_kmer_cov), sizeof(size_t));
            if (stream_out.good()) stream_out.write(reinterpret_cast<const char*>(&connected_comp_id), sizeof(size_t));

            if (stream_out.good()) read_ids.write(stream_out);
            if (stream_out.good()) ambiguity_ids.write(stream_out);

            return stream_out.good();
        }

        inline bool read(istream& stream_in){

            if (stream_in.good()) stream_in.read(reinterpret_cast<char*>(&cycle_kmer_cov), sizeof(size_t));
            if (stream_in.good()) stream_in.read(reinterpret_cast<char*>(&connected_comp_id), sizeof(size_t));

            if (stream_in.good()) read_ids.read(stream_in);
            if (stream_in.good()) ambiguity_ids.read(stream_in);

            return stream_in.good();
        }

        inline void clear_ambiguity_char(){

            ambiguity_ids.clear();     
        }

        inline void add_ambiguity_char(const size_t pos, const char c){

            ambiguity_ids.add((pos << 4) + getAmbiguityIndex(c));
        }

        vector<pair<size_t, char>> get_ambiguity_char(const const_UnitigMap<UnitigData>& um) const {

            vector<pair<size_t, char>> v, v_tmp = get_ambiguity_char();

            const size_t sz = um.len + um.getGraph()->getK() - 1;
            const size_t end = um.dist + sz;

            if (um.strand){

                for (const auto p : v_tmp){

                    if ((p.first >= um.dist) && (p.first < end)) v.push_back({p.first - um.dist, p.second});
                }
            }
            else {

                for (auto it = v_tmp.rbegin(); it != v_tmp.rend(); ++it) {

                    if ((it->first >= um.dist) && (it->first < end)) v.push_back({sz - (it->first - um.dist) - 1, reverse_complement(it->second)});
                }
            }

            return v;
        }

        inline bool has_ambiguity_char() const {

            return !ambiguity_ids.isEmpty();
        }

        inline void runOptimizeAmbiguityChar() {

            ambiguity_ids.runOptimize();
        }

        inline void print(const const_UnitigMap<UnitigData>& um) const {

            const_UnitigMap<UnitigData> l_um(um);

            l_um.dist = 0;
            l_um.len = l_um.size - l_um.getGraph()->getK() + 1;
            l_um.strand = true;

            const vector<pair<size_t, char>> v = get_ambiguity_char(l_um);

            cout << "- Unitig: " << l_um.referenceUnitigToString() << endl;
            cout << "- Unitig length: " << l_um.size << endl;
            cout << "- Coverage: " << getCoverage() << endl;
            cout << "- K-mer coverage: " << getKmerCoverage(l_um) << endl;
            cout << "- Read mapping: " << get_readID().cardinality() << endl;
            cout << "- Short cycle: " << isInCycle() << endl;

            if (v.empty()) cout << "- SNPs: None" << endl;
            else {

                cout << "- SNPs:" << endl;

                for (const auto& p : v) cout << p.first << " " << p.second << endl;
            }

            cout << " " << endl;
        }

    private:

        struct compare_ambiguity {

            bool operator()(const pair<size_t, char>& a, const pair<size_t, char>& b) const {

                return (a.first < b.first) || ((a.first == b.first) && (a.second < b.second));
            }
        };

        vector<pair<size_t, char>> get_ambiguity_char() const {

            vector<pair<size_t, char>> v;

            for (const uint32_t id : ambiguity_ids) v.push_back({id >> 4, getAmbiguityIndexRev(id & 0xF)});
            
            sort(v.begin(), v.end(), compare_ambiguity());

            return v;
        }

        size_t cycle_kmer_cov;
        size_t connected_comp_id;

        PairID read_ids;
        PairID ambiguity_ids;
};

#endif