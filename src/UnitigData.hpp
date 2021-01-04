#ifndef RATATOSK_UNITIGDATA_HPP
#define RATATOSK_UNITIGDATA_HPP

#include <iostream>

#include "CompactedDBG.hpp"

#include "PairID.hpp"
#include "Common.hpp"

class UnitigData : public CDBG_Data_t<UnitigData> {

    public:

        UnitigData() {

            clear();
        }

        UnitigData(const UnitigData& o) : kmCov_cardBranches(o.kmCov_cardBranches), connected_comp_id(o.connected_comp_id), read_ids(o.read_ids), ambiguity_ids(o.ambiguity_ids), hap_ids(o.hap_ids) {}

        UnitigData(UnitigData&& o) : kmCov_cardBranches(o.kmCov_cardBranches), connected_comp_id(o.connected_comp_id), read_ids(move(o.read_ids)), ambiguity_ids(move(o.ambiguity_ids)), hap_ids(move(o.hap_ids)) {

            o.clear();
        }

        UnitigData& operator=(const UnitigData& o){

            if (this != &o){

                kmCov_cardBranches = o.kmCov_cardBranches;
                connected_comp_id = o.connected_comp_id;
                read_ids = o.read_ids;
                ambiguity_ids = o.ambiguity_ids;
                hap_ids = o.hap_ids;
            }

            return *this;
        }

        UnitigData& operator=(UnitigData&& o){

            if (this != &o){

                kmCov_cardBranches = o.kmCov_cardBranches;
                connected_comp_id = o.connected_comp_id;
                read_ids = move(o.read_ids);
                ambiguity_ids = move(o.ambiguity_ids);
                hap_ids = move(o.hap_ids);

                o.clear();
            }

            return *this;
        }

        void clear(){

            connected_comp_id = 0;
            kmCov_cardBranches = 0;

            read_ids.clear();
            ambiguity_ids.clear();
            hap_ids.clear();
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

            read_ids = pli_dest->get_readID() | pli_src->get_readID();
            hap_ids = pli_dest->get_hapID() | pli_src->get_hapID();

            read_ids.runOptimize();
            hap_ids.runOptimize();

            kmCov_cardBranches = 0;
            connected_comp_id = max(pli_src->getConnectedComp(), pli_dest->getConnectedComp());

            //increaseCoverage(pli_dest->getCoverage());
            //increaseCoverage(pli_src->getCoverage());
            increasePhasedCoverage(pli_dest->getPhasedCoverage() + pli_src->getPhasedCoverage());
            increaseUnphasedCoverage(pli_dest->getUnphasedCoverage() + pli_src->getUnphasedCoverage());

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
            hap_ids |= pli_src->get_hapID();

            read_ids.runOptimize();
            hap_ids.runOptimize();

            connected_comp_id = max(pli_src->getConnectedComp(), pli_dest->getConnectedComp());

            //increaseCoverage(pli_src->getKmerCoverage(um_src) * um_src.len);
            increasePhasedCoverage(pli_src->getPhasedKmerCoverage(um_src) * um_src.len);
            increaseUnphasedCoverage(pli_src->getUnphasedKmerCoverage(um_src) * um_src.len);
        }

        void extract(const UnitigMap<UnitigData>& um_src, bool last_extraction) {

            const UnitigData* pli_src = um_src.getData();
            const vector<pair<size_t, char>> v_amb = pli_src->get_ambiguity_char(um_src);

            for (const auto p : v_amb) add_ambiguity_char(p.first, p.second);

            connected_comp_id = pli_src->getConnectedComp();
            read_ids = pli_src->get_readID();
            hap_ids = pli_src->get_hapID();

            read_ids.runOptimize();
            ambiguity_ids.runOptimize();
            hap_ids.runOptimize();

            //increaseCoverage(pli_src->getKmerCoverage(um_src) * um_src.len);
            increasePhasedCoverage(pli_src->getPhasedKmerCoverage(um_src) * um_src.len);
            increaseUnphasedCoverage(pli_src->getUnphasedKmerCoverage(um_src) * um_src.len);
        }

        string serialize(const const_UnitigMap<UnitigData, void>& um_src) const {

            return string(  "KC:Z:" + std::to_string(um_src.getData()->getKmerCoverage(um_src)) + '\t' + 
                            "RC:Z:" + std::to_string(um_src.getData()->get_readID().cardinality()) + '\t' +
                            "PU:Z:" + std::to_string(um_src.getData()->getConnectedComp())
                            );
        }

        inline bool write(ostream& stream_out) const {

            if (stream_out.good()) stream_out.write(reinterpret_cast<const char*>(&kmCov_cardBranches), sizeof(size_t));
            if (stream_out.good()) stream_out.write(reinterpret_cast<const char*>(&connected_comp_id), sizeof(size_t));

            if (stream_out.good()) read_ids.write(stream_out);
            if (stream_out.good()) ambiguity_ids.write(stream_out);
            if (stream_out.good()) hap_ids.write(stream_out);

            return stream_out.good();
        }

        inline bool read(istream& stream_in){

            if (stream_in.good()) stream_in.read(reinterpret_cast<char*>(&kmCov_cardBranches), sizeof(size_t));
            if (stream_in.good()) stream_in.read(reinterpret_cast<char*>(&connected_comp_id), sizeof(size_t));

            if (stream_in.good()) read_ids.read(stream_in);
            if (stream_in.good()) ambiguity_ids.read(stream_in);
            if (stream_in.good()) hap_ids.read(stream_in);

            return stream_in.good();
        }

        inline void print(const const_UnitigMap<UnitigData>& um) const {

            const_UnitigMap<UnitigData> l_um(um);

            l_um.dist = 0;
            l_um.len = l_um.size - l_um.getGraph()->getK() + 1;
            l_um.strand = true;

            const vector<pair<size_t, char>> v = get_ambiguity_char(l_um);

            cout << "- Unitig: " << l_um.referenceUnitigToString() << endl;
            cout << "- Unitig length: " << l_um.size << endl;
            cout << "- Successors: " << l_um.getSuccessors().cardinality() << endl;
            cout << "- Predecessors: " << l_um.getPredecessors().cardinality() << endl;
            cout << "- Coverage: " << getPhasedCoverage() << " (phased) + " << getUnphasedCoverage() << " (unphased) = " << getCoverage() << endl;
            cout << "- Kmer coverage: " << getPhasedKmerCoverage(l_um) << " (phased) + " << getUnphasedKmerCoverage(l_um) << " (unphased) = " << getKmerCoverage(l_um) << endl;
            cout << "- Read mapping: " << get_readID().cardinality() << endl;

            if (v.empty()) cout << "- SNPs: None" << endl;
            else {

                cout << "- SNPs:" << endl;

                for (const auto& p : v) cout << p.first << " " << p.second << endl;
            }

            cout << " " << endl;

            if (hap_ids.isEmpty()) cout << "Phasing: None" << endl;
            else {

                cout << "- Phasing:" << endl;

                for (const uint32_t hapid : hap_ids) cout << (hapid >> 1) << " " << (hapid & 0x1ULL) << endl;
            }

                cout << " " << endl;
        }

        inline void resetCoverage() {

            kmCov_cardBranches &= 0xc000000000000000ULL;
        }

        inline void increasePhasedCoverage(const size_t count){

            const size_t no_change = kmCov_cardBranches & 0x7fffffff80000000ULL;
            const size_t phased_count = kmCov_cardBranches & 0x000000007fffffffULL;

            if ((phased_count + count) <= 0x000000007fffffffULL) kmCov_cardBranches = no_change | (phased_count + count);
            else kmCov_cardBranches = no_change | 0x000000007fffffffULL;
        }

        inline void increaseUnphasedCoverage(const size_t count){

            const size_t no_change = kmCov_cardBranches & 0xc00000007fffffffULL;
            const size_t unphased_count = (kmCov_cardBranches >> 31) & 0x000000007fffffffULL;

            if ((unphased_count + count) <= 0x000000007fffffffULL) kmCov_cardBranches = no_change | ((unphased_count + count) << 31);
            else kmCov_cardBranches = no_change | (0x000000007fffffffULL << 31);
        }

        inline void decreasePhasedCoverage(const size_t count){

            const size_t no_change = kmCov_cardBranches & 0x7fffffff80000000ULL;
            const size_t phased_count = kmCov_cardBranches & 0x000000007fffffffULL;

            if (phased_count < count) kmCov_cardBranches = no_change;
            else kmCov_cardBranches = no_change | (phased_count - count);
        }

        inline void decreaseUnphasedCoverage(const size_t count){

            const size_t no_change = kmCov_cardBranches & 0xc00000007fffffffULL;
            const size_t unphased_count = (kmCov_cardBranches >> 31) & 0x000000007fffffffULL;

            if (unphased_count < count) kmCov_cardBranches = no_change;
            else kmCov_cardBranches = no_change | ((unphased_count - count) << 31);
        }

        inline size_t getPhasedCoverage() const {

            return (kmCov_cardBranches & 0x000000007fffffffULL);
        }

        inline size_t getUnphasedCoverage() const {

            return (kmCov_cardBranches >> 31) & 0x000000007fffffffULL;
        }

        inline size_t getCoverage() const {

            return getPhasedCoverage() + getUnphasedCoverage();
        }

        inline double getPhasedKmerCoverage(const const_UnitigMap<UnitigData>& um) const {

            return round(static_cast<double>(getPhasedCoverage()) / static_cast<double>(um.size - um.getGraph()->getK() + 1));
        }

        inline double getUnphasedKmerCoverage(const const_UnitigMap<UnitigData>& um) const {

            return round(static_cast<double>(getUnphasedCoverage()) / static_cast<double>(um.size - um.getGraph()->getK() + 1));
        }

        inline double getKmerCoverage(const const_UnitigMap<UnitigData>& um) const {

            return round(static_cast<double>(getCoverage()) / static_cast<double>(um.size - um.getGraph()->getK() + 1));
        }

        inline void setBranching(const bool isBranching) {

            kmCov_cardBranches &= 0x7fffffffffffffffULL;
            kmCov_cardBranches |= static_cast<size_t>(isBranching) << 63;
        }

        inline bool isBranching() const {

            return static_cast<bool>(kmCov_cardBranches >> 63);
        }

        inline const PairID& get_readID() const {

            return read_ids;
        }   

        inline PairID& get_readID() {

            return read_ids;
        }

        inline const PairID& get_hapID() const {

            return hap_ids;
        }   

        inline PairID& get_hapID() {

            return hap_ids;
        }

        inline size_t getConnectedComp() const {

            return connected_comp_id;
        }

        inline void setConnectedComp(const size_t id) {

            connected_comp_id = id;
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

        size_t kmCov_cardBranches;
        size_t connected_comp_id;

        PairID read_ids;
        PairID ambiguity_ids;
        PairID hap_ids;
};

#endif