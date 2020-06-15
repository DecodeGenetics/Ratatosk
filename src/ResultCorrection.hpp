#ifndef RATATOSK_RESULTCORRECTION_HPP
#define RATATOSK_RESULTCORRECTION_HPP

#include <iostream>

class ResultCorrection {

	public: 

		ResultCorrection(const size_t seq_len) : old_seq_len(seq_len), min_cov_vertex(0) {}

		ResultCorrection(const ResultCorrection& o) : 	pos_corrected_old_seq(o.pos_corrected_old_seq), seq(o.seq), qual(o.qual), old_seq_len(o.old_seq_len),
														min_cov_vertex(o.min_cov_vertex), r_s(o.r_s), r_w(o.r_w), r_t(o.r_t) {}

		ResultCorrection(ResultCorrection&& o) :	pos_corrected_old_seq(move(o.pos_corrected_old_seq)), seq(move(o.seq)), qual(move(o.qual)), old_seq_len(o.old_seq_len),
													min_cov_vertex(o.min_cov_vertex), r_s(move(o.r_s)), r_w(move(o.r_w)), r_t(move(o.r_t)) { o.clear(); }

		ResultCorrection operator=(const ResultCorrection& o){

			if (&o != this){

				pos_corrected_old_seq = o.pos_corrected_old_seq;
				old_seq_len = o.old_seq_len;
				min_cov_vertex = o.min_cov_vertex;

				seq = o.seq;
				qual = o.qual;

				r_s = o.r_s;
				r_w = o.r_w;
				r_t = o.r_t;
			}

			return *this;
		}

		ResultCorrection operator=(ResultCorrection&& o){

			if (&o != this){

				pos_corrected_old_seq = move(o.pos_corrected_old_seq);
				old_seq_len = o.old_seq_len;
				min_cov_vertex = o.min_cov_vertex;

				seq = move(o.seq);
				qual = move(o.qual);

				r_s = move(o.r_s);
				r_w = move(o.r_w);
				r_t = move(o.r_t);

				o.clear();
			}

			return *this;
		}

		~ResultCorrection() {

			clear();
		}

		inline void clear() {

			pos_corrected_old_seq = Roaring();

			seq.clear();
			qual.clear();

			r_s.clear();
			r_w.clear();
			r_t.clear();

			old_seq_len = 0;
			min_cov_vertex = 0;
		}

		inline ResultCorrection reverseComplement() {

			if (seq.length() != 0){

				Roaring tmp;

				for (const uint32_t pos : pos_corrected_old_seq) tmp.add(old_seq_len - pos - 1);

				pos_corrected_old_seq = move(tmp);

				seq = reverse_complement(seq);

				reverse(qual.begin(), qual.end()); 
			}

			return *this;
		}

		inline void setSequence(const string& s) { seq = s; }
		inline void setSequence(string&& s) { seq = move(s); }

		inline void setQuality(const string& q) { qual = q; }
		inline void setQuality(string&& q) { qual = move(q); }

		inline const string& getSequence() const { return seq; }
		inline const string& getQuality() const { return qual; }

		inline size_t getLengthOldSequence() const { return old_seq_len; }

		inline void addCorrectedPosOldSeq(const uint32_t pos) { pos_corrected_old_seq.add(pos); }

		inline bool isCorrectedPosOldSeq(const uint32_t pos) const { return pos_corrected_old_seq.contains(pos); }

		inline void addCorrectedPosOldSeq(const uint64_t pos_start, const uint64_t pos_end) { pos_corrected_old_seq.addRange(pos_start, pos_end); }

		inline size_t getNbCorrectedPosOldSeq() const { return pos_corrected_old_seq.cardinality(); }

		inline bool isCorrected() const { return (pos_corrected_old_seq.cardinality() == old_seq_len); }

		inline size_t getLastCorrectedPosOldSeq() const { return pos_corrected_old_seq.maximum(); }

		inline void clearCorrectedPosOldSeq() { pos_corrected_old_seq = Roaring(); }

		inline size_t getLengthCorrectedRegion(const size_t pos) const {

			size_t sum = 0;

			for (size_t i = pos; i < old_seq_len; ++i){

				if (pos_corrected_old_seq.contains(i)) ++sum;
				else break;
			}

			return sum;
		}

		inline size_t getLengthUncorrectedRegion(const size_t pos) const {

			size_t sum = 0;

			for (size_t i = pos; i < old_seq_len; ++i){

				if (!pos_corrected_old_seq.contains(i)) ++sum;
				else break;
			}

			return sum;
		}

		inline void setSourcePairID(const PairID& p_id) { r_s = p_id; }
		inline void setTargetPairID(const PairID& p_id) { r_t = p_id; }
		inline void setWeakPairID(const PairID& p_id) { r_w = p_id; }
		inline void setMinCovVertex(const size_t v){ min_cov_vertex = v; }
		inline void setPartitions(const Roaring& r_part){ partitions = r_part; }

		inline const PairID& getSourcePairID() const { return r_s; }
		inline const PairID& getTargetPairID() const { return r_t; }
		inline const PairID& getWeakPairID() const { return r_w; }
		inline size_t getMinCovVertex() const { return min_cov_vertex; }
		inline const Roaring& getPartitions() const { return partitions; }

	//private:

		Roaring pos_corrected_old_seq;

		string seq;
		string qual;

		PairID r_s;
		PairID r_w;
		PairID r_t;

		Roaring partitions;

		size_t old_seq_len;
		size_t min_cov_vertex;
};

#endif