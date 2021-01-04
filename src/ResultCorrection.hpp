#ifndef RATATOSK_RESULTCORRECTION_HPP
#define RATATOSK_RESULTCORRECTION_HPP

#include <iostream>

class ResultCorrection {

	public: 

		ResultCorrection(const size_t seq_len) : old_seq_len(seq_len), is_corrected(false) {}

		ResultCorrection(const ResultCorrection& o) : 	pos_corrected_old_seq(o.pos_corrected_old_seq), seq(o.seq), qual(o.qual), old_seq_len(o.old_seq_len),
														is_corrected(o.is_corrected), w_pids(o.w_pids) {}

		ResultCorrection(ResultCorrection&& o) :	pos_corrected_old_seq(move(o.pos_corrected_old_seq)), seq(move(o.seq)), qual(move(o.qual)), old_seq_len(o.old_seq_len),
													is_corrected(o.is_corrected), w_pids(move(o.w_pids)) { o.clear(); }

		ResultCorrection operator=(const ResultCorrection& o){

			if (&o != this){

				pos_corrected_old_seq = o.pos_corrected_old_seq;
				old_seq_len = o.old_seq_len;

				is_corrected = o.is_corrected;

				seq = o.seq;
				qual = o.qual;
				w_pids = o.w_pids;
			}

			return *this;
		}

		ResultCorrection operator=(ResultCorrection&& o){

			if (&o != this){

				pos_corrected_old_seq = move(o.pos_corrected_old_seq);
				old_seq_len = o.old_seq_len;

				is_corrected = o.is_corrected;

				seq = move(o.seq);
				qual = move(o.qual);
				w_pids = move(o.w_pids);

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
			w_pids.clear();

			old_seq_len = 0;

			is_corrected = false;
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

		inline bool isCorrected() const { return is_corrected; }

		inline void setCorrected() { is_corrected = true; }

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

		inline double getMeanQualityScore(const size_t start, const size_t end) const {

			double score = 0.0;

			if (start < end) {

				for (size_t pos = start; pos < end; ++pos) score += getScore(qual[pos]);

				score /= static_cast<double>(end - start);
			}

			return score;
		}

		inline const Roaring& getPartitions() const { return partitions; }

		inline void setWeightsPairID(const WeightsPairID& w_pids_) { w_pids = w_pids_; }
		inline const WeightsPairID& getWeightsPairID() const { return w_pids; }

	//private:

		Roaring pos_corrected_old_seq;

		string seq;
		string qual;

		WeightsPairID w_pids;

		Roaring partitions;

		size_t old_seq_len;

		bool is_corrected;
};

#endif