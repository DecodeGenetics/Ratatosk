#ifndef RATATOSK_PATH_HPP
#define RATATOSK_PATH_HPP

#include <iostream>
#include <set>

#include "UnitigMap.hpp"

#include "Common.hpp"

using namespace std;

template<typename U> class Path;

template<typename Unitig_data_t = void>
class Path {

	typedef Unitig_data_t U;

	public:

		class PathOut {

			template<typename U> friend class Path;

			public:

				PathOut() : path(nullptr) {}

				PathOut(const PathOut& o) : path_seq_str(o.path_seq_str), path_qual_str(o.path_qual_str), v_path_um(o.v_path_um), path(o.path){}

				PathOut(PathOut&& o) : path_seq_str(move(o.path_seq_str)), path_qual_str(move(o.path_qual_str)), v_path_um(move(o.v_path_um)), path(o.path) {

					o.path = nullptr;
				}

			    PathOut& operator=(const PathOut& o) {

			        clear();

					path_seq_str = o.path_seq_str;
					path_qual_str = o.path_qual_str;
					v_path_um = o.v_path_um;
					path = o.path;

			        return *this;
			    }

			    PathOut& operator=(PathOut&& o) {

			    	if (this != &o) {

				        path_seq_str = move(o.path_seq_str);
				        path_qual_str = move(o.path_qual_str);
						v_path_um = move(o.v_path_um);
						path = o.path;

						o.path = nullptr;
					}

			        return *this;
			    }

			    void clear() {

			    	path_seq_str.clear();
			    	path_qual_str.clear();
			    	v_path_um.clear();

			    	path = nullptr;
			    }

			    inline const string& toString() const {

			    	return path_seq_str;
			    }

			    inline const string& toQualityString() const {

			    	return path_qual_str;
			    }

			    inline const vector<const_UnitigMap<U>>& toVector() const {

			    	return v_path_um;
			    }

			private:

				PathOut(const Path<U>& path_) : path(&path_) {}

				PathOut(const string& path_seq_str_, const string& path_qual_str_,
						const vector<const_UnitigMap<U>>& v_path_um_, const Path<U>& path_) :	path_seq_str(path_seq_str_), path_qual_str(path_qual_str_), v_path_um(v_path_um_), path(&path_) {}

				PathOut(string&& path_seq_str_, string&& path_qual_str_,
						vector<const_UnitigMap<U>>&& v_path_um_, const Path<U>& path_) : path_seq_str(move(path_seq_str_)), path_qual_str(move(path_qual_str_)), v_path_um(move(v_path_um_)), path(&path_) {}

				string path_seq_str;
				string path_qual_str;

				vector<const_UnitigMap<U>> v_path_um;
				
				const Path<U>* path; 
			
		};

		Path() : l(0) {}

		Path(const Path<U>& o) : start(o.start), end(o.end), succ(o.succ), qual(o.qual), l(o.l) {}

		Path(Path<U>&& o) : start(o.start), end(o.end), succ(move(o.succ)), qual(move(o.qual)), l(o.l) {

			o.clear();
		}

		~Path() {

			clear();
		}

	    Path& operator=(const Path<U>& o) {

	        clear();

			start = o.start;
			end = o.end;
			succ = o.succ;
			qual = o.qual;
			l = o.l;

	        return *this;
	    }

	    Path& operator=(Path<U>&& o) {

	    	if (this != &o) {

		        clear();

				start = o.start;
				end = o.end;
				succ = move(o.succ);
				qual = move(o.qual);
				l = o.l;

				o.clear();
			}

	        return *this;
	    }

		void clear() {

			start = const_UnitigMap<U>();
			end = const_UnitigMap<U>();

			succ.clear();
			qual.clear();

			l = 0;
		}

		inline const const_UnitigMap<U>& front() const { return start; }
		inline const const_UnitigMap<U>& back() const { return (end.isEmpty ? start : end); }

		bool replace_back(const const_UnitigMap<U>& um, const double score = -1.0) {

			if (um.isEmpty) return false;

			if (end.isEmpty){

				if ((score >= 0.0) && (score <= 1.0)){

					if (start.isEmpty) qual = string(1, getQual(score));
					else if (!qual.empty()) qual[0] = getQual(score);
				}

				l = um.len + um.getGraph()->getK() - 1;
				start = um;
			}
			else {

				l -= end.len;
				l += um.len;

				end = um;

				if ((score >= 0.0) && (score <= 1.0) && !qual.empty()) qual[qual.length()-1] = getQual(score);
			}

			return true;
		}

		bool extend(const const_UnitigMap<U>& um, const double score = -1.0) {

			if (um.isEmpty) return false;

			if (start.isEmpty){

				start = um;
				l = um.len + um.getGraph()->getK() - 1;

				if ((score >= 0.0) && (score <= 1.0)) qual += getQual(score);
			}
			else {

				if (!end.isEmpty) {

					if (end.strand) succ.append(1, end.getUnitigHead().toString()[um.getGraph()->getK() - 1]);
					else succ.append(1, end.getUnitigTail().twin().toString()[um.getGraph()->getK() - 1]);
				}

				end = um;
				l += um.len;

				if ((score >= 0.0) && (score <= 1.0) && !qual.empty()) qual += getQual(score);
			}

			return true;
		}

		// When merging P1 to P2, the last k-mer of P1 must be the same as the first kmer of P2
		bool merge(const Path<U>& o){

			if (o.length() == 0) return true;

			if (length() == 0){

				*this = o;
				return true;
			}

			// Make sure those two are merge-compatible
			if (qual.empty() != o.qual.empty()) return false;
			if (end.isEmpty && ((start.getUnitigHead() != o.start.getUnitigHead()) || (start.strand != o.start.strand))) return false;
			if (!end.isEmpty && ((end.getUnitigHead() != o.start.getUnitigHead()) || (end.strand != o.start.strand))) return false;

			const size_t k = start.getGraph()->getK();

			if (end.isEmpty){

				if (!start.strand) start.dist = o.start.dist;

				start.len += o.start.len - 1;
				end = o.end;
	            succ = o.succ;
			}
			else {

				if (!end.strand) end.dist = o.start.dist;

				end.len += o.start.len - 1;

				if ((o.succ.length() != 0) || !o.end.isEmpty) {

					if (end.strand) succ.append(1, end.getUnitigHead().toString()[k - 1]);
					else succ.append(1, end.getUnitigTail().twin().toString()[k - 1]);

	            	succ.append(o.succ);

	            	end = o.end;
				}
			}

			l += o.l - k;

			if (!qual.empty() && !o.qual.empty()) qual.append(o.qual.substr(1));

            return true;
		}

		inline size_t size() const {

			return (!start.isEmpty) + (!end.isEmpty) + succ.length();
		}

		inline size_t length() const {

			return l;
		}

		inline void setQuality(const char c) {

			qual = string(size(), c);
		}

		inline bool replaceQuality(const size_t pos_start, const size_t length, const char c) {

			if ((pos_start + length) <= qual.length()){

				qual.replace(pos_start, length, length, c);

				return true;
			}

			return false;
		}

		inline void rmQualityScore() {

			qual.clear();
		}

		string toString() const {

			if (start.isEmpty) return string();
			if (end.isEmpty) return start.mappedSequenceToString();

			const size_t k = start.getGraph()->getK();

			const_UnitigMap<U> curr = start;

			string path_str;

			path_str.reserve(l);
			path_str.append(start.mappedSequenceToString());

			for (const char c : succ){

				const Kmer fw(curr.getMappedTail().forwardBase(c));

				curr = curr.getGraph()->find(fw, true);

				curr.dist = 0;
				curr.len = curr.size - k + 1;

				const string mapped_str = curr.mappedSequenceToString();

				path_str.append(mapped_str.substr(k-1, mapped_str.length() - k + 1));
			}

			if (!end.isEmpty){

				const string end_str = end.mappedSequenceToString();

				path_str.append(end_str.substr(k-1, end_str.length() - k + 1));
			}

			return path_str;
		}

		string toQualityString() const {

			if (qual.empty()) return string();
			if (start.isEmpty) return string();

			const size_t k = start.getGraph()->getK();

			if (end.isEmpty) return string(start.len + k - 1, qual[0]);

			const_UnitigMap<U> curr = start;

			size_t i = 0;

			string path_qual(start.len + k - 1, qual[0]);

			path_qual.reserve(l);

			for (const char c : succ){

				const Kmer fw(curr.getMappedTail().forwardBase(c));

				curr = curr.getGraph()->find(fw, true);

				curr.dist = 0;
				curr.len = curr.size - k + 1;

				path_qual.append(string(curr.len, qual[++i]));
			}

			if (!end.isEmpty) path_qual.append(string(end.len, qual[++i]));

			return path_qual;
		}

		vector<char> toQualityVector() const {

			if (qual.empty() || start.isEmpty) return vector<char>();

			return vector<char>(qual.begin(), qual.end()); 
		}

		PathOut toStringVector() const {

			if (start.isEmpty) return PathOut(*this);

			if (end.isEmpty){

				const string path_str = start.mappedSequenceToString();
				const string path_qual = qual.empty() ? string() : string(path_str.length(), qual[0]);

				return PathOut(path_str, path_qual, vector<const_UnitigMap<U>>(1, start), *this);
			}

			const size_t k = start.getGraph()->getK();
			const string start_str = start.mappedSequenceToString();

			const_UnitigMap<U> curr = start;

			vector<const_UnitigMap<U>> v_um(1, start);

			string path_str;
			string path_qual;

			v_um.reserve(size());

			path_str.reserve(l);
			path_str.append(start_str);

			size_t i = 1;

			for (const char c : succ){

				const Kmer fw(curr.getMappedTail().forwardBase(c));

				curr = curr.getGraph()->find(fw, true);

				curr.dist = 0;
				curr.len = curr.size - k + 1;

				const string mapped_str = curr.mappedSequenceToString();

				path_str.append(mapped_str.substr(k-1, mapped_str.length() - k + 1));
				v_um.push_back(curr);

				++i;
			}

			if (!end.isEmpty){

				const string end_str = end.mappedSequenceToString();

				path_str.append(end_str.substr(k-1, end_str.length() - k + 1));
				v_um.push_back(end);
			}

			if (!qual.empty()){

				path_qual = string(v_um[0].len + k - 1, qual[0]);
				path_qual.reserve(l);

				for (i = 1; i < v_um.size(); ++i) path_qual.append(string(v_um[i].len, qual[i]));
			}

			return PathOut(move(path_str), move(path_qual), move(v_um), *this);
		}

		PathOut toStringVector(const PathOut& po) const {

			if (start.isEmpty) return PathOut(*this);

			if (end.isEmpty){

				const string path_str = start.mappedSequenceToString();
				const string path_qual = qual.empty() ? string() : string(path_str.length(), qual[0]);

				return PathOut(path_str, path_qual, vector<const_UnitigMap<U>>(1, start), *this);
			}

			size_t len_pref_vertices = 0;
			size_t len_pref_str = 0;

			const size_t k = start.getGraph()->getK();
			const size_t succ_len = succ.length();

			const Path<U>& o = *(po.path);

			// Compute the Longest Common Prefix of the two paths in terms of vertices
			if (!o.start.isEmpty && (start == o.start)){

				++len_pref_vertices;

				if (!succ.empty() && !o.succ.empty()) len_pref_vertices += cstrMatch(succ.c_str(), o.succ.c_str());
			}

			const bool no_prefix = (len_pref_vertices == 0);

			const_UnitigMap<U> curr = no_prefix ? start : po.v_path_um[len_pref_vertices - 1];

			vector<const_UnitigMap<U>> v_um(po.v_path_um.begin(), po.v_path_um.begin() + len_pref_vertices);

			string path_str;
			string path_qual;

			path_str.reserve(l);
			v_um.reserve(size());
			
			if (no_prefix){

				path_str.append(start.mappedSequenceToString());
				v_um.push_back(start);
			}
			else {

				len_pref_str = po.v_path_um[0].len + k - 1;

				for (size_t i = 1; i < len_pref_vertices; ++i) len_pref_str += po.v_path_um[i].len;

				path_str.append(po.path_seq_str.substr(0, len_pref_str));
			}

			for (size_t i = no_prefix ? 0 : len_pref_vertices - 1; i < succ_len; ++i){

				const char c = succ[i];
				const Kmer fw(curr.getMappedTail().forwardBase(c));

				curr = curr.getGraph()->find(fw, true);

				curr.dist = 0;
				curr.len = curr.size - k + 1;

				const string mapped_str = curr.mappedSequenceToString();

				path_str.append(mapped_str.substr(k-1, mapped_str.length() - k + 1));
				v_um.push_back(curr);
			}

			if (!end.isEmpty){

				const string end_str = end.mappedSequenceToString();

				path_str.append(end_str.substr(k-1, end_str.length() - k + 1));
				v_um.push_back(end);
			}

			if (!qual.empty()){

				path_qual = string(v_um[0].len + k - 1, qual[0]);
				path_qual.reserve(l);

				for (size_t i = 1; i < v_um.size(); ++i) path_qual.append(string(v_um[i].len, qual[i]));
			}

			return PathOut(move(path_str), move(path_qual), move(v_um), *this);
		}

		vector<const_UnitigMap<U>> toVector() const {

			if (start.isEmpty) return vector<const_UnitigMap<U>>();
			if (end.isEmpty) return vector<const_UnitigMap<U>>(1, start);

			const size_t k = start.getGraph()->getK();

			const_UnitigMap<U> curr = start;

			vector<const_UnitigMap<U>> v_um(1, start);

			v_um.reserve(size());

			for (const char c : succ){

				const Kmer fw(curr.getMappedTail().forwardBase(c));

				curr = curr.getGraph()->find(fw, true);

				curr.dist = 0;
				curr.len = curr.size - k + 1;

				v_um.push_back(curr);
			}

			if (!end.isEmpty) v_um.push_back(end);

			return v_um;
		}

	private:

		const_UnitigMap<U> start;
		const_UnitigMap<U> end;

		string succ;
		string qual;

		size_t l;
};

#endif
