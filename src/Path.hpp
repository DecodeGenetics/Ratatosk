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

		Path(const const_UnitigMap<U>& um_start, const string& ext, const const_UnitigMap<U>& um_end) {

			l = 0; // Empty string

			if (um_start.isEmpty) return;

			const size_t k = um_start.getGraph()->getK();

			if (um_end.isEmpty && (ext.length() == 0)) {

				start = um_start;
				l = um_start.len + k - 1;

				return;
			}

			const_UnitigMap<U> curr = um_start;

			size_t len = um_start.len + k - 1;

			for (const char c : ext){

				const Kmer fw(curr.getMappedTail().forwardBase(c));

				curr = curr.getGraph()->find(fw, true);

				if (curr.isEmpty) return;

				curr.dist = 0;
				curr.len = curr.size - k + 1;

				len += curr.len;
			}

			if (!um_end.isEmpty) len += um_end.len;
			else return;

			start = um_start;
			end = um_end;

			succ = ext;

			l = len;
		}

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

		Path<U> rev_comp() const {

			Path<U> out;

			if (!start.isEmpty){

				out.l = l; // Length is the same

				// Quality score are the same but in reverse order
				{
					out.qual = qual;
					reverse(out.qual.begin(), out.qual.end());
				}

				if (end.isEmpty) {

					out.start = start;
					out.start.strand = !out.start.strand;
				}
				else {

					out.start = end;
					out.end = start;

					out.start.strand = !out.start.strand;
					out.end.strand = !out.end.strand;

					if (succ.length() != 0) {

						const size_t k = start.getGraph()->getK();

						const_UnitigMap<U> curr = start;

						out.succ.reserve(succ.length());

						for (const char c : succ){

							const Kmer fw(curr.getMappedTail().forwardBase(c));

							curr = curr.getGraph()->find(fw, true);

							curr.dist = 0;
							curr.len = curr.size - k + 1;

							if (curr.strand) out.succ.append(1, curr.getUnitigTail().toString()[0]);
							else out.succ.append(1, curr.getUnitigHead().twin().toString()[0]);
						}

						out.succ = reverse_complement(out.succ);
					}
				}
			}

			return out;
		}

		inline const const_UnitigMap<U>& front() const {

			return start;
		}

		inline const const_UnitigMap<U>& back() const {

			return (end.isEmpty ? start : end);
		}

		bool replace_back(const const_UnitigMap<U>& um, const string& qual_s) {

			if (um.isEmpty || start.isEmpty) return false;

			const size_t k = um.getGraph()->getK();
			const size_t um_subtr_len = um.len + k - 1;

			if (end.isEmpty){

				l = um_subtr_len;
				start = um;

				if (qual_s.length() == um_subtr_len) qual = qual_s;
				else return false;
			}
			else {

				l -= end.len;
				l += um.len;

				if (qual_s.length() == um_subtr_len) {

					const size_t end_substr_len = end.len + k - 1;

					qual.replace(qual.length() - end_substr_len, end_substr_len, qual_s, 0, qual_s.length());
				}
				else return false;

				end = um;
			}

			return true;
		}

		bool extend(const const_UnitigMap<U>& um) {

			if (um.isEmpty) return false;

			if (start.isEmpty){

				start = um;
				l = um.len + um.getGraph()->getK() - 1;
			}
			else {

				if (!end.isEmpty) {

					if (end.strand) succ.append(1, end.getUnitigHead().toString()[um.getGraph()->getK() - 1]);
					else succ.append(1, end.getUnitigTail().twin().toString()[um.getGraph()->getK() - 1]);
				}

				end = um;
				l += um.len;
			}

			return true;
		}

		bool extend(const const_UnitigMap<U>& um, const string& qual_s) {

			if (um.isEmpty) return false;

			const size_t k = um.getGraph()->getK();
			const size_t um_subtr_len = um.len + k - 1;

			if (start.isEmpty){

				start = um;
				l = um_subtr_len;

				if (qual_s.length() == um_subtr_len) qual = qual_s;
				else return false;
			}
			else {

				if (!end.isEmpty) {

					if (end.strand) succ.append(1, end.getUnitigHead().toString()[um.getGraph()->getK() - 1]);
					else succ.append(1, end.getUnitigTail().twin().toString()[um.getGraph()->getK() - 1]);
				}

				end = um;
				l += um.len;

				if (qual_s.length() == um_subtr_len) qual += qual_s.substr(k-1);
				else return false;
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

			//if (o.qual.length() != 0) qual.append(o.qual.substr(1));
			if (o.qual.length() != 0) qual.append(o.qual.substr(k));

            return true;
		}

		inline size_t size() const {

			return static_cast<size_t>(!start.isEmpty) + static_cast<size_t>(!end.isEmpty) + succ.length();
		}

		inline size_t length() const {

			return l;
		}

		/*inline void setQuality(const char c) {

			qual = string(size(), c);
		}*/

		inline void setQuality(const string& q) {

			if (q.length() == l) qual = q;
		}

		/*inline bool setQuality(const string& q) {

			if (q.length() == l) {

				if (q.length() != toString().length()) return false;

				qual = q;
			}
			else return false;

			return true;
		}*/

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

				path_str.append(mapped_str.substr(k-1));
			}

			if (!end.isEmpty){

				const string end_str = end.mappedSequenceToString();

				path_str.append(end_str.substr(k-1));
			}

			return path_str;
		}

		void prunePrefix(const size_t len) {

			if (start.isEmpty || (l == 0) || (len >= l)) return;

			const size_t k = start.getGraph()->getK();
			
			if (end.isEmpty) {

				if (!start.strand) start.dist += (l - len);

				start.len -= (l - len);
			}
			else if (succ.empty()) {

				if ((start.len + k - 1) >= len) {

					l = start.len + k - 1;
					end = const_UnitigMap<U>();

					if (!start.strand) start.dist += l - len;

					start.len -= l - len;
				}
				else {

					if (!end.strand) end.dist += l - len;

					end.len -= l - len;
				}
			}
			else {

				if ((start.len + k - 1) >= len) {

					l = start.len + k - 1;
					end = const_UnitigMap<U>();

					if (!start.strand) start.dist += l - len;

					start.len -= l - len;
				}
				else if (len > (l - end.len)) {

					if (!end.strand) end.dist += l - len;

					end.len -= l - len;
				}
				else {

					const_UnitigMap<U> curr_um = start;

					const string old_succ = move(succ);

					l = start.len + k - 1;

					for (const char c : old_succ) {

						const Kmer fw(curr_um.getMappedTail().forwardBase(c));

						curr_um = curr_um.getGraph()->find(fw, true);

						curr_um.dist = 0;
						curr_um.len = curr_um.size - k + 1;

						l += curr_um.len;

						if (l < len) succ += c;
						else {

							end = curr_um;

							if (!end.strand) end.dist += l - len;

							end.len -= l - len;

							break;
						}
					}
				}
			}

			l = len;

			if (qual.length() != 0) qual = qual.substr(0, l);
		}

		/*string toQualityString() const {

			if (qual.empty()) return string();
			if (start.isEmpty) return string();

			const size_t k = start.getGraph()->getK();

			if (end.isEmpty) return string(start.len + k - 1, qual[0]);

			const_UnitigMap<U> curr = start;

			size_t i = 1;

			string path_qual(start.len + k - 1, qual[0]);

			path_qual.reserve(l);

			for (const char c : succ){

				const Kmer fw(curr.getMappedTail().forwardBase(c));

				curr = curr.getGraph()->find(fw, true);

				curr.dist = 0;
				curr.len = curr.size - k + 1;

				path_qual.append(curr.len, qual[i++]);
			}

			if (!end.isEmpty) path_qual.append(end.len, qual[i++]);

			return path_qual;
		}*/

		inline const string& toQualityString() const {

			return qual;
		}

		PathOut toStringVector() const {

			if (start.isEmpty) return PathOut(*this);

			if (end.isEmpty){

				const string path_str = start.mappedSequenceToString();
				//const string path_qual = qual.empty() ? string() : string(path_str.length(), qual[0]);

				return PathOut(path_str, /*path_qual*/qual, vector<const_UnitigMap<U>>(1, start), *this);
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

				path_str.append(mapped_str.substr(k-1));
				v_um.push_back(curr);

				++i;
			}

			if (!end.isEmpty){

				const string end_str = end.mappedSequenceToString();

				path_str.append(end_str.substr(k-1));
				v_um.push_back(end);
			}

			/*if (!qual.empty()){

				path_qual = string(v_um[0].len + k - 1, qual[0]);
				path_qual.reserve(l);

				for (i = 1; i < v_um.size(); ++i) path_qual.append(v_um[i].len, qual[i]);
			}*/
			path_qual = qual;

			return PathOut(move(path_str), move(path_qual), move(v_um), *this);
		}

		PathOut toStringVector(const PathOut& po) const {

			if (start.isEmpty) return PathOut(*this);

			if (end.isEmpty){

				const string path_str = start.mappedSequenceToString();
				//const string path_qual = qual.empty() ? string() : string(path_str.length(), qual[0]);

				return PathOut(path_str, /*path_qual*/qual, vector<const_UnitigMap<U>>(1, start), *this);
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

				path_str.append(mapped_str.substr(k-1));
				v_um.push_back(curr);
			}

			if (!end.isEmpty){

				const string end_str = end.mappedSequenceToString();

				path_str.append(end_str.substr(k-1));
				v_um.push_back(end);
			}

			/*if (!qual.empty()){

				path_qual = string(v_um[0].len + k - 1, qual[0]);
				path_qual.reserve(l);

				for (size_t i = 1; i < v_um.size(); ++i) path_qual.append(v_um[i].len, qual[i]);
			}*/
			path_qual = qual;

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

		/*vector<char> toQualityVector() const {

			if (qual.empty() || start.isEmpty) return vector<char>();

			return vector<char>(qual.begin(), qual.end()); 
		}*/

		string getMiddleCompactedPath() const {

			return succ;
		}

	private:

		const_UnitigMap<U> start;
		const_UnitigMap<U> end;

		string succ;
		string qual;

		size_t l;
};

#endif
