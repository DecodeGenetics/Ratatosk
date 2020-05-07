#include "GraphTraversal.hpp"

pair<vector<Path<UnitigData>>, bool> explorePathsBFS(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const TinyBloomFilter<uint32_t>& bf, const PairID& r,
														const const_UnitigMap<UnitigData>& um_s,
														const size_t min_cov_vertex, const bool long_read_correct){

	auto resizeVector = [&](vector<Path<UnitigData>>& v) {

		if (v.size() <= 1) return;

		vector<const Path<UnitigData>*> v_ptr;

		for (const auto& p : v) v_ptr.push_back(&p);

		const pair<int, int> p = selectBestPrefixAlignment(ref, ref_len, v_ptr);

		v = vector<Path<UnitigData>>(1, v[p.first]);
	};

	auto resizeQueue = [&](queue<Path<UnitigData>>& q) {

		if (q.empty()) return;

		vector<Path<UnitigData>> v;

	    while (!q.empty()){

	    	v.push_back(move(q.front()));
	    	q.pop();
	    }

	    resizeVector(v);

	    for (auto& p : v) q.push(move(p));
	};

	auto computeScore = [&](const Path<UnitigData>& p) {

		const Path<UnitigData>::PathOut p_str_vum = p.toStringVector();
		const string& path_curr = p_str_vum.toString();
    	const vector<const_UnitigMap<UnitigData>>& v_um = p_str_vum.toVector();

		EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		EdlibAlignResult align = edlibAlign(path_curr.c_str(), path_curr.length(), ref, ref_len, config);

    	PairID r_um = v_um[0].getData()->get_readID();

    	for (size_t j = 1; j < v_um.size(); ++j){

    		if (!v_um[j].isSameReferenceUnitig(v_um[j-1])) r_um |= v_um[j].getData()->get_readID();
    	}

		const double score_kmer = 1.0 - (static_cast<double>(align.editDistance) / path_curr.length());
		const double score_read = static_cast<double>(getNumberSharedPairID(r, r_um)) / r.cardinality();
		const double score = (score_read * score_kmer) / ((score_read * score_kmer) + (1.0 - score_read) * (1.0 - score_kmer));

		edlibFreeAlignResult(align);

		return (std::isnan(score) ? 0.0 : score);
	};

	auto explore = [&](const const_UnitigMap<UnitigData>& um, const size_t level) -> vector<Path<UnitigData>> {

		vector<Path<UnitigData>> best_successors;
		vector<Path<UnitigData>> over_represented;

		double best_score = 0;

		if (long_read_correct){

			exploreSubGraphLong(opt, bf, r, ref, ref_len, um, const_UnitigMap<UnitigData>(), Path<UnitigData>(), um.getData()->get_readID(),
								min_cov_vertex, best_score, best_successors, over_represented);
		}
		else {

			exploreSubGraph(opt, bf, r, ref, ref_len, um, const_UnitigMap<UnitigData>(), Path<UnitigData>(), um.getData()->get_readID(),
							level - 1, min_cov_vertex, best_score, best_successors, over_represented);
		}

		if (!best_successors.empty()) best_successors = vector<Path<UnitigData>>(1, best_successors[selectBestSubstringAlignment(ref, ref_len, best_successors).first]);

		return best_successors;
	};

	vector<Path<UnitigData>> v, v_tmp;

	if (!um_s.isEmpty){

		const size_t level = 3;

		const double seq_entropy = getEntropy(ref, ref_len);

		const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

		const size_t max_len_path = max((ref_len - opt.k) * opt.weak_region_len_factor, 1.0) + opt.k;
		const size_t min_len_path = max((ref_len - opt.k) / opt.weak_region_len_factor, 1.0) + opt.k;

		const size_t max_paths = getMaxPaths(seq_entropy, max_len_path, opt.k);

		size_t max_branches = 2 * getMaxBranch(seq_entropy, max_len_path, opt.k);
		size_t max_sz_stck = 512;

		const size_t dec_sz_stck = max_sz_stck / (((max_len_path - opt.k + 1) / max_branches) + 1);

		KmerHashTable<vector<Path<UnitigData>>> km_h;

		queue<Path<UnitigData>> q, q_long;
	    
	    const_UnitigMap<UnitigData> um_start_tmp(um_s); // Create a non-const local copy of unitig given in parameter
	    Path<UnitigData> p_tmp;

	    auto start_time = std::chrono::steady_clock::now();

		if (um_start_tmp.strand){

			um_start_tmp.dist += um_start_tmp.len - 1;
			um_start_tmp.len = um_start_tmp.size - um_start_tmp.dist - opt.k + 1;
		}
	    else {

	    	um_start_tmp.len = um_s.dist + 1;
	    	um_start_tmp.dist = 0;
	    }

		if ((um_start_tmp.len + opt.k - 1) >= min_len_path) {

			const_UnitigMap<UnitigData> back(um_start_tmp);

			if ((back.len + opt.k - 1) > max_len_path){

	    		if (!back.strand) back.dist = back.len - (max_len_path - opt.k + 1);

			    back.len = max_len_path - opt.k + 1;
			}

		    if (out_qual) p_tmp.extend(back, 1.0);
		    else p_tmp.extend(back);

    		v.push_back(p_tmp);
    		p_tmp.clear();
		}

	    if (out_qual) p_tmp.extend(um_start_tmp, 1.0);
	    else p_tmp.extend(um_start_tmp);

	    q.push(move(p_tmp));

	    while (v.empty()){

		    while (!q.empty()){

		    	if (q.size() >= max_sz_stck) resizeQueue(q);

		        Path<UnitigData> p = move(q.front()); // Get the unitig on top of the stack
		        const const_UnitigMap<UnitigData>& um = p.back();

		        q.pop(); // Delete unitig on the top of the stack

		        if (p.length() < max_len_path){

		        	if (p.size() < max_branches){

						const Kmer head = um.strand ? um.getUnitigHead() : um.getUnitigHead().twin();
		        		const KmerHashTable<vector<Path<UnitigData>>>::const_iterator it_km_h = km_h.find(head);

		        		vector<Path<UnitigData>> exploreFW;

		        		if (it_km_h == km_h.end()){

		        			exploreFW = explore(um, level);
		        			km_h.insert(head, exploreFW);
		        		}
		        		else exploreFW = *it_km_h;

		        		bool extended = false;

			        	for (const auto& path : exploreFW){

			        		Path<UnitigData> p_tmp(p);

			        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();

			        		if (out_qual){

				        		const vector<char> v_qual = path.toQualityVector();

				        		for (size_t i = 0; i < v_path.size(); ++i){

				        			p_tmp.extend(v_path[i], getScore(v_qual[i]));

				        			if ((p_tmp.length() >= min_len_path) && (p_tmp.length() <= max_len_path)){

				        				v_tmp.push_back(p_tmp);
				        				extended = true;
				        			}
				        		}
			        		}
			        		else {

				        		for (size_t i = 0; i < v_path.size(); ++i){

				        			p_tmp.extend(v_path[i]);

				        			if ((p_tmp.length() >= min_len_path) && (p_tmp.length() <= max_len_path)){

				        				v_tmp.push_back(p_tmp);
				        				extended = true;
				        			}
				        		}
			        		}

			        		if ((!long_read_correct && (path.size() == level)) || (long_read_correct && (path.length() >= opt.large_k))){

				        		q.push(move(p_tmp));
				        		extended = true;
				        	}
		        		}

		        		if (!extended && (p.length() <= max_len_path) && (computeScore(p) >= 0.25)) v_tmp.push_back(p);

		        		if (v_tmp.size() >= max_paths){

		        			for (auto& p_tmp : v_tmp){

								const_UnitigMap<UnitigData> back = p_tmp.back();

								if (p_tmp.length() > max_len_path){

						    		if (!back.strand) back.dist = back.len - (max_len_path - opt.k + 1);
								    
								    back.len = max_len_path - opt.k + 1;

								    p_tmp.replace_back(back);
								}

					    		v.push_back(move(p_tmp));
							}

							v_tmp.clear();
						}
					}
					else {

	        			q_long.push(p);

	        			if (q_long.size() >= max_sz_stck) resizeQueue(q_long);
	        		}
		        }
		    }

			for (auto& p_tmp : v_tmp){

				const_UnitigMap<UnitigData> back = p_tmp.back();

				if (p_tmp.length() > max_len_path){

		    		if (!back.strand) back.dist = back.len - (max_len_path - opt.k + 1);
				    
				    back.len = max_len_path - opt.k + 1;

				    p_tmp.replace_back(back);
				}

	    		v.push_back(move(p_tmp));
			}

			if (v.empty() && !q_long.empty()){

				max_branches += max_branches;
				max_sz_stck = max(max_sz_stck - dec_sz_stck, static_cast<size_t>(2));

				if (max_branches >= (max_len_path - opt.k + 1)) break;

				resizeQueue(q_long);

				q = move(q_long);

            	if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_time).count() >= opt.max_time) break;
			}
			else break;
		}
	}

	return {v, true};
}

pair<vector<Path<UnitigData>>, bool> explorePathsBFS2(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const TinyBloomFilter<uint32_t>& bf, const PairID& r,
														const const_UnitigMap<UnitigData>& um_s, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_e,
														const size_t min_cov_vertex, const bool long_read_correct){

	KmerHashTable<const_UnitigMap<UnitigData>> um_h;

	auto customSort = [](const pair<size_t, size_t>& p1, const pair<size_t, size_t>& p2) {

		return (p1.first > p2.first);
	};

	auto resizeVector = [&](vector<Path<UnitigData>>& v) { // Select 2 "best" paths based on alignment and colors

		if (v.size() <= 1) return;

		vector<const Path<UnitigData>*> v_ptr;

		for (const auto& p : v) v_ptr.push_back(&p);

		const pair<int, int> p = selectBestPrefixAlignment(ref, ref_len, v_ptr);

		v = vector<Path<UnitigData>>(1, v[p.first]);
	};

	auto resizeQueue = [&](queue<Path<UnitigData>>& q) { // Select 2 "best" paths based on alignment and colors

		if (q.empty()) return;

		vector<Path<UnitigData>> v;

	    while (!q.empty()){

	    	v.push_back(move(q.front()));
	    	q.pop();
	    }

	    resizeVector(v);

	    for (auto& p : v) q.push(move(p));
	};

	auto explore = [&](const const_UnitigMap<UnitigData>& um, const size_t level) -> pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> {

		vector<Path<UnitigData>> best_successors;
		vector<Path<UnitigData>> over_represented;

		double best_score = 0;

		if (long_read_correct){

			exploreSubGraphLong(opt, bf, r, ref, ref_len, um, um_h, Path<UnitigData>(), um.getData()->get_readID(),
								min_cov_vertex, best_score, best_successors, over_represented);
		}
		else {

			exploreSubGraph(opt, bf, r, ref, ref_len, um, um_h, Path<UnitigData>(), um.getData()->get_readID(),
							level - 1, min_cov_vertex, best_score, best_successors, over_represented);
		}

		if (!best_successors.empty()) best_successors = vector<Path<UnitigData>>(1, best_successors[selectBestSubstringAlignment(ref, ref_len, best_successors).first]);

    	return {over_represented, best_successors};
	};

	vector<Path<UnitigData>> v, v_tmp;

	if (!um_s.isEmpty){

		const size_t level = 3;

		const double seq_entropy = getEntropy(ref, ref_len);

		const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

		const size_t max_len_path = max((ref_len - opt.k) * opt.weak_region_len_factor, 1.0) + opt.k;
		const size_t min_len_path = max((ref_len - opt.k) / opt.weak_region_len_factor, 1.0) + opt.k;

		const size_t max_paths = getMaxPaths(seq_entropy, max_len_path, opt.k);

		size_t max_branches = 2 * getMaxBranch(seq_entropy, max_len_path, opt.k);
		size_t max_sz_stck = 512;

		const size_t dec_sz_stck = max_sz_stck / (((max_len_path - opt.k + 1) / max_branches) + 1);

		KmerHashTable<pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>>> km_h;

		queue<Path<UnitigData>> q, q_long;
	    
	    const_UnitigMap<UnitigData> um_start_tmp(um_s); // Create a non-const local copy of unitig given in parameter
	    Path<UnitigData> p_start_tmp;

	    auto start_time = std::chrono::steady_clock::now();

	    for (const auto& p : v_um_e){

	    	um_h.insert(p.second.strand ? p.second.getUnitigHead() : p.second.getUnitigHead().twin(), p.second);
	    }

		if (um_start_tmp.strand){

			um_start_tmp.dist += um_start_tmp.len - 1;
			um_start_tmp.len = um_start_tmp.size - um_start_tmp.dist - opt.k + 1;
		}
	    else {

	    	um_start_tmp.len = um_s.dist + 1;
	    	um_start_tmp.dist = 0;
	    }

		const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it = um_h.find(um_start_tmp.strand ? um_start_tmp.getUnitigHead() : um_start_tmp.getUnitigHead().twin());

		if (it != um_h.end()){

			const const_UnitigMap<UnitigData>& um_e = *it;

			if (um_start_tmp.dist <= um_e.dist){

				const size_t len = (um_start_tmp.len + opt.k - 1) - (um_e.strand ? um_e.size - um_e.dist - opt.k : um_e.dist);

				if ((len >= min_len_path) && (len <= max_len_path)) {

		    		const_UnitigMap<UnitigData> back_tmp(um_start_tmp);

		    		if (back_tmp.strand) back_tmp.len = um_e.dist - back_tmp.dist + 1;
				    else {

				    	back_tmp.dist = um_e.dist;
				    	back_tmp.len -= um_e.dist;
				    }

				    if (out_qual) p_start_tmp.extend(back_tmp, 1.0);
				    else p_start_tmp.extend(back_tmp);

	        		v.push_back(p_start_tmp);
	        		p_start_tmp.clear();
				}
			}
		}

		if (out_qual) p_start_tmp.extend(um_start_tmp, 1.0);
		else p_start_tmp.extend(um_start_tmp);

	    q.push(move(p_start_tmp));

	    while (v.empty()){

		    while (!q.empty()){

		    	if (q.size() >= max_sz_stck) resizeQueue(q);

		        Path<UnitigData> p = move(q.front()); // Get the unitig on top of the stack
		        const const_UnitigMap<UnitigData>& um = p.back();
		        const Kmer head = um.strand ? um.getUnitigHead() : um.getUnitigHead().twin();

		        q.pop(); // Delete unitig on the top of the stack

		        if (p.length() < max_len_path){

		        	if (p.size() < max_branches) {

						pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> exploreFW;

		        		const KmerHashTable<pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>>>::const_iterator it_km_h = km_h.find(head);

		        		if (it_km_h == km_h.end()){

		        			exploreFW = explore(um, level);

		        			km_h.insert(head, exploreFW);
		        		}
		        		else exploreFW = *it_km_h;

			        	for (const auto& path : exploreFW.first){

			        		Path<UnitigData> p_tmp(p);

			        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();

			        		if (out_qual){

				        		const vector<char> v_qual = path.toQualityVector();

				        		for (size_t i = 0; i < v_path.size(); ++i) p_tmp.extend(v_path[i], getScore(v_qual[i]));
			        		}
			        		else {

			        			for (size_t i = 0; i < v_path.size(); ++i) p_tmp.extend(v_path[i]);
			        		}

			        		v_tmp.push_back(move(p_tmp));
		        		}

			        	for (const auto& path : exploreFW.second){

			        		if ((!long_read_correct && (path.size() == level)) || (long_read_correct && (path.length() >= opt.large_k))){

				        		Path<UnitigData> p_tmp(p);

				        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();

				        		if (out_qual){

					        		const vector<char> v_qual = path.toQualityVector();

					        		for (size_t i = 0; i < v_path.size(); ++i) p_tmp.extend(v_path[i], getScore(v_qual[i]));
				        		}
				        		else {

				        			for (size_t i = 0; i < v_path.size(); ++i) p_tmp.extend(v_path[i]);
				        		}

			        			q.push(move(p_tmp));
			        		}
		        		}
		        		
		        		if (v_tmp.size() >= max_paths){

							for (auto& p_tmp : v_tmp){

								const const_UnitigMap<UnitigData>& back = p_tmp.back();
								const Kmer head_back = back.strand ? back.getUnitigHead() : back.getUnitigHead().twin();
								const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it_um_h = um_h.find(head_back);

								if (it_um_h != um_h.end()){

									const const_UnitigMap<UnitigData>& um_e = *it_um_h;

						    		if (back.dist <= um_e.dist){

						    			const size_t len = p_tmp.length() - (um_e.strand ? um_e.size - um_e.dist - opt.k : um_e.dist);

						    			if ((len >= min_len_path) && (len <= max_len_path)) {

								    		const_UnitigMap<UnitigData> back_tmp(back);

								    		if (back_tmp.strand) back_tmp.len = um_e.dist - back_tmp.dist + 1;
										    else {

										    	back_tmp.dist = um_e.dist;
										    	back_tmp.len -= um_e.dist;
										    }

										    if (v.size() + 1 >= max_paths) resizeVector(v);

							        		v.push_back(move(p_tmp));
							        		v.back().replace_back(back_tmp);
						    			}
						    		}
					    		}
							}

		        			v_tmp.clear();
		        		}
	        		}
	        		else if (max_branches < (max_len_path - opt.k + 1)){

	        			q_long.push(p);

	        			if (q_long.size() >= max_sz_stck) resizeQueue(q_long);
	        		}
		        }

		        const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it_um_h = um_h.find(head);

		        if (it_um_h != um_h.end()){

		        	const const_UnitigMap<UnitigData>& um_e = *it_um_h;

		    		if (um.dist <= um_e.dist){

		    			const size_t len = p.length() - (um_e.strand ? um_e.size - um_e.dist - opt.k : um_e.dist);

		    			if ((len >= min_len_path) && (len <= max_len_path)) {

				    		const_UnitigMap<UnitigData> back_tmp(um);

				    		if (back_tmp.strand) back_tmp.len = um_e.dist - back_tmp.dist + 1;
						    else {

						    	back_tmp.dist = um_e.dist;
						    	back_tmp.len -= um_e.dist;
						    }

						    if (v.size() + 1 >= max_paths) resizeVector(v);

			        		v.push_back(p);
			        		v.back().replace_back(back_tmp);
		    			}
		    		}
		    	}
		    }

			for (auto& p_tmp : v_tmp){

				const const_UnitigMap<UnitigData>& back = p_tmp.back();
				const Kmer head_back = back.strand ? back.getUnitigHead() : back.getUnitigHead().twin();
				const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it_um_h = um_h.find(head_back);

				if (it_um_h != um_h.end()){

					const const_UnitigMap<UnitigData>& um_e = *it_um_h;

		    		if (back.dist <= um_e.dist){

		    			const size_t len = p_tmp.length() - (um_e.strand ? um_e.size - um_e.dist - opt.k : um_e.dist);

		    			if ((len >= min_len_path) && (len <= max_len_path)) {

				    		const_UnitigMap<UnitigData> back_tmp(back);

				    		if (back_tmp.strand) back_tmp.len = um_e.dist - back_tmp.dist + 1;
						    else {

						    	back_tmp.dist = um_e.dist;
						    	back_tmp.len -= um_e.dist;
						    }

						    if (v.size() + 1 >= max_paths) resizeVector(v);

			        		v.push_back(move(p_tmp));
			        		v.back().replace_back(back_tmp);
		    			}
		    		}
	    		}
			}

			v_tmp.clear();

			if (v.empty() && !q_long.empty()){

				max_branches += max_branches;
				max_sz_stck = max(max_sz_stck - dec_sz_stck, static_cast<size_t>(2));

				if (max_branches >= (max_len_path - opt.k + 1)) break;

				resizeQueue(q_long);

				q = move(q_long);

				if (std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_time).count() >= opt.max_time) break;
			}
			else break;
		}
	}

	return {v, true};
}

void exploreSubGraph(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
					const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e,
					const Path<UnitigData>& p, const PairID& r, const size_t level, const size_t min_cov_vertex,
					double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& over_represented) {

	struct info_traversal {

		Path<UnitigData> p;
		PairID r;

		size_t l;

		info_traversal(const Path<UnitigData>& p_, const PairID& r_, const size_t l_) : p(p_), r(r_), l(l_) {};
	};

	auto checkContiguity = [](const vector<const_UnitigMap<UnitigData>>& v_um, const size_t min_cov){

		PairID p_inter;

		for (const auto& um : v_um){

			const PairID& p_id_um = um.getData()->get_readID();

			if (p_id_um.cardinality() == 0) return true;
			else p_inter &= p_id_um;
		}

		if (p_inter.cardinality() < min_cov) return false;

		return true;	
	};

	auto getScorePath = [&checkContiguity, &p_read, &ref, ref_len](const Path<UnitigData>& path, const PairID& p_id_path){

		const Path<UnitigData>::PathOut path_out = path.toStringVector();

		if ((path.length() <= 50) && !checkContiguity(path_out.toVector(), 3)) return std::numeric_limits<double>::min();

		const string& path_curr = path_out.toString();

		const EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

		EdlibAlignResult align = edlibAlign(path_curr.c_str(), path_curr.length(), ref, ref_len, config);

		const double score_kmer = 1.0 - (static_cast<double>(align.editDistance) / path_curr.length());
		const double score_read = static_cast<double>(getNumberSharedPairID(p_read, p_id_path)) / p_read.cardinality();
		const double score = (score_read * score_kmer) / ((score_read * score_kmer) + (1.0 - score_read) * (1.0 - score_kmer));

		edlibFreeAlignResult(align);

		return (std::isnan(score) ? 0.0 : score);
	};

	const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

	stack<info_traversal> stck;

	stck.push(info_traversal(p, r, level));

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* pli = um_start.getData();

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			//if (pli->get_count_successor(um_start, succ) >= opt.min_cov_edges) {

				const UnitigData* pli_succ = succ.getData();
				const PairID& r_succ = pli_succ->get_readID();

				if (!um_e.isEmpty && (succ.isSameReferenceUnitig(um_e)) && (um_e.strand == succ.strand) && (succ.dist <= um_e.dist)){

					Path<UnitigData> path(i_t.p);

					path.extend(succ);

					if (!out_qual){

						if (path.length() <= 50) {

							if (checkContiguity(path.toVector(), 3)) over_represented.push_back(move(path));
						}
						else over_represented.push_back(move(path));
					}
					else {

						const PairID r_curr = i_t.r | r_succ;
						const double score = getScorePath(path, r_curr);

						if (score != std::numeric_limits<double>::min()){

							path.setQuality(getQual(score, opt.out_qual));

							over_represented.push_back(move(path));
						}
					}
				}

				size_t l_count = 0;

				if (r_succ.isEmpty()) l_count = min_cov_vertex;
				else l_count = min_cov_vertex & (static_cast<size_t>(!hasEnoughSharedPairID(bf_read, p_read, r_succ, min_cov_vertex)) - 1);

				if (l_count >= min_cov_vertex){

					Path<UnitigData> path(i_t.p);

					path.extend(succ);

					const PairID r_curr = i_t.r | r_succ;
					
					if (i_t.l != 0) stck.push(info_traversal(path, r_curr, i_t.l - 1));
					else {

						const double score = getScorePath(path, r_curr);

						if (score > best_score){

							best_score = score;
							best_successors.clear();
						}
						
						if (score >= best_score){

							if (out_qual) path.setQuality(getQual(score, opt.out_qual));

							best_successors.push_back(move(path));
						}
					}
				}
			//}
		}
	}
}

void exploreSubGraph(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
					const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e,
					const Path<UnitigData>& p, const PairID& r, const size_t level, const size_t min_cov_vertex,
					double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& over_represented) {

	struct info_traversal {

		Path<UnitigData> p;
		PairID r;

		size_t l;

		info_traversal(const Path<UnitigData>& p_, const PairID& r_, const size_t l_) : p(p_), r(r_), l(l_) {};
	};

	auto checkContiguity = [](const vector<const_UnitigMap<UnitigData>>& v_um, const size_t min_cov){

		PairID p_inter;

		for (const auto& um : v_um){

			const PairID& p_id_um = um.getData()->get_readID();

			if (p_id_um.cardinality() == 0) return true;
			else p_inter &= p_id_um;
		}

		if (p_inter.cardinality() < min_cov) return false;

		return true;	
	};

	auto getScorePath = [&checkContiguity, &p_read, &ref, ref_len](const Path<UnitigData>& path, const PairID& p_id_path){

		const Path<UnitigData>::PathOut path_out = path.toStringVector();

		if ((path.length() <= 50) && !checkContiguity(path_out.toVector(), 3)) return std::numeric_limits<double>::min();

		const string& path_curr = path_out.toString();

		const EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

		EdlibAlignResult align = edlibAlign(path_curr.c_str(), path_curr.length(), ref, ref_len, config);

		const double score_kmer = 1.0 - (static_cast<double>(align.editDistance) / path_curr.length());
		const double score_read = static_cast<double>(getNumberSharedPairID(p_read, p_id_path)) / p_read.cardinality();
		const double score = (score_read * score_kmer) / ((score_read * score_kmer) + (1.0 - score_read) * (1.0 - score_kmer));

		edlibFreeAlignResult(align);

		return (std::isnan(score) ? 0.0 : score);
	};

	const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

	stack<info_traversal> stck;

	stck.push(info_traversal(p, r, level));

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* pli = um_start.getData();

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			//if (pli->get_count_successor(um_start, succ) >= opt.min_cov_edges) {

				const UnitigData* pli_succ = succ.getData();
				const PairID& r_succ = pli_succ->get_readID();
				const Kmer head = succ.strand ? succ.getUnitigHead() : succ.getUnitigHead().twin();
				const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it_h_um_e = h_um_e.find(head);

				if (it_h_um_e != h_um_e.end()){

					const const_UnitigMap<UnitigData>& um_e = *it_h_um_e;

					if (succ.dist <= um_e.dist){

						Path<UnitigData> path(i_t.p);

						path.extend(succ);

						if (!out_qual){

							if (path.length() <= 50) {

								if (checkContiguity(path.toVector(), 3)) over_represented.push_back(move(path));
							}
							else over_represented.push_back(move(path));
						}
						else {

							const PairID r_curr = i_t.r | r_succ;
							const double score = getScorePath(path, r_curr);

							if (score != std::numeric_limits<double>::min()){

								path.setQuality(getQual(score, opt.out_qual));

								over_represented.push_back(move(path));
							}
						}
					}
				}

				size_t l_count = 0;

				if (r_succ.isEmpty()) l_count = min_cov_vertex;
				else l_count = min_cov_vertex & (static_cast<size_t>(!hasEnoughSharedPairID(bf_read, p_read, r_succ, min_cov_vertex)) - 1);

				if (l_count >= min_cov_vertex){

					Path<UnitigData> path(i_t.p);

					path.extend(succ);

					const PairID r_curr = i_t.r | r_succ;
					
					if (i_t.l != 0) stck.push(info_traversal(path, r_curr, i_t.l - 1));
					else {

						const double score = getScorePath(path, r_curr);

						if (score > best_score){

							best_score = score;
							best_successors.clear();
						}
						
						if (score >= best_score){

							if (out_qual) path.setQuality(getQual(score, opt.out_qual));

							best_successors.push_back(move(path));
						}
					}
				}
			//}
		}
	}
}

void exploreSubGraphLong(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
						const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e,
						const Path<UnitigData>& p, const PairID& r, const size_t min_cov_vertex,
						double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& semi_weak_successors) {

	struct info_traversal {

		Path<UnitigData> p;
		PairID r;

		bool c;

		info_traversal(const Path<UnitigData>& p_, const PairID& r_, const bool c_) : p(p_), r(r_), c(c_) {};
	};

	auto getScorePath = [&p_read, &ref, ref_len](const Path<UnitigData>& path, const PairID& p_id_path){

		const string path_curr = path.toString();

		EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		EdlibAlignResult align = edlibAlign(path_curr.c_str(), path_curr.length(), ref, ref_len, config);

		const double score_kmer = 1.0 - (static_cast<double>(align.editDistance) / path_curr.length());
		const double score_read = static_cast<double>(getNumberSharedPairID(p_read, p_id_path)) / p_read.cardinality();
		const double score = (score_read * score_kmer) / ((score_read * score_kmer) + (1.0 - score_read) * (1.0 - score_kmer));

		edlibFreeAlignResult(align);

		return (std::isnan(score) ? 0.0 : score);
	};

	const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

	stack<info_traversal> stck;

	stck.push(info_traversal(p, r, false));

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* pli = um_start.getData();

		const bool short_cycle = (/*i_t.c ||*/ um_start.getData()->isInCycle() || um_start.getData()->get_readID().isEmpty());

		Path<UnitigData> best_p;
		PairID best_r;

		double l_best_score = -1;

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			//if (pli->get_count_successor(um_start, succ) >= opt.min_cov_edges) {

				const UnitigData* pli_succ = succ.getData();
				const PairID& r_succ = pli_succ->get_readID();

				if (!um_e.isEmpty && (succ.isSameReferenceUnitig(um_e)) && (um_e.strand == succ.strand) && (succ.dist <= um_e.dist)){

					Path<UnitigData> path(i_t.p);

					path.extend(succ);

					if (out_qual){

						const PairID r_curr = i_t.r | r_succ;
						const double score = getScorePath(path, r_curr);

						path.setQuality(getQual(score, opt.out_qual));
					}

					semi_weak_successors.push_back(move(path));
				}

				PairID r_curr = i_t.r;

				size_t l_count = 0;

				if (r_succ.isEmpty() && !r_curr.isEmpty()) l_count = min_cov_vertex;
				else {

					if (um.getData()->get_readID().isEmpty() && r_curr.isEmpty()) r_curr = r_succ;
					else r_curr &= r_succ;

					l_count = min_cov_vertex & (static_cast<size_t>(!hasEnoughSharedPairID(bf_read, p_read, r_curr, min_cov_vertex)) - 1);
				}

				if (l_count == min_cov_vertex){

					Path<UnitigData> path(i_t.p);

					path.extend(succ);
					
					if (!short_cycle && (path.length() < opt.large_k)) stck.push(info_traversal(path, r_curr, false));
					else {

						const double score = getScorePath(path, r_curr);

						if (path.length() >= opt.large_k){

							if (score > best_score){

								best_score = score;
								best_successors.clear();
							}
							
							if (score >= best_score){

								if (out_qual) path.setQuality(getQual(score, opt.out_qual));

								best_successors.push_back(path);
							}
						}
						
						if (short_cycle) {

							if (score > l_best_score) l_best_score = score;

							if (score >= l_best_score) {

								best_p = path;
								best_r = r_curr;
							}
						}
					}
				}
			//}
		}

		if (short_cycle && (l_best_score >= 0) && (best_p.length() < opt.large_k)) stck.push(info_traversal(best_p, best_r, true));
	}
}

void exploreSubGraphLong(const Correct_Opt& opt, const TinyBloomFilter<uint32_t>& bf_read, const PairID& p_read, const char* ref, const size_t ref_len,
						const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e,
						const Path<UnitigData>& p, const PairID& r, const size_t min_cov_vertex,
						double& best_score, vector<Path<UnitigData>>& best_successors, vector<Path<UnitigData>>& semi_weak_successors) {

	struct info_traversal {

		Path<UnitigData> p;
		PairID r;
		bool c;

		info_traversal(const Path<UnitigData>& p_, const PairID& r_, const bool c_) : p(p_), r(r_), c(c_) {};
	};

	auto getScorePath = [&p_read, &ref, ref_len](const Path<UnitigData>& path, const PairID& p_id_path){

		const string path_curr = path.toString();

		EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		EdlibAlignResult align = edlibAlign(path_curr.c_str(), path_curr.length(), ref, ref_len, config);

		const double score_kmer = 1.0 - (static_cast<double>(align.editDistance) / path_curr.length());
		const double score_read = static_cast<double>(getNumberSharedPairID(p_read, p_id_path)) / p_read.cardinality();
		const double score = (score_read * score_kmer) / ((score_read * score_kmer) + (1.0 - score_read) * (1.0 - score_kmer));

		edlibFreeAlignResult(align);

		return (std::isnan(score) ? 0.0 : score);
	};

	const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

	stack<info_traversal> stck;

	stck.push(info_traversal(p, r, false));

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* pli = um_start.getData();

		const bool short_cycle = (/*i_t.c ||*/ um_start.getData()->isInCycle() || um_start.getData()->get_readID().isEmpty());

		Path<UnitigData> best_p;
		PairID best_r;

		double l_best_score = -1;

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			//if (pli->get_count_successor(um_start, succ) >= opt.min_cov_edges) {

				const UnitigData* pli_succ = succ.getData();
				const PairID& r_succ = pli_succ->get_readID();
				const Kmer head = succ.strand ? succ.getUnitigHead() : succ.getUnitigHead().twin();
				const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it_h_um_e = h_um_e.find(head);

				if (it_h_um_e != h_um_e.end()){

					const const_UnitigMap<UnitigData>& um_e = *it_h_um_e;

					if (succ.dist <= um_e.dist){

						Path<UnitigData> path(i_t.p);

						path.extend(succ);

						if (out_qual){

							const PairID r_curr = i_t.r | r_succ;
							const double score = getScorePath(path, r_curr);

							path.setQuality(getQual(score, opt.out_qual));
						}

						semi_weak_successors.push_back(move(path));
					}
				}

				PairID r_curr = i_t.r;

				size_t l_count = 0;

				if (r_succ.isEmpty() && !r_curr.isEmpty()) l_count = min_cov_vertex;
				else {

					if (um.getData()->get_readID().isEmpty() && r_curr.isEmpty()) r_curr = r_succ;
					else r_curr &= r_succ;

					l_count = min_cov_vertex & (static_cast<size_t>(!hasEnoughSharedPairID(bf_read, p_read, r_curr, min_cov_vertex)) - 1);
				}

				if (l_count == min_cov_vertex){

					Path<UnitigData> path(i_t.p);

					path.extend(succ);
					
					if (!short_cycle && (path.length() < opt.large_k)) stck.push(info_traversal(path, r_curr, false));
					else {

						const double score = getScorePath(path, r_curr);

						if (path.length() >= opt.large_k){

							if (score > best_score){

								best_score = score;
								best_successors.clear();
							}
							
							if (score >= best_score){

								if (out_qual) path.setQuality(getQual(score, opt.out_qual));

								best_successors.push_back(path);
							}
						}
						
						if (short_cycle) {

							if (score > l_best_score) l_best_score = score;

							if (score >= l_best_score) {

								best_p = path;
								best_r = r_curr;
							}
						}
					}
				}
			//}
		}

		if (short_cycle && (l_best_score >= 0) && (best_p.length() < opt.large_k)) stck.push(info_traversal(best_p, best_r, true));
	}
}