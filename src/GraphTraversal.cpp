#include "GraphTraversal.hpp"

pair<vector<Path<UnitigData>>, bool> explorePathsBFS(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const WeightsPairID& w_pid, const const_UnitigMap<UnitigData>& um_s,
														const bool long_read_correct, const uint64_t hap_id){

	const uint64_t undetermined_hap_id = 0xffffffffffffffffULL;

	vector<Path<UnitigData>> v, v_tmp;

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

	auto explore = [&](const const_UnitigMap<UnitigData>& um, const Path<UnitigData>& path, const size_t level) -> vector<Path<UnitigData>> {

		vector<Path<UnitigData>> terminal_paths;
		vector<Path<UnitigData>> non_terminal_paths;

		size_t end_pos_ref = 0;

		if ((path.length() > (um.len + opt.k - 1)) && !um.isEmpty){

			EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0);
			EdlibAlignResult align = edlibAlign(path.toString().c_str(), path.length() - um.len - opt.k + 1, ref, ref_len, config);

			end_pos_ref = align.endLocations[0] + 1;

			edlibFreeAlignResult(align);
		}

		if ((ref_len - end_pos_ref) != 0){

			pair<double, double> scores;

			if (long_read_correct) scores = exploreSubGraphLong(opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, um, const_UnitigMap<UnitigData>(), terminal_paths, non_terminal_paths, hap_id);
			else scores = exploreSubGraph(opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, um, const_UnitigMap<UnitigData>(), level - 1, terminal_paths, non_terminal_paths, hap_id);

			if (non_terminal_paths.size() > 1) {

				if (scores.second != 0.0) non_terminal_paths = selectMostContiguous(non_terminal_paths, w_pid);
				else {

					non_terminal_paths = selectMostContiguous(non_terminal_paths, opt.min_cov_vertices, hap_id);

					if (!non_terminal_paths.empty()) non_terminal_paths = vector<Path<UnitigData>>(1, non_terminal_paths[selectBestSubstringAlignment(ref + end_pos_ref, ref_len - end_pos_ref, non_terminal_paths).first]);
				}
			}
		}

		return non_terminal_paths;
	};

	if (!um_s.isEmpty){

		const size_t level = 4;

		const double seq_entropy = getEntropy(ref, ref_len);

		const size_t min_len_path = getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).first + opt.k;
		const size_t max_len_path = getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).second + opt.k;

		const size_t max_len_subpath = opt.k * opt.large_k_factor;

		//const size_t max_paths = getMaxPaths(seq_entropy, max_len_path, opt.k);
		const size_t max_paths = 1024;

		size_t max_branches = 2 * getMaxBranch(seq_entropy, max_len_path, opt.k);
		size_t max_sz_stck = 512;

		const size_t dec_sz_stck = max_sz_stck / (((max_len_path - opt.k + 1) / max_branches) + 1);

		KmerHashTable<vector<Path<UnitigData>>> km_h;

		queue<Path<UnitigData>> q, q_long;
	    
	    const_UnitigMap<UnitigData> um_start_tmp(um_s); // Create a non-const local copy of unitig given in parameter
	    Path<UnitigData> p_tmp;

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

		    p_tmp.extend(back, 1.0);
    		v.push_back(p_tmp);
    		p_tmp.clear();
		}

	    p_tmp.extend(um_start_tmp, 1.0);
	    q.push(move(p_tmp));

	    while (v.empty()){

		    while (!q.empty()){

		    	if (q.size() >= max_sz_stck) resizeQueue(q);

		        Path<UnitigData> p = move(q.front()); // Get the unitig on top of the stack
		        const const_UnitigMap<UnitigData>& um = p.back();

		        q.pop(); // Delete unitig on the top of the stack

		        if (p.length() < max_len_path){

		        	if (p.size() < max_branches){

		        		const vector<Path<UnitigData>> exploreFW = explore(um, p, level);

		        		bool extended = false;

			        	for (const auto& path : exploreFW){

			        		Path<UnitigData> p_ext(p);

			        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();
			        		const vector<char> v_qual = path.toQualityVector();

			        		for (size_t i = 0; i < v_path.size(); ++i){

			        			p_ext.extend(v_path[i], getScore(v_qual[i]));

			        			if ((p_ext.length() >= min_len_path) && (p_ext.length() <= max_len_path)){

			        				v_tmp.push_back(p_ext);
			        				extended = true;
			        			}
			        		}

			        		if ((!long_read_correct && (path.size() == level)) || (long_read_correct && (path.length() >= max_len_subpath))){

				        		q.push(move(p_ext));
				        		extended = true;
				        	}
		        		}

		        		if (!extended && (p.length() <= max_len_path)) v_tmp.push_back(p);

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
			}
			else break;
		}
	}

	return {v, true};
}

pair<vector<Path<UnitigData>>, bool> explorePathsBFS2(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const WeightsPairID& w_pid, const const_UnitigMap<UnitigData>& um_s, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_e,
														const bool long_read_correct, const uint64_t hap_id){

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

	auto explore = [&](const const_UnitigMap<UnitigData>& um, const Path<UnitigData>& path, const size_t level) -> pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> {

		vector<Path<UnitigData>> terminal_paths;
		vector<Path<UnitigData>> non_terminal_paths;

		size_t end_pos_ref = 0;

		if ((path.length() > (um.len + opt.k - 1)) && !um.isEmpty){

			EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0);
			EdlibAlignResult align = edlibAlign(path.toString().c_str(), path.length() - um.len - opt.k + 1, ref, ref_len, config);

			end_pos_ref = align.endLocations[0] + 1;

			edlibFreeAlignResult(align);
		}

		if ((ref_len - end_pos_ref) != 0) {

			pair<double, double> scores;

			if (long_read_correct) scores = exploreSubGraphLong(opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, um, um_h, terminal_paths, non_terminal_paths, hap_id);
			else scores = exploreSubGraph(opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, um, um_h, level - 1, terminal_paths, non_terminal_paths, hap_id);

			if (non_terminal_paths.size() > 1) {

				if (scores.second != 0.0) non_terminal_paths = selectMostContiguous(non_terminal_paths, w_pid);
				else {

					non_terminal_paths = selectMostContiguous(non_terminal_paths, opt.min_cov_vertices, hap_id);
					
					if (!non_terminal_paths.empty()) non_terminal_paths = vector<Path<UnitigData>>(1, non_terminal_paths[selectBestSubstringAlignment(ref + end_pos_ref, ref_len - end_pos_ref, non_terminal_paths).first]);
				}
			}
		}

    	return {terminal_paths, non_terminal_paths};
	};

	vector<Path<UnitigData>> v, v_tmp;

	if (!um_s.isEmpty){

		const size_t level = 4;

		const double seq_entropy = getEntropy(ref, ref_len);

		const size_t min_len_path = getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).first + opt.k;
		const size_t max_len_path = getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).second + opt.k;

		const size_t max_len_subpath = opt.k * opt.large_k_factor;

		const size_t max_paths = 1024;

		size_t max_branches = 2 * getMaxBranch(seq_entropy, max_len_path, opt.k);
		size_t max_sz_stck = 512;

		const size_t dec_sz_stck = max_sz_stck / (((max_len_path - opt.k + 1) / max_branches) + 1);

		KmerHashTable<pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>>> km_h;

		queue<Path<UnitigData>> q, q_long;
	    
	    const_UnitigMap<UnitigData> um_start_tmp(um_s); // Create a non-const local copy of unitig given in parameter
	    Path<UnitigData> p_start_tmp;

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

				    p_start_tmp.extend(back_tmp, 1.0);

	        		v.push_back(p_start_tmp);
	        		p_start_tmp.clear();
				}
			}
		}

		p_start_tmp.extend(um_start_tmp, 1.0);
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

						const pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> exploreFW = explore(um, p, level);

			        	for (const auto& path : exploreFW.first){

			        		Path<UnitigData> p_ext(p);

			        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();
			        		const vector<char> v_qual = path.toQualityVector();

			        		for (size_t i = 0; i < v_path.size(); ++i) p_ext.extend(v_path[i], getScore(v_qual[i]));

			        		v_tmp.push_back(move(p_ext));
		        		}

			        	for (const auto& path : exploreFW.second){

			        		if ((!long_read_correct && (path.size() == level)) || (long_read_correct && (path.length() >= max_len_subpath))){

				        		Path<UnitigData> p_ext(p);

				        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();
				        		const vector<char> v_qual = path.toQualityVector();

				        		for (size_t i = 0; i < v_path.size(); ++i) p_ext.extend(v_path[i], getScore(v_qual[i]));

			        			q.push(move(p_ext));
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
			}
			else break;
		}
	}

	return {v, true};
}

pair<double, double> exploreSubGraph(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
										const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e, const size_t level,
										vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id) {

	const uint64_t undetermined_hap_id = 0xffffffffffffffffULL;

	double score_t = 0.0, score_nt = 0.0;
	double score_align_t = 0.0, score_align_nt = 0.0;

	stack<info_traversal> stck;

	stck.push(info_traversal(level));

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* pli = um_start.getData();
		const PairID& r_start = pli->get_readID();

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			const UnitigData* ud_succ = succ.getData();
			const PairID& r_succ = ud_succ->get_readID();
			const PairID& hap_succ = ud_succ->get_hapID();

			if (hasEnoughSharedPairID(r_start, r_succ, opt.min_cov_vertices)){

				Path<UnitigData> path(i_t.p);

				path.extend(succ);

				if (!um_e.isEmpty && (succ.isSameReferenceUnitig(um_e)) && (um_e.strand == succ.strand) && (succ.dist <= um_e.dist)){

					const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

					if (scores.first >= score_t) {

						path.setQuality(getQual(scores.first, opt.out_qual));

						if (scores.first > score_t) terminal_paths.clear();
						terminal_paths.push_back(path);

						score_t = scores.first;
						score_align_t = scores.second;
					}
				}
				
				if (i_t.l != 0) stck.push(info_traversal(path, i_t.l - 1));
				else {

					const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

					if (scores.first >= score_nt) {

						path.setQuality(getQual(scores.first, opt.out_qual));

						if (scores.first > score_nt) non_terminal_paths.clear();

						non_terminal_paths.push_back(path);

						score_nt = scores.first;
						score_align_nt = scores.second;
					}
				}
			}
		}
	}

	return {score_t, score_nt};
}

pair<double, double> exploreSubGraph(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
										const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e, const size_t level,
										vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id) {

	const uint64_t undetermined_hap_id = 0xffffffffffffffffULL;

	double score_t = 0.0, score_nt = 0.0;
	double score_align_t = 0.0, score_align_nt = 0.0;

	stack<info_traversal> stck;

	stck.push(info_traversal(level));

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* pli = um_start.getData();
		const PairID& r_start = pli->get_readID();

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			const UnitigData* ud_succ = succ.getData();
			const PairID& r_succ = ud_succ->get_readID();
			const PairID& hap_succ = ud_succ->get_hapID();

			if (hasEnoughSharedPairID(r_start, r_succ, opt.min_cov_vertices)){

				const Kmer head = succ.strand ? succ.getUnitigHead() : succ.getUnitigHead().twin();
				const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it_h_um_e = h_um_e.find(head);

				Path<UnitigData> path(i_t.p);

				path.extend(succ);

				if (it_h_um_e != h_um_e.end()){

					const const_UnitigMap<UnitigData>& um_e = *it_h_um_e;

					if (succ.dist <= um_e.dist){

						const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

						if (scores.first >= score_t) {

							path.setQuality(getQual(scores.first, opt.out_qual));

							if (scores.first > score_t) terminal_paths.clear();
							terminal_paths.push_back(path);

							score_t = scores.first;
							score_align_t = scores.second;
						}
					}
				}
				
				if (i_t.l != 0) stck.push(info_traversal(path, i_t.l - 1));
				else {

					const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

					if (scores.first >= score_nt) {

						path.setQuality(getQual(scores.first, opt.out_qual));

						if (scores.first > score_nt) non_terminal_paths.clear();

						non_terminal_paths.push_back(path);

						score_nt = scores.first;
						score_align_nt = scores.second;
					}
				}
			}
		}
	}

	return {score_t, score_nt};
}

pair<double, double> exploreSubGraphLong(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
											const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e,
											vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id) {

	const size_t max_len_subpath = opt.k * opt.large_k_factor;

	double score_t = 0.0, score_nt = 0.0;
	double score_align_t = 0.0, score_align_nt = 0.0;

	stack<info_traversal> stck;

	stck.push(info_traversal());

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* ud = um_start.getData();
		const PairID& r_start = ud->get_readID();
		const size_t km_cov_start = ud->getKmerCoverage(um_start);

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			const UnitigData* ud_succ = succ.getData();
			const PairID& r_succ = ud_succ->get_readID();
			const size_t km_cov_succ = ud_succ->getKmerCoverage(succ);

			if ((km_cov_start >= opt.max_cov_vertices) || (km_cov_succ >= opt.max_cov_vertices) || hasEnoughSharedPairID(r_start, r_succ, opt.min_cov_vertices)){

				Path<UnitigData> path(i_t.p);

				path.extend(succ);

				if (!um_e.isEmpty && (succ.isSameReferenceUnitig(um_e)) && (um_e.strand == succ.strand) && (succ.dist <= um_e.dist)){

					const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

					if (scores.first >= score_t) {

						path.setQuality(getQual(scores.first, opt.out_qual));

						if (scores.first > score_t) terminal_paths.clear();
						terminal_paths.push_back(path);

						score_t = scores.first;
						score_align_t = scores.second;
					}
				}
				
				if (path.length() < max_len_subpath) stck.push(info_traversal(path));
				else {

					const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

					if (scores.first >= score_nt) {

						path.setQuality(getQual(scores.first, opt.out_qual));

						if (scores.first > score_nt) non_terminal_paths.clear();

						non_terminal_paths.push_back(path);

						score_nt = scores.first;
						score_align_nt = scores.second;
					}
				}
			}
		}
	}

	return {score_t, score_nt};
}

pair<double, double> exploreSubGraphLong(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len,
											const const_UnitigMap<UnitigData>& um, const KmerHashTable<const_UnitigMap<UnitigData>>& h_um_e,
											vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id) {

	const size_t max_len_subpath = opt.k * opt.large_k_factor;

	double score_t = 0.0, score_nt = 0.0;
	double score_align_t = 0.0, score_align_nt = 0.0;

	stack<info_traversal> stck;

	stck.push(info_traversal());

	while (!stck.empty()) {

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* ud = um_start.getData();
		const PairID& r_start = ud->get_readID();
		const size_t km_cov_start = ud->getKmerCoverage(um_start);

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			const UnitigData* ud_succ = succ.getData();
			const PairID& r_succ = ud_succ->get_readID();
			const size_t km_cov_succ = ud_succ->getKmerCoverage(succ);

			if ((km_cov_start >= opt.max_cov_vertices) || (km_cov_succ >= opt.max_cov_vertices) || hasEnoughSharedPairID(r_start, r_succ, opt.min_cov_vertices)){

				const Kmer head = succ.strand ? succ.getUnitigHead() : succ.getUnitigHead().twin();
				const KmerHashTable<const_UnitigMap<UnitigData>>::const_iterator it_h_um_e = h_um_e.find(head);

				Path<UnitigData> path(i_t.p);

				path.extend(succ);

				if (it_h_um_e != h_um_e.end()){

					const const_UnitigMap<UnitigData>& um_e = *it_h_um_e;

					if (succ.dist <= um_e.dist){

						const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

						if (scores.first >= score_t) {

							path.setQuality(getQual(scores.first, opt.out_qual));

							if (scores.first > score_t) terminal_paths.clear();
							terminal_paths.push_back(path);

							score_t = scores.first;
							score_align_t = scores.second;
						}
					}
				}
				
				if (path.length() < max_len_subpath) stck.push(info_traversal(path));
				else {

					const pair<double, double> scores = getScorePath(path, w_pid, hap_id, ref, ref_len);

					if (scores.first >= score_nt) {

						path.setQuality(getQual(scores.first, opt.out_qual));

						if (scores.first > score_nt) non_terminal_paths.clear();

						non_terminal_paths.push_back(path);

						score_nt = scores.first;
						score_align_nt = scores.second;
					}
				}
			}
		}
	}

	return {score_t, score_nt};
}

pair<double, double> getScorePath(const Path<UnitigData>& path, const WeightsPairID& w_pid, const uint64_t hap_id, const char* ref, const size_t ref_len){

	double score = 0.0;
	double score_align = 0.0;
	double score_pids = 0.0;

	if (path.length() != 0){

		const Path<UnitigData>::PathOut path_str_um = path.toStringVector();
		const string& path_str = path_str_um.toString();
		const vector<const_UnitigMap<UnitigData>>& v_um = path_str_um.toVector();
		
		EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
		EdlibAlignResult align;

		if (path_str.length() >= ref_len){

			align = edlibAlign(ref, ref_len, path_str.c_str(), path_str.length(), config);
			score_align = 1.0 - (static_cast<double>(align.editDistance) / ref_len);
		}
		else {

			align = edlibAlign(path_str.c_str(), path_str.length(), ref, ref_len, config);
			score_align = 1.0 - (static_cast<double>(align.editDistance) / path_str.length());
		}

		edlibFreeAlignResult(align);

		if (w_pid.weighted_pids.isEmpty() && w_pid.noWeight_pids.isEmpty()) score = 0.0;
		else {

			for (const auto& um : v_um){

				const PairID& pids = um.getData()->get_readID();

				if (!pids.isEmpty()) {

					const double shared_weighted = static_cast<double>(getNumberSharedPairID(pids, w_pid.weighted_pids));
					const double shared_noWeight = static_cast<double>(getNumberSharedPairID(pids, w_pid.noWeight_pids));

					score_pids += (shared_weighted * w_pid.weight) + shared_noWeight;
				}
			}

			score_pids /= (w_pid.sum_pids_weights * v_um.size());

			score = score_pids * score_align; // Conflate the two probabilities (1/2)
			score /= score + ((1.0 - score_pids) * (1.0 - score_align)); // Conflate the two probabilities (2/2)
			score = (std::isnan(score) ? 0.0 : score);
		}
	}

	return pair<double, double>(score, score_align);
}

vector<Path<UnitigData>> selectMostContiguous(const vector<Path<UnitigData>>& v_paths, const WeightsPairID& w_pid) {

	double score_contiguity = 0;

	vector<Path<UnitigData>> v_paths_out;

	for (const auto& path : v_paths) {

		const vector<const_UnitigMap<UnitigData>> v_um = path.toVector();

		PairID pid_weighted = (v_um.front().getData()->get_readID() & w_pid.weighted_pids);
		PairID pid_no_weights = (v_um.front().getData()->get_readID() & w_pid.noWeight_pids);

		for (const auto& um : v_um){

			if (!pid_weighted.isEmpty()) pid_weighted &= (um.getData()->get_readID() & w_pid.weighted_pids);
			if (!pid_no_weights.isEmpty()) pid_no_weights &= (um.getData()->get_readID() & w_pid.noWeight_pids);
		}

		const double l_score_contiguity = static_cast<double>(pid_weighted.cardinality()) * w_pid.weight + static_cast<double>(pid_no_weights.cardinality());

		if (l_score_contiguity >= score_contiguity) {

			if (l_score_contiguity > score_contiguity) v_paths_out.clear();

			v_paths_out.push_back(path);

			score_contiguity = l_score_contiguity;
		}
	}

	return v_paths_out;
}

vector<Path<UnitigData>> selectMostContiguous(const vector<Path<UnitigData>>& v_paths, const size_t min_cov, const uint64_t hap_id) {

	vector<Path<UnitigData>> v_paths_out;

	for (const auto& path : v_paths) {

		const vector<const_UnitigMap<UnitigData>> v_um = path.toVector();

		PairID pid = v_um.front().getData()->get_readID();

		for (const auto& um : v_um){

			if (pid.cardinality() >= min_cov) pid &= um.getData()->get_readID();
		}

		if (pid.cardinality() >= min_cov) v_paths_out.push_back(path);
	}

	return v_paths_out;
}

vector<pair<size_t, char>> getAmbiguityVector(const vector<const_UnitigMap<UnitigData>>& v_um, const size_t k) {

	vector<pair<size_t, char>> v_amb;

	size_t prev_l = 0;
	size_t pos_prev_l = 0;

	for (const auto& um : v_um){

		const vector<pair<size_t, char>> v_amb_um = um.getData()->get_ambiguity_char(um);

		vector<pair<size_t, char>> v_amb_tmp;

		auto it_prev = v_amb.begin() + pos_prev_l;
		auto it_curr = v_amb_um.begin();

		while ((it_prev != v_amb.end()) && (it_curr != v_amb_um.end()) && (it_curr->first < k-1)){

			const size_t it_curr_pos = it_curr->first + prev_l;

			if (it_prev->first < it_curr_pos){

				v_amb_tmp.push_back(*it_prev);
				++it_prev;
			}
			else if (it_prev->first > it_curr_pos){

				v_amb_tmp.push_back({it_curr_pos, it_curr->second});
				++it_curr;
			}
			else {

				bool nuc_a_s, nuc_c_s, nuc_g_s, nuc_t_s;
				bool nuc_a_v, nuc_c_v, nuc_g_v, nuc_t_v;

				getAmbiguityRev(it_prev->second, nuc_a_s, nuc_c_s, nuc_g_s, nuc_t_s);
				getAmbiguityRev(it_curr->second, nuc_a_v, nuc_c_v, nuc_g_v, nuc_t_v);

				v_amb_tmp.push_back({it_prev->first, getAmbiguity((nuc_a_s || nuc_a_v), (nuc_c_s || nuc_c_v), (nuc_g_s || nuc_g_v), (nuc_t_s || nuc_t_v))});

				++it_prev;
				++it_curr;
			}
		}

		while (it_prev != v_amb.end()){

			v_amb_tmp.push_back(*it_prev);
			++it_prev;
		}

		while (it_curr != v_amb_um.end()){

			v_amb_tmp.push_back({it_curr->first + prev_l, it_curr->second});
			++it_curr;
		}

		prev_l += um.len;

		v_amb.erase(v_amb.begin() + pos_prev_l, v_amb.end());

		for (auto& p : v_amb_tmp){

			v_amb.push_back(move(p));

			pos_prev_l += static_cast<size_t>(v_amb.back().first < prev_l);
		}
	}

	return v_amb;
}

vector<pair<size_t, char>> getAmbiguityVector(const Path<UnitigData>::PathOut& path, const size_t k, const bool rm_solid_amb) {

	vector<pair<size_t, char>> v1 = getAmbiguityVector(path.toVector(), k), v2;

	if (rm_solid_amb) {

		const size_t len = path.toString().length() - k;

		for (const auto p : v1) {

			if ((p.first >= k) && (p.first < len)) v2.push_back(p);
		}

		v1 = move(v2);
	}

	return v1;
}

void annotateConnectedComponents(CompactedDBG<UnitigData>& dbg) {

	size_t cc_id = 1;

	for (auto& um : dbg){

		if (um.getData()->getConnectedComp() == 0) {

			unordered_set<Kmer, KmerHash> km_s;
			vector<UnitigMap<UnitigData>> v_um;

			v_um.push_back(um);
			km_s.insert(um.getUnitigHead());

			while (!v_um.empty()) {

				UnitigMap<UnitigData> um_pop = v_um.back();

				v_um.pop_back();
				um_pop.getData()->setConnectedComp(cc_id);

				for (const auto& um_pred: um_pop.getPredecessors()){

					if (km_s.insert(um_pred.getUnitigHead()).second) v_um.push_back(um_pred);
				}

				for (const auto& um_succ: um_pop.getSuccessors()){

					if (km_s.insert(um_succ.getUnitigHead()).second) v_um.push_back(um_succ);
				}
			}
		}
	}
}

bool isValidSNPcandidate(local_graph_traversal& lgt, const const_UnitigMap<UnitigData>& um_a, const const_UnitigMap<UnitigData>& um_b, const size_t limit_length_path, const size_t limit_sz_stack){

	const PairID& p_ids_a = um_a.getData()->get_readID();
	const PairID& p_ids_b = um_b.getData()->get_readID();

	const size_t k = um_a.getGraph()->getK();
	const size_t max_len_path = limit_length_path + (um_a.size - k + 1);

	auto exploreLocalGraphFW = [&]() {

		if (lgt.m_km_fw.empty()) {

			lgt.q_um_fw.push({um_a, {0, 0}});
			lgt.m_km_fw.insert({(um_a.strand ? um_a.getUnitigHead() : um_a.getUnitigTail().twin()), 0});

			lgt.p_ids_fw = p_ids_a;
		}

		while (!lgt.q_um_fw.empty()) {

			const pair<const_UnitigMap<UnitigData>, pair<size_t, size_t>> elem_pop = lgt.q_um_fw.front();

			const const_UnitigMap<UnitigData>& um_pop = elem_pop.first;

			const Kmer km_pop = um_pop.strand ? um_pop.getUnitigHead() : um_pop.getUnitigTail().twin();

			const size_t nb_nodes_before = elem_pop.second.second;
			const size_t nb_nodes_after = nb_nodes_before + 1;

			const size_t dist_before = elem_pop.second.first;
			const size_t dist_after = dist_before + (um_pop.size - k + 1);

			lgt.q_um_fw.pop();

			if (dist_before <= lgt.m_km_fw.find(km_pop)->second) {

				for (const auto& um : um_pop.getSuccessors()){

					pair<unordered_map<Kmer, size_t, KmerHash>::iterator, bool> p_km_m = lgt.m_km_fw.insert({um.getMappedHead(), dist_after});

					if (p_km_m.second){ // Unitig hasn't been visited so far

						const PairID& p_ids_succ = um.getData()->get_readID();

						if (hasEnoughSharedPairID(p_ids_succ, p_ids_a, 1)){ // Unitig shares enough reads with start unitig

							lgt.p_ids_fw |= p_ids_succ;

							if (dist_after < max_len_path) lgt.q_um_fw.push({um, {dist_after, nb_nodes_after}});
							if (hasEnoughSharedPairID(p_ids_succ, p_ids_b, 1)) return true;
						}
						else p_km_m.first->second = 0xffffffffffffffffULL; // Unitig is visited but is invalid
					}
					else if ((p_km_m.first->second != 0xffffffffffffffffULL) && (dist_after < p_km_m.first->second)) { // Unitig is valid, has been visited already but distance is shorter

						p_km_m.first->second = dist_after;

						if (dist_after < max_len_path) lgt.q_um_fw.push({um, {dist_after, nb_nodes_after}});
					}
				}

				if (lgt.m_km_fw.size() >= limit_sz_stack) return true;
			}
		}

		return false;
	};

	auto exploreLocalGraphBW = [&]() {

		if (lgt.m_km_bw.empty()) {

			lgt.q_um_bw.push({um_a, {0, 0}});
			lgt.m_km_bw.insert({(um_a.strand ? um_a.getUnitigHead() : um_a.getUnitigTail().twin()), 0});

			lgt.p_ids_bw = p_ids_a;
		}

		while (!lgt.q_um_bw.empty()) {

			const pair<const_UnitigMap<UnitigData>, pair<size_t, size_t>> elem_pop = lgt.q_um_bw.front();

			const const_UnitigMap<UnitigData>& um_pop = elem_pop.first;

			const Kmer km_pop = um_pop.strand ? um_pop.getUnitigHead() : um_pop.getUnitigTail().twin();

			const size_t nb_nodes_before = elem_pop.second.second;
			const size_t nb_nodes_after = nb_nodes_before + 1;

			const size_t dist_before = elem_pop.second.first;
			const size_t dist_after = dist_before + (um_pop.size - k + 1);

			lgt.q_um_bw.pop();

			if (dist_before <= lgt.m_km_bw.find(km_pop)->second) {

				for (const auto& um : um_pop.getPredecessors()){

					pair<unordered_map<Kmer, size_t, KmerHash>::iterator, bool> p_km_m = lgt.m_km_bw.insert({um.getMappedHead(), dist_after});

					if (p_km_m.second){ // Unitig hasn't been visited so far

						const PairID& p_ids_succ = um.getData()->get_readID();

						if (hasEnoughSharedPairID(p_ids_succ, p_ids_a, 1)){ // Unitig shares enough reads with start unitig

							lgt.p_ids_bw |= p_ids_succ;

							if (dist_after < limit_length_path) lgt.q_um_bw.push({um, {dist_after, nb_nodes_after}});
							if (hasEnoughSharedPairID(p_ids_succ, p_ids_b, 1)) return true;
						}
						else p_km_m.first->second = 0xffffffffffffffffULL; // Unitig is visited but is invalid
					}
					else if ((p_km_m.first->second != 0xffffffffffffffffULL) && (dist_after < p_km_m.first->second)) { // Unitig is valid, has been visited already but distance is shorter

						p_km_m.first->second = dist_after;

						if (dist_after < limit_length_path) lgt.q_um_bw.push({um, {dist_after, nb_nodes_after}});
					}
				}

				if (lgt.m_km_bw.size() >= limit_sz_stack) return true;
			}
		}

		return false;
	};

	if (hasEnoughSharedPairID(lgt.p_ids_fw, p_ids_b, 1) || hasEnoughSharedPairID(lgt.p_ids_bw, p_ids_b, 1)) return true;
	else if ((lgt.m_km_fw.size() >= limit_sz_stack) || (lgt.m_km_bw.size() >= limit_sz_stack)) return true;
	
	return (exploreLocalGraphFW() || exploreLocalGraphBW());
}