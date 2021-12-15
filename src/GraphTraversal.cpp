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

	auto explore = [&](	const const_UnitigMap<UnitigData>& um, const Path<UnitigData>& path, const size_t max_len_path, const size_t level,
						unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid) -> vector<Path<UnitigData>> {

		const size_t path_len = path.length();

		const bool non_empty_path = ((path_len > (um.len + opt.k - 1)) && !um.isEmpty);

		const size_t path_len_prefix = (non_empty_path ? (path_len - um.len - opt.k + 1) : 0);

		vector<Path<UnitigData>> terminal_paths;
		vector<Path<UnitigData>> non_terminal_paths;

		size_t end_pos_ref = 0;

		if (non_empty_path){

			EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
			EdlibAlignResult align = edlibAlign(path.toString().c_str(), path_len_prefix, ref, ref_len, config);

			end_pos_ref = align.endLocations[0] + 1;

			edlibFreeAlignResult(align);
		}

		// Suffix of the ref. seq. uncovered by corrected seq. and path to extend has not reached max. possible length yet.
		if (((ref_len - end_pos_ref) != 0) && (path_len < max_len_path)) {

			pair<double, double> scores;

			const size_t l_max_len_path = max_len_path - path_len_prefix;

			if (long_read_correct) {

				scores = exploreSubGraphLong(	opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, l_max_len_path,
												um, const_UnitigMap<UnitigData>(), terminal_paths, non_terminal_paths, hap_id,
												m_pid);
			}
			else {

				scores = exploreSubGraph(	opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, l_max_len_path,
											um, const_UnitigMap<UnitigData>(), level - 1, terminal_paths, non_terminal_paths, hap_id,
											m_pid);
			}

			// Multi-round correction
			if (!terminal_paths.empty() && (scores.first < opt.min_score)) terminal_paths.clear();
			if (!non_terminal_paths.empty() && (scores.second < opt.min_score)) non_terminal_paths.clear();

			if (non_terminal_paths.size() > 1) non_terminal_paths = vector<Path<UnitigData>>(1, non_terminal_paths[selectBestSubstringAlignment(ref + end_pos_ref, ref_len - end_pos_ref, non_terminal_paths).first]);
		}

		return non_terminal_paths;
	};

	if (!um_s.isEmpty && um_s.getData()->hasSharedPids()){

		const size_t level = 4;

		const size_t min_len_path = getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).first + opt.k;
		const size_t max_len_path = max(getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).second, 10UL) + opt.k;

		const size_t max_len_subpath = opt.k * opt.large_k_factor;

		const size_t max_paths = 1024;
		const size_t max_sz_stck = 512;

		unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr> m_pid;

		queue<Path<UnitigData>> q;
	    
	    // Push start unitig suffix in the graph traversal queue
	    {
		    const_UnitigMap<UnitigData> um_start_tmp(um_s);
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

			    p_tmp.extend(back, string(back.len + opt.k - 1, getQual(1.0)));
	    		v.push_back(move(p_tmp));
			}

		    p_tmp.extend(um_start_tmp, string(um_start_tmp.len + opt.k - 1, getQual(1.0)));
		    q.push(move(p_tmp));
		}

	    while (!q.empty()){

	        Path<UnitigData> p = move(q.front()); // Get the unitig on top of the stack

	        const const_UnitigMap<UnitigData>& um = p.back();

	        q.pop(); // Delete unitig on the top of the stack

	        if (p.length() < max_len_path){

        		const vector<Path<UnitigData>> exploreFW = explore(um, p, max_len_path, level, m_pid);

	        	for (const auto& path : exploreFW){

	        		Path<UnitigData> p_ext(p);

	        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();
	        		const string& s_qual = path.toQualityString();

	        		for (size_t i = 0, j = 0; i < v_path.size(); ++i){

	        			p_ext.extend(v_path[i], s_qual.substr(j, v_path[i].len + opt.k - 1));

	        			if ((p_ext.length() >= min_len_path) && (p_ext.length() <= max_len_path)) v_tmp.push_back(p_ext);

	        			j += v_path[i].len;
	        		}

	        		if ((!long_read_correct && (path.size() == level)) || (long_read_correct && (path.length() >= max_len_subpath))){

		        		q.push(move(p_ext));

		        		if (q.size() >= max_sz_stck) resizeQueue(q);
		        	}
        		}

        		if (v_tmp.size() >= max_paths){

        			for (auto& p_tmp : v_tmp){

						p_tmp.prunePrefix(max_len_path);
			    		v.push_back(move(p_tmp));
					}

					v_tmp.clear();
				}
	        }
	    }

		for (auto& p_tmp : v_tmp){

			p_tmp.prunePrefix(max_len_path);
    		v.push_back(move(p_tmp));
		}
	}

	if (!v.empty()) {

		if (v.size() > 1) v = vector<Path<UnitigData>>(1, v[selectBestAlignment(v, ref, ref_len).first]);
		
		v = fixRepeats(opt, v, ref, ref_len); // Fix repeats
	}

	return {v, true};
}

pair<vector<Path<UnitigData>>, bool> explorePathsBFS2(	const Correct_Opt& opt, const char* ref, const size_t ref_len,
														const WeightsPairID& w_pid, const const_UnitigMap<UnitigData>& um_s, const const_UnitigMap<UnitigData>& um_e,
														const bool long_read_correct, const uint64_t hap_id){

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

	auto explore = [&](	const const_UnitigMap<UnitigData>& um, const Path<UnitigData>& path, const size_t max_len_path, const size_t level,
						unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid) -> pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> {

		const size_t path_len = path.length();

		const bool non_empty_path = ((path_len > (um.len + opt.k - 1)) && !um.isEmpty);

		const size_t path_len_prefix = (non_empty_path ? (path_len - um.len - opt.k + 1) : 0);

		vector<Path<UnitigData>> terminal_paths;
		vector<Path<UnitigData>> non_terminal_paths;

		size_t end_pos_ref = 0;

		if (non_empty_path){

			EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
			EdlibAlignResult align = edlibAlign(path.toString().c_str(), path_len_prefix, ref, ref_len, config);

			end_pos_ref = align.endLocations[0] + 1;

			edlibFreeAlignResult(align);
		}

		// Suffix of the ref. seq. uncovered by corrected seq. and path to extend has not reached max. possible length yet.
		if (((ref_len - end_pos_ref) != 0) && (path_len < max_len_path)) {

			pair<double, double> scores;

			const size_t l_max_len_path = max_len_path - path_len_prefix;

			if (long_read_correct) {

				scores = exploreSubGraphLong(	opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, l_max_len_path,
												um, um_e, terminal_paths, non_terminal_paths, hap_id, m_pid);
			}
			else {

				scores = exploreSubGraph(	opt, w_pid, ref + end_pos_ref, ref_len - end_pos_ref, l_max_len_path,
											um, um_e, level - 1, terminal_paths, non_terminal_paths, hap_id, m_pid);
			}

			// Multi-round correction
			if (!terminal_paths.empty() && (scores.first < opt.min_score)) terminal_paths.clear();
			if (!non_terminal_paths.empty() && (scores.second < opt.min_score)) non_terminal_paths.clear();

			if (non_terminal_paths.size() > 1) {

				non_terminal_paths = vector<Path<UnitigData>>(1, non_terminal_paths[selectBestSubstringAlignment(ref + end_pos_ref, ref_len - end_pos_ref, non_terminal_paths).first]);
			}
		}

    	return {terminal_paths, non_terminal_paths};
	};

	vector<Path<UnitigData>> v, v_tmp;

	if (!um_s.isEmpty && !um_e.isEmpty && um_s.getData()->hasSharedPids() && um_e.getData()->hasSharedPids()){

		const size_t level = 4;

		const size_t min_len_path = getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).first + opt.k;
		const size_t max_len_path = max(getMinMaxLength(ref_len - opt.k, opt.weak_region_len_factor).second, 10UL) + opt.k;

		const size_t max_len_subpath = opt.k * opt.large_k_factor;

		const size_t max_paths = 1024;
		const size_t max_sz_stck = 512;

		unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr> m_pid;

		queue<Path<UnitigData>> q;
	    
	    {
		    const_UnitigMap<UnitigData> um_start_tmp(um_s); // Create a non-const local copy of unitig given in parameter

		    Path<UnitigData> p_start_tmp;

			if (um_start_tmp.strand){

				um_start_tmp.dist += um_start_tmp.len - 1;
				um_start_tmp.len = um_start_tmp.size - um_start_tmp.dist - opt.k + 1;
			}
		    else {

		    	um_start_tmp.len = um_s.dist + 1;
		    	um_start_tmp.dist = 0;
		    }

			if (um_s.isSameReferenceUnitig(um_e) && (um_s.strand == um_e.strand) && (um_start_tmp.dist <= um_e.dist)) {

				const size_t len = (um_start_tmp.len + opt.k - 1) - (um_e.strand ? um_e.size - um_e.dist - opt.k : um_e.dist);

				if ((len >= min_len_path) && (len <= max_len_path)) {

		    		const_UnitigMap<UnitigData> back_tmp(um_start_tmp);

		    		if (back_tmp.strand) back_tmp.len = um_e.dist - back_tmp.dist + 1;
				    else {

				    	back_tmp.dist = um_e.dist;
				    	back_tmp.len -= um_e.dist;
				    }

				    p_start_tmp.extend(back_tmp, string(back_tmp.len + opt.k - 1, getQual(1.0)));
	        		v.push_back(move(p_start_tmp));
				}
			}

			p_start_tmp.extend(um_start_tmp, string(um_start_tmp.len + opt.k - 1, getQual(1.0)));
		    q.push(move(p_start_tmp));
		}

		while (!q.empty()){

	        Path<UnitigData> p = move(q.front()); // Get the unitig on top of the stack

	        const const_UnitigMap<UnitigData>& um = p.back();
	        const Kmer head = um.strand ? um.getUnitigHead() : um.getUnitigHead().twin();

	        q.pop(); // Delete unitig on the top of the stack

	        if (p.length() < max_len_path){

				const pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> exploreFW = explore(um, p, max_len_path, level, m_pid);

	        	for (const auto& path : exploreFW.first){

	        		Path<UnitigData> p_ext(p);

	        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();
	        		const string& s_qual = path.toQualityString();

	        		for (size_t i = 0, j = 0; i < v_path.size(); ++i){

	        			p_ext.extend(v_path[i], s_qual.substr(j, v_path[i].len + opt.k - 1));
	        			j += v_path[i].len;
	        		}

	        		v_tmp.push_back(move(p_ext));
        		}

	        	for (const auto& path : exploreFW.second){

	        		if ((!long_read_correct && (path.size() == level)) || (long_read_correct && (path.length() >= max_len_subpath))){

		        		Path<UnitigData> p_ext(p);

		        		const vector<const_UnitigMap<UnitigData>> v_path = path.toVector();
		        		const string& s_qual = path.toQualityString();

		        		for (size_t i = 0, j = 0; i < v_path.size(); ++i){

		        			p_ext.extend(v_path[i], s_qual.substr(j, v_path[i].len + opt.k - 1));
		        			j += v_path[i].len;
		        		}

		        		q.push(move(p_ext));

	        			if (q.size() >= max_sz_stck) resizeQueue(q);
	        		}
        		}
        		
        		if (v_tmp.size() >= max_paths){

				    for (auto& p_tmp : v_tmp) {
				    	
						if ((p_tmp.length() >= min_len_path) && (p_tmp.length() <= max_len_path)) {

							if (v.size() + 1 >= max_paths) resizeVector(v);

							v.push_back(move(p_tmp));
						}
					}

					v_tmp.clear();
        		}
	        }
	    }

	    {
		    for (auto& p_tmp : v_tmp){

				if ((p_tmp.length() >= min_len_path) && (p_tmp.length() <= max_len_path)) {

					if (v.size() + 1 >= max_paths) resizeVector(v);

					v.push_back(move(p_tmp));
				}
			}

			v_tmp.clear();
		}
	}

	if (!v.empty()) {

		if (v.size() > 1) v = vector<Path<UnitigData>>(1, v[selectBestAlignment(v, ref, ref_len).first]);

		v = fixRepeats(opt, v, ref, ref_len); // Fix repeats
	}

	return {v, true};
}

pair<double, double> exploreSubGraph(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len, const size_t max_len_path,
										const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e, const size_t level,
										vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id,
										unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid) {

	const uint64_t undetermined_hap_id = 0xffffffffffffffffULL;

	double score_t1 = 0.0, score_nt1 = 0.0;
	double score_t2 = 0.0, score_nt2 = 0.0;

	stack<info_traversal> stck;

	unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr> l_m_pid = m_pid;

	stck.push(info_traversal(level));

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* ud = um_start.getData();

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			const UnitigData* ud_succ = succ.getData();
			const SharedPairID& spid_succ = ud_succ->getPairID();

			const pair<unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>::iterator, bool> it_l_m_pid = l_m_pid.insert(pair<const SharedPairID*, pair<double, bool>>(&spid_succ, pair<double, bool>(-1.0, false)));

			if (it_l_m_pid.second) it_l_m_pid.first->second.second = (w_pid.all_pids.isEmpty() || (getNumberSharedPairID(spid_succ, w_pid.all_pids, opt.min_cov_vertices) >= opt.min_cov_vertices));

			if (ud->getSharedPids(um_start.strand, succ.getMappedHead().toString()[opt.k-1]) && it_l_m_pid.first->second.second) {

				// Terminal
				{
					if (!um_e.isEmpty && (succ.isSameReferenceUnitig(um_e)) && (um_e.strand == succ.strand)){

						Path<UnitigData> path(i_t.p);
		    			const_UnitigMap<UnitigData> succ_pref(succ);

			    		if (succ_pref.strand) {

			    			succ_pref.dist = 0;
			    			succ_pref.len = um_e.dist + 1;
			    		}
					    else {

					    	succ_pref.dist = um_e.dist;
					    	succ_pref.len = succ_pref.size - um_e.dist - opt.k + 1;
					    }

					    path.extend(succ_pref); // Score is fake and will be changed later

					    if (path.length() <= max_len_path) { //Terminal path is not too large

							const pair<double, double> scores = getScorePath(opt, path, ref, ref_len, true, w_pid/*, l_m_pid*/);

							if (scores.first >= score_t1) {

								if (scores.first > score_t1) terminal_paths.clear();

								terminal_paths.push_back(move(path));

								score_t2 = score_t1;
								score_t1 = scores.first;
							}
							else if (scores.first > score_t2) score_t2 = scores.first;
						}
					}
				}

				// Non-terminal
				{
					Path<UnitigData> path(i_t.p);

					path.extend(succ); // Score is fake and will be changed later
				
					if (i_t.l != 0) stck.push(info_traversal(path, i_t.l - 1));
					else if (succ.getSuccessors().cardinality() > 0) {

						const pair<double, double> scores = getScorePath(opt, path, ref, ref_len, false, w_pid/*, l_m_pid*/);

						if (scores.first >= score_nt1) {

							if (scores.first > score_nt1) non_terminal_paths.clear();

							non_terminal_paths.push_back(move(path));

							score_nt2 = score_nt1;
							score_nt1 = scores.first;
						}
						else if (scores.first > score_nt2) score_nt2 = scores.first;
					}
				}
			}
		}
	}

	for (auto& path : terminal_paths) {

		path.setQuality(getScorePath(opt, path, ref, ref_len, score_t1, score_t2));

		for (const auto& um : path.toVector()) {

			if (um.getData()->isShortCycle()) {

				const unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>::const_iterator it = l_m_pid.find(&um.getData()->getPairID());

				if (it != l_m_pid.end()) m_pid.insert(*it);
			}
		}
	}

	for (auto& path : non_terminal_paths) {

		path.setQuality(getScorePath(opt, path, ref, ref_len, score_nt1, score_nt2));

		for (const auto& um : path.toVector()) {

			if (um.getData()->isShortCycle()) {

				const unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>::const_iterator it = l_m_pid.find(&um.getData()->getPairID());

				if (it != l_m_pid.end()) m_pid.insert(*it);
			}
		}
	}

	return {score_t1, score_nt1};
}

pair<double, double> exploreSubGraphLong(	const Correct_Opt& opt, const WeightsPairID& w_pid, const char* ref, const size_t ref_len, const size_t max_len_path,
											const const_UnitigMap<UnitigData>& um, const const_UnitigMap<UnitigData>& um_e, 
											vector<Path<UnitigData>>& terminal_paths, vector<Path<UnitigData>>& non_terminal_paths, const uint64_t hap_id,
											unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid) {

	const size_t max_len_subpath = opt.k * opt.large_k_factor;

	double score_t1 = 0.0, score_nt1 = 0.0;
	double score_t2 = 0.0, score_nt2 = 0.0;

	stack<info_traversal> stck;

	unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr> l_m_pid = m_pid;

	stck.push(info_traversal());

	while (!stck.empty()){

		const info_traversal i_t = stck.top();
		const const_UnitigMap<UnitigData>& um_start = (i_t.p.length() == 0) ? um : i_t.p.back();
		const UnitigData* ud = um_start.getData();

		stck.pop();

		for (const auto& succ : um_start.getSuccessors()){

			const UnitigData* ud_succ = succ.getData();
			const SharedPairID& spid_succ = ud_succ->getPairID();

			const pair<unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>::iterator, bool> it_l_m_pid = l_m_pid.insert(pair<const SharedPairID*, pair<double, bool>>(&spid_succ, pair<double, bool>(-1.0, false)));

			if (it_l_m_pid.second) it_l_m_pid.first->second.second = (w_pid.all_pids.isEmpty() || (getNumberSharedPairID(spid_succ, w_pid.all_pids, opt.min_cov_vertices) >= opt.min_cov_vertices));

			if (ud->getSharedPids(um_start.strand, succ.getMappedHead().toString()[opt.k-1]) && it_l_m_pid.first->second.second) {

				// Terminal
				{
					if (!um_e.isEmpty && (succ.isSameReferenceUnitig(um_e)) && (um_e.strand == succ.strand)){

						Path<UnitigData> path(i_t.p);
		    			const_UnitigMap<UnitigData> succ_pref(succ);

			    		if (succ_pref.strand) {

			    			succ_pref.dist = 0;
			    			succ_pref.len = um_e.dist + 1;
			    		}
					    else {

					    	succ_pref.dist = um_e.dist;
					    	succ_pref.len = succ_pref.size - um_e.dist - opt.k + 1;
					    }

					    path.extend(succ_pref); // Score is fake and will be changed later

					    if (path.length() <= max_len_path) { //Terminal path is not too large

							const pair<double, double> scores = getScorePath(opt, path, ref, ref_len, true, w_pid/*, l_m_pid*/);

							if (scores.first >= score_t1) {

								if (scores.first > score_t1) terminal_paths.clear();

								terminal_paths.push_back(move(path));

								score_t2 = score_t1;
								score_t1 = scores.first;
							}
							else if (scores.first > score_t2) score_t2 = scores.first;
						}
					}
				}
				
				// Non-terminal
				{
					Path<UnitigData> path(i_t.p);

					path.extend(succ);  // Score is fake and will be changed later

					if (path.length() < max_len_subpath) stck.push(info_traversal(path));
					else if (succ.getSuccessors().cardinality() > 0) {

						const pair<double, double> scores = getScorePath(opt, path, ref, ref_len, false, w_pid/*, l_m_pid*/);

						if (scores.first >= score_nt1) {

							if (scores.first > score_nt1) non_terminal_paths.clear();

							non_terminal_paths.push_back(move(path));

							score_nt2 = score_nt1;
							score_nt1 = scores.first;
						}
						else if (scores.first > score_nt2) score_nt2 = scores.first;
					}
				}
			}
		}
	}

	for (auto& path : terminal_paths) {

		path.setQuality(getScorePath(opt, path, ref, ref_len, score_t1, score_t2));

		for (const auto& um : path.toVector()) {

			if (um.getData()->isShortCycle()) {

				const unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>::const_iterator it = l_m_pid.find(&um.getData()->getPairID());

				if (it != l_m_pid.end()) m_pid.insert(*it);
			}
		}
	}

	for (auto& path : non_terminal_paths) {

		path.setQuality(getScorePath(opt, path, ref, ref_len, score_nt1, score_nt2));

		for (const auto& um : path.toVector()) {

			if (um.getData()->isShortCycle()) {

				const unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>::const_iterator it = l_m_pid.find(&um.getData()->getPairID());

				if (it != l_m_pid.end()) m_pid.insert(*it);
			}
		}
	}

	return {score_t1, score_nt1};
}

string getScorePath(const Correct_Opt& opt, const Path<UnitigData>& path, const char* ref, const size_t ref_len, const double score_best, const double score_second_best) {

	const double score_comp = score_best * ((score_best == 0.0) ? 0.0 : (1.0-(score_second_best/score_best)));
	const string path_str = path.toString();

	const EdlibAlignConfig config_ambiguous = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, edlib_iupac_alpha, sz_edlib_iupac_alpha);
	EdlibAlignResult align = edlibAlign(path_str.c_str(), path_str.length(), ref, ref_len, config_ambiguous);

	char* cigar = edlibAlignmentToCigar(align.alignment, align.alignmentLength, EDLIB_CIGAR_STANDARD);

	const size_t cigar_len = strlen(cigar);

	const char c_best_score = getQual(score_best);

	string qual_out(path_str.length(), getQual(score_comp, opt.out_qual));

	size_t cigar_pos = 0;
	size_t prev_cigar_pos = 0;
	size_t query_pos = 0;
	size_t ref_pos = 0;

	while (cigar_pos != cigar_len){

		if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

			const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

			if (cigar[cigar_pos] == 'M') { //match

				for (size_t i = 0; i < cigar_l; ++i) {

					if (path_str[query_pos + i] == ref[ref_pos + i]) qual_out[query_pos + i] = c_best_score;
				}

				query_pos += cigar_l;
				ref_pos += cigar_l;
			}
			else if ((cigar[cigar_pos] == 'I') || (cigar[cigar_pos] == 'S')) query_pos += cigar_l; //insertion or soft-clipping
			else if (cigar[cigar_pos] == 'D') ref_pos += cigar_l;  //deletion

			prev_cigar_pos = cigar_pos + 1;
		}

		++cigar_pos;
	}

	free(cigar);
	edlibFreeAlignResult(align);

	return qual_out;
}

pair<double, double> getScorePath(	const Correct_Opt& opt, const Path<UnitigData>& path, const char* ref, const size_t ref_len, const bool terminal,
									const WeightsPairID& w_pid, unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>& m_pid){

	double score = 0.0;
	double score_align = 0.0;
	double score_pids = 0.0;

	if (path.length() != 0){

		const Path<UnitigData>::PathOut path_str_um = path.toStringVector();
		const string& path_str = path_str_um.toString();
		const vector<const_UnitigMap<UnitigData>>& v_um = path_str_um.toVector();

		// Compute alignment score
		{
			EdlibAlignResult align;

			if (terminal) {

				const EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));

				align = edlibAlign(path_str.c_str(), path_str.length(), ref, ref_len, config);
				score_align = 1.0 - (static_cast<double>(align.editDistance) / path_str.length());
			}
			else {

				const EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));

				if (path_str.length() >= ref_len){

					align = edlibAlign(ref, ref_len, path_str.c_str(), path_str.length(), config);
					score_align = 1.0 - (static_cast<double>(align.editDistance) / ref_len);
				}
				else {

					const size_t l_ref_len = min(ref_len, static_cast<size_t>(path_str.length() * (1.0 + opt.weak_region_len_factor)));

					align = edlibAlign(path_str.c_str(), path_str.length(), ref, l_ref_len, config);
					score_align = 1.0 - (static_cast<double>(align.editDistance) / path_str.length());
				}
			}

			edlibFreeAlignResult(align);

			score_align = min(max(score_align, 0.0), 1.0);
		}

		if (w_pid.all_pids.cardinality() < opt.min_cov_vertices) score = 0.0;
		else  { // Final score is max. bounded by the alignment or PairID scores

			for (size_t i = 0; i < v_um.size(); ++i){

				const const_UnitigMap<UnitigData>& um = v_um[i];
				const SharedPairID& pids = um.getData()->getPairID();

				if (!pids.isEmpty()) {

					const pair<unordered_map<const SharedPairID*, pair<double, bool>, HashSharedPairIDptr>::iterator, bool> it_m_pid = m_pid.insert(pair<const SharedPairID*, pair<double, bool>>(&pids, pair<double, bool>(-1.0, false)));

					if (it_m_pid.second || (it_m_pid.first->second.first == -1.0)){

						const double shared_weighted = static_cast<double>(getNumberSharedPairID(pids, w_pid.weighted_pids));
						const double shared_noWeight = static_cast<double>(getNumberSharedPairID(pids, w_pid.noWeight_pids));

						it_m_pid.first->second.first = (shared_weighted * w_pid.weight) + shared_noWeight;
						it_m_pid.first->second.second = (static_cast<size_t>(shared_weighted + shared_noWeight) >= opt.min_cov_vertices);
					}

					score_pids += it_m_pid.first->second.first;
				}
			}

			score_pids /= (w_pid.sum_pids_weights * v_um.size());

			// Conflation
			//if (score_pids != 0.0) {
			//
			//	score = score_pids * score_align; // Conflate the two probabilities (1/2)
			//	score /= score + ((1.0 - score_pids) * (1.0 - score_align)); // Conflate the two probabilities (2/2)
			//
			//	score = (std::isnan(score) ? 0.0 : score); // Conflation can return nan (not a number) -> replace by 0.0
			//	score = min(max(score, 0.0), 1.0); // Make sure score is >= 0.0 and <= 1.0 because of double precision side effect (1.00000...0001)
			//}
			//else score = 0.0;

			//Mean
			score = (score_pids + score_align) / 2.0;
		}
	}

	return pair<double, double>(score, score_align);
}

pair<double, double> getScorePath(const Correct_Opt& opt, const Path<UnitigData>& path, const char* ref, const size_t ref_len, const bool terminal, const WeightsPairID& w_pid) {

	double score = 0.0;

	if (path.length() != 0){

		const string& path_str = path.toString();

		// Compute alignment score
		EdlibAlignResult align;

		if (terminal) {

			const EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));

			align = edlibAlign(path_str.c_str(), path_str.length(), ref, ref_len, config);
			score = 1.0 - (static_cast<double>(align.editDistance) / path_str.length());
		}
		else {

			const EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));

			if (path_str.length() >= ref_len){

				align = edlibAlign(ref, ref_len, path_str.c_str(), path_str.length(), config);
				score = 1.0 - (static_cast<double>(align.editDistance) / ref_len);
			}
			else {

				const size_t l_ref_len = min(ref_len, static_cast<size_t>(path_str.length() * (1.0 + opt.weak_region_len_factor)));

				align = edlibAlign(path_str.c_str(), path_str.length(), ref, l_ref_len, config);
				score = 1.0 - (static_cast<double>(align.editDistance) / path_str.length());
			}
		}

		edlibFreeAlignResult(align);

		score = min(max(score, 0.0), 1.0);
	}

	return {score, score};
}

vector<Path<UnitigData>> selectMostContiguous(const vector<Path<UnitigData>>& v_paths, const WeightsPairID& w_pid) {

	double score_contiguity = 0.0;

	vector<Path<UnitigData>> v_paths_out;

	for (const auto& path : v_paths) {

		const vector<const_UnitigMap<UnitigData>> v_um = path.toVector();

		PairID pid_weighted = (v_um.front().getData()->getPairID().toPairID() & w_pid.weighted_pids);
		PairID pid_no_weights = (v_um.front().getData()->getPairID().toPairID() & w_pid.noWeight_pids);

		for (const auto& um : v_um){

			if (!pid_weighted.isEmpty()) pid_weighted &= (um.getData()->getPairID().toPairID() & w_pid.weighted_pids);
			if (!pid_no_weights.isEmpty()) pid_no_weights &= (um.getData()->getPairID().toPairID() & w_pid.noWeight_pids);
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

		PairID pid = v_um.front().getData()->getPairID().toPairID();

		for (const auto& um : v_um){

			if (pid.cardinality() >= min_cov) pid &= um.getData()->getPairID().toPairID();
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

bool isValidSNPcandidate(	local_graph_traversal& lgt_fw, local_graph_traversal& lgt_bw,
							const const_UnitigMap<UnitigData>& um_a, const const_UnitigMap<UnitigData>& um_b,
							const size_t min_cov, const size_t limit_sz_stack){

	auto exploreLocalGraph = [min_cov, limit_sz_stack](local_graph_traversal& lgt, const const_UnitigMap<UnitigData>& um_a, const const_UnitigMap<UnitigData>& um_b) {

		const SharedPairID& spid_a = um_a.getData()->getPairID();
		const SharedPairID& spid_b = um_b.getData()->getPairID();

		const size_t k = um_a.getGraph()->getK();

		if ((spid_a.cardinality() < min_cov) || (spid_b.cardinality() < min_cov)) return false;

		if (lgt.m_km.empty()) {

			lgt.q_um.push(um_a);
			lgt.m_km.insert({um_a.strand ? um_a.getUnitigHead() : um_a.getUnitigTail().twin(), &spid_a});
		}
		else if (lgt.m_km.size() >= limit_sz_stack) return true;

		while (!lgt.q_um.empty()) {

			const const_UnitigMap<UnitigData> um = lgt.q_um.front();

			lgt.q_um.pop();

			for (const auto& um_neigh : um.getSuccessors()){

				if (um.getData()->getSharedPids(um.strand, um_neigh.getMappedHead().toString()[k-1])) {

					const SharedPairID& spid = um_neigh.getData()->getPairID();

					if (lgt.m_km.insert({um_neigh.getMappedHead(), &spid}).second){ // Unitig hasn't been visited so far

						if (getNumberSharedPairID(spid, spid_a, min_cov) >= min_cov){ // Unitig shares enough reads with start unitig

							if (getNumberSharedPairID(spid, spid_b, min_cov) >= min_cov) return true;
							else lgt.q_um.push(um_neigh);
						}
					}
				}
			}

			if (lgt.m_km.size() >= limit_sz_stack) return true;
		}

		return false;
	};

	const SharedPairID& spid_b = um_b.getData()->getPairID();

	bool isValidFW = false, isValidBW = false;

	{
		for (const auto& it : lgt_fw.m_km) {

			if (getNumberSharedPairID(*it.second, spid_b, min_cov) >= min_cov) {

				isValidFW = true;
				break;
			}
		}

		if (!isValidFW) isValidFW = exploreLocalGraph(lgt_fw, um_a, um_b);
	}

	if (isValidFW) {

		for (const auto& it : lgt_bw.m_km) {

			if (getNumberSharedPairID(*it.second, spid_b, min_cov) >= min_cov) {

				isValidBW = true;
				break;
			}
		}

		if (!isValidBW) {

			const_UnitigMap<UnitigData> um_a_rev = um_a;
			const_UnitigMap<UnitigData> um_b_rev = um_b;

			um_a_rev.strand = !um_a.strand;
			um_b_rev.strand = !um_b.strand;

			isValidBW = exploreLocalGraph(lgt_bw, um_a_rev, um_b_rev);
		}
	}

	return (isValidFW && isValidBW);
}

vector<Path<UnitigData>> fixRepeats(const Correct_Opt& opt, const vector<Path<UnitigData>>& v_path, const char* ref, const size_t ref_len) {

	vector<Path<UnitigData>> v_path_out;

	unordered_map<Kmer, vector<Path<UnitigData>>, KmerHash> m_cycles;

	for (auto path : v_path) {

		vector<const_UnitigMap<UnitigData>> v_um_path = path.toVector();
		string s_qual = path.toQualityString();

		EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		EdlibAlignResult align = edlibAlign(path.toString().c_str(), path.length(), ref, ref_len, config);

		int64_t editDistance = align.editDistance;

		edlibFreeAlignResult(align);

	    auto evaluatePath = [&](const Path<UnitigData>& repeat, const size_t pos_v_um, const bool strand) {

	    	Path<UnitigData> path_ext;

	    	size_t len_prefix = 0;

			const vector<const_UnitigMap<UnitigData>> v_um_repeat = repeat.toVector();

			for (size_t j = 0; j < pos_v_um; ++j) {

				path_ext.extend(v_um_path[j]); // Copy prefix
				len_prefix += v_um_path[j].len;
			}

			for (size_t j = 0; j < v_um_repeat.size(); ++j) path_ext.extend(v_um_repeat[j]); // Copy repeat
			for (size_t j = pos_v_um+1; j < v_um_path.size(); ++j) path_ext.extend(v_um_path[j]); // Copy suffix

			string l_s_qual = s_qual;

			l_s_qual.replace(len_prefix, v_um_path[pos_v_um].len + opt.k - 1, string(repeat.length(), getQual(1.0)), 0, repeat.length());

			path_ext.setQuality(l_s_qual);

			EdlibAlignConfig config = edlibNewAlignConfig(editDistance, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
			EdlibAlignResult align = edlibAlign(path_ext.toString().c_str(), path_ext.length(), ref, ref_len, config);

			const int64_t repeat_editDistance = align.editDistance;

			edlibFreeAlignResult(align);

			return pair<Path<UnitigData>, int64_t>(path_ext, repeat_editDistance);
	    };

		for (size_t i = 0; i < v_um_path.size(); ++i){

			const const_UnitigMap<UnitigData>& um_path = v_um_path[i];

			if (um_path.getData()->isShortCycle()) {

			    Path<UnitigData> best_path_ext;

			    const_UnitigMap<UnitigData> um_start(um_path); // Create a non-const local copy of unitig given in parameter
			    const_UnitigMap<UnitigData> um_end(um_path); // Create a non-const local copy of unitig given in parameter

			    pair<unordered_map<Kmer, vector<Path<UnitigData>>, KmerHash>::iterator, bool> m_it = {m_cycles.end(), false};

			    // Set prefix and suffix correctly
			    {
					um_start.len = um_path.size - um_path.dist - opt.k + 1;

	    			um_end.dist = 0;
	    			um_end.len = um_path.dist + um_path.len;

	    			// If start/end unitig is in bw (rev-comp) direction, we rev-comp the cycle as cycles are always stored in unitigs in fw direction
					um_start.strand = true;
			    	um_end.strand = true;
    			}

    			if ((i == 0) || (i == (v_um_path.size() - 1))) {

					const vector<const char*> v_compact_cycles = um_path.getData()->getCompactCycles();

					for (const auto& str_compact_cycle : v_compact_cycles) {

						const Path<UnitigData> path_extended(um_start, string(str_compact_cycle), um_end);

						pair<Path<UnitigData>, int64_t> p_eval = evaluatePath(um_path.strand ? path_extended : path_extended.rev_comp(), i, um_path.strand);

					    if ((p_eval.second >= 0) && (p_eval.second < editDistance)) {

							editDistance = p_eval.second;
							best_path_ext = move(p_eval.first);
					    }
					}
    			}
    			else {

    				const Kmer head = um_path.strand ? um_path.getUnitigHead() : um_path.getUnitigTail().twin();

	    			m_it = m_cycles.insert({head, vector<Path<UnitigData>>()});

	    			if (m_it.second) {

	    				const vector<const char*> v_compact_cycles = um_path.getData()->getCompactCycles();

						for (const auto& str_compact_cycle : v_compact_cycles) {

							Path<UnitigData> path_extended(um_start, string(str_compact_cycle), um_end);

							if (!um_path.strand) path_extended = path_extended.rev_comp();

							m_it.first->second.push_back(path_extended);

							pair<Path<UnitigData>, int64_t> p_eval = evaluatePath(path_extended, i, um_path.strand);

						    if ((p_eval.second >= 0) && (p_eval.second < editDistance)) {

								editDistance = p_eval.second;
								best_path_ext = move(p_eval.first);
						    }
						}
					}
					else {

						for (const auto& path_compact_cycle : m_it.first->second) {

							pair<Path<UnitigData>, int64_t> p_eval = evaluatePath(path_compact_cycle, i, um_path.strand);

						    if ((p_eval.second >= 0) && (p_eval.second < editDistance)) {

								editDistance = p_eval.second;
								best_path_ext = move(p_eval.first);
						    }
						}
					}
				}

				if (best_path_ext.length() != 0) { // Found a better aligning path

					const size_t diff_v_um = best_path_ext.size() - path.size();

					path = move(best_path_ext);

					v_um_path = path.toVector();
					s_qual = path.toQualityString();

					i += diff_v_um - 1;
				}
				else {

					for (size_t j = i+1; j < v_um_path.size(); ++j){

						if (um_path.isSameReferenceUnitig(v_um_path[j])) ++i;
						else break;
					}
				}
				/*else if (m_it.second) { // Skip the next UnitigMap in the sequence if they belong to cycles we have already tested (and failed to improve the alignment)

					unordered_set<Kmer, KmerHash> us_km;

					for (const auto& path_compact_cycle : m_it.first->second) {

						const vector<const_UnitigMap<UnitigData>> v_um_repeat = path_compact_cycle.toVector();

						for (const auto& um : v_um_repeat) {

							const Kmer head = um.strand ? um.getUnitigHead() : um.getUnitigTail().twin();

							us_km.insert(head);
						}
					}

					for (size_t j = i+1; j < v_um_path.size(); ++j){

						const Kmer head = v_um_path[j].strand ? v_um_path[j].getUnitigHead() : v_um_path[j].getUnitigTail().twin();

						if (us_km.find(head) != us_km.end()) ++i;
						else break;
					}
				}*/
			}
		}

		v_path_out.push_back(path);
	}

	return v_path_out;
};