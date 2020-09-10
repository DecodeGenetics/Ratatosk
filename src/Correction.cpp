#include "Correction.hpp"

pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> extractSemiWeakPaths(	const Correct_Opt& opt, const string& s,
																				const TinyBloomFilter<uint32_t>& bf, const PairID& r,
																				const pair<size_t, const_UnitigMap<UnitigData>>& um_solid_start,
																				const pair<size_t, const_UnitigMap<UnitigData>>& um_solid_end, 
																				const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_weak, size_t i_weak,
																				const size_t min_cov_vertex, const bool long_read_correct) {

	auto customSort = [](const pair<Path<UnitigData>, size_t>& p1, const pair<Path<UnitigData>, size_t>& p2){

		return (p1.first.back().mappedSequenceToString() < p2.first.back().mappedSequenceToString());
	};

	pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> paths;
	vector<pair<Path<UnitigData>, size_t>> paths1, paths2;

	const bool no_end = um_solid_end.second.isEmpty;
	const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

	const size_t step_sz = opt.k;
    const size_t pos_um_solid2 = (no_end ? s.length() - opt.k : um_solid_end.first);
    const size_t len_weak_region = (pos_um_solid2 - um_solid_start.first) + opt.k;

    const size_t max_len_weak_region = long_read_correct ? opt.max_len_weak_region2 : opt.max_len_weak_region1;
	const size_t max_paths = 512;

	size_t next_weak_pos;

    bool begin = true;
    bool end = false;

    Path<UnitigData> tmp;

    if (out_qual) tmp.extend(um_solid_start.second, 1.0);
	else tmp.extend(um_solid_start.second);

	paths1.push_back({tmp, um_solid_start.first});

	while ((i_weak < v_um_weak.size()) && (v_um_weak[i_weak].first < um_solid_start.first)) ++i_weak;

	if (i_weak < v_um_weak.size()) next_weak_pos = v_um_weak[i_weak].first;

	while (!paths1.empty() && !end) {

		pair<vector<Path<UnitigData>>, bool> g_prev_path, g_prev_path_ext;

		if (i_weak < v_um_weak.size()){

			while ((i_weak < v_um_weak.size()) && (v_um_weak[i_weak].first < (pos_um_solid2 - opt.k)) && ((v_um_weak[i_weak].first < next_weak_pos) || !hasUniquePosition(v_um_weak, i_weak))) ++i_weak;
		}
		else i_weak = v_um_weak.size();

        sort(paths1.begin(), paths1.end(), customSort);

        end = ((i_weak == v_um_weak.size()) || (v_um_weak[i_weak].first >= (pos_um_solid2 - opt.k)));

        if (end){ // Next semi-solid hit is after the target solid kmer in the query

    		for (size_t i = 0; i < paths1.size(); ++i){

    			const pair<Path<UnitigData>, size_t>& p = paths1[i];
    			const size_t l_len_weak_region = (pos_um_solid2 - p.second) + opt.k;

            	if ((i == 0) || (paths1[i].first.back() != paths1[i-1].first.back())){

            		g_prev_path = {vector<Path<UnitigData>>(), false};
            		g_prev_path_ext = {vector<Path<UnitigData>>(), false};

            		if (l_len_weak_region <= max_len_weak_region){

						if (no_end){

		        			g_prev_path = explorePathsBFS(	opt, s.c_str() + p.second, l_len_weak_region, bf, r, (begin ? um_solid_start.second : p.first.back()),
			        										min_cov_vertex, long_read_correct);
		        		}
		        		else {

		        			const vector<pair<size_t, const_UnitigMap<UnitigData>>> v_end(1, um_solid_end);

		        			g_prev_path = explorePathsBFS2(	opt, s.c_str() + p.second, l_len_weak_region, bf, r, (begin ? um_solid_start.second : p.first.back()), v_end,
			        										min_cov_vertex, long_read_correct);
		        		}
	        		}
            	}

            	if (g_prev_path.second && !g_prev_path.first.empty()){

            		for (auto& p_tmp : g_prev_path.first){

            			tmp = p.first;

            			tmp.merge(p_tmp);

            			pair<Path<UnitigData>, size_t> toPush = {move(tmp), pos_um_solid2};

            			paths2.push_back(move(toPush));
            		}
            	}
            	else if (!no_end){

					if (!g_prev_path_ext.second && g_prev_path_ext.first.empty() && (l_len_weak_region <= max_len_weak_region)){

            			g_prev_path_ext = explorePathsBFS(	opt, s.c_str() + p.second, l_len_weak_region, bf, r, (begin ? um_solid_start.second : p.first.back()),
            												min_cov_vertex, long_read_correct);
            		}

	            	if (g_prev_path_ext.second && !g_prev_path_ext.first.empty()){

	            		for (const auto& p_tmp : g_prev_path_ext.first){

	            			tmp = p.first;

	            			tmp.merge(p_tmp);
	            			paths.second.push_back(move(tmp));
	            		}
	            	}
	            	else paths.second.push_back(p.first);
            	}
            	else paths.second.push_back(p.first);
            }
        }
        else {

    		for (size_t i = 0; i < paths1.size(); ++i){

    			const pair<Path<UnitigData>, size_t>& p = paths1[i];
	            const size_t l_len_weak_region = (v_um_weak[i_weak].first - p.second) + opt.k;

	            if ((i == 0) || (paths1[i].first.back() != paths1[i-1].first.back())){

		        	size_t i_weak_end = i_weak;

		        	while ((i_weak_end < v_um_weak.size()) && (v_um_weak[i_weak_end].first < (pos_um_solid2 - opt.k)) && (v_um_weak[i_weak_end].first <= v_um_weak[i_weak].first)) ++i_weak_end;

            		g_prev_path = {vector<Path<UnitigData>>(), false};
            		g_prev_path_ext = {vector<Path<UnitigData>>(), false};

	        		const vector<pair<size_t, const_UnitigMap<UnitigData>>> v_um_weak_tmp(v_um_weak.begin() + i_weak, v_um_weak.begin() + i_weak_end);

	        		if (l_len_weak_region <= max_len_weak_region){

	        			g_prev_path = explorePathsBFS2(	opt, s.c_str() + p.second, l_len_weak_region, bf, r, (begin ? um_solid_start.second : p.first.back()), v_um_weak_tmp,
	            										min_cov_vertex, long_read_correct);
            		}
	            }

            	if (g_prev_path.second && !g_prev_path.first.empty()){

            		for (const auto& p_tmp : g_prev_path.first){

            			tmp = p.first;

            			tmp.merge(p_tmp);

            			pair<Path<UnitigData>, size_t> toPush = {move(tmp), v_um_weak[i_weak].first};

            			paths2.push_back(move(toPush));
            		}
            	}
            	else {

            		if (!g_prev_path_ext.second && g_prev_path_ext.first.empty() && (l_len_weak_region <= max_len_weak_region)){

            			g_prev_path_ext = explorePathsBFS(	opt, s.c_str() + p.second, l_len_weak_region, bf, r, (begin ? um_solid_start.second : p.first.back()),
            												min_cov_vertex, long_read_correct);
            		}

            		if (g_prev_path_ext.second && !g_prev_path_ext.first.empty()){

	            		for (const auto& p_tmp : g_prev_path_ext.first){

	            			tmp = p.first;

	            			tmp.merge(p_tmp);
	            			paths.second.push_back(move(tmp));
	            		}
	            	}
	            	else paths.second.push_back(p.first);
	            }
            }

            next_weak_pos = v_um_weak[i_weak].first + opt.k;
		}

		begin = false;
        paths1 = move(paths2);

		if (!end && (paths1.size() > max_paths)) {

			size_t best_id = 0;
			size_t max_len_ref = 0;

			bool large_paths = true;

			for (const auto& p : paths1){

				max_len_ref = max(max_len_ref, p.second - um_solid_start.first + opt.k);

				if (p.first.length() < 500) large_paths = false;
			}

			if (large_paths){

				vector<const Path<UnitigData>*> v_ptr;

				for (const auto& p : paths1) v_ptr.push_back(&(p.first));

				best_id = selectBestPrefixAlignment(s.c_str() + um_solid_start.first, max_len_ref, v_ptr).first;
			}
			else {

				double best_score = 0;

				Path<UnitigData>::PathOut p_str_vum_prev = paths1[0].first.toStringVector(), p_str_vum;

			    for (size_t i = 0; i < paths1.size(); ++i){

					p_str_vum = paths1[i].first.toStringVector(p_str_vum_prev);

					const string& path_curr = p_str_vum.toString();
			    	const vector<const_UnitigMap<UnitigData>>& v_um = p_str_vum.toVector();

					EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
					EdlibAlignResult align = edlibAlign(path_curr.c_str(), path_curr.length(), s.c_str() + um_solid_start.first, max_len_ref, config);

					const double score_kmer = 1.0 - (static_cast<double>(align.editDistance) / path_curr.length());

			    	PairID r_um = v_um[0].getData()->get_readID();

			    	for (size_t j = 1; j < v_um.size(); ++j){

			    		if (!v_um[j].isSameReferenceUnitig(v_um[j-1])) r_um |= v_um[j].getData()->get_readID();
			    	}

					const double score_read = static_cast<double>(getNumberSharedPairID(r, r_um)) / r.cardinality();

					double score = (score_read * score_kmer) / ((score_read * score_kmer) + (1.0 - score_read) * (1.0 - score_kmer));

					if (std::isnan(score)) score = 0.0;

					if (score >= best_score){

						best_score = score;
						best_id = i;
					}

					p_str_vum_prev = move(p_str_vum);

					edlibFreeAlignResult(align);
			    }
			}

		    if (!paths1.empty()) paths2.push_back(paths1[best_id]);

		    paths1 = move(paths2);
		}
	}

	for (auto& p : paths1) paths.first.push_back(move(p.first));

	return paths;
}

pair<string, string> correctSequence(	const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt,
										const string& s_fw, const string& q_fw,
										const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_solid, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_um_weak,
										const bool long_read_correct, const Roaring* all_partitions) {

	const bool out_qual = static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual);

	if ((s_fw.length() <= dbg.getK()) || v_um_solid.empty() || (v_um_solid.size() == s_fw.length() - opt.k + 1)) {

		if (out_qual) return {s_fw, long_read_correct ? q_fw : string(s_fw.length(), getQual(0.0))};

		return {s_fw, string()};
	}

	const size_t seq_len = s_fw.length();

	const string s_bw(reverse_complement(s_fw));

	const size_t max_len_weak_anchors = long_read_correct ? opt.max_len_weak_region2 : opt.max_len_weak_region1;

	const double max_norm_edit_distance = (opt.weak_region_len_factor - 1.0);

	string q_bw = q_fw;

	size_t prev_pos = v_um_solid[0].first;
	size_t i_solid = 0;
	size_t i_weak = 0;

	stringstream corrected_s;
	stringstream corrected_q;

	vector<pair<size_t, const_UnitigMap<UnitigData>>> v_um_solid_rev(v_um_solid);
	vector<pair<size_t, const_UnitigMap<UnitigData>>> v_um_weak_rev(v_um_weak);

	reverse(v_um_solid_rev.begin(), v_um_solid_rev.end());
	reverse(v_um_weak_rev.begin(), v_um_weak_rev.end());
	reverse(q_bw.begin(), q_bw.end());

	for (auto& p : v_um_solid_rev){

		p.first = seq_len - p.first - opt.k;
		p.second.strand = !p.second.strand;
	}

	for (auto& p : v_um_weak_rev){

		p.first = seq_len - p.first - opt.k;
		p.second.strand = !p.second.strand;
	}

	auto getAmbiguityVector = [](const vector<const_UnitigMap<UnitigData>>& v_um, const size_t k){

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
	};

	auto correct = [&] (const string& s, const string& q,
						const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_s, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v_w,
						const size_t i_s, const size_t i_w, const ResultCorrection* rc) -> ResultCorrection {

		const bool has_end_pt = ((i_s + 1) < v_s.size());
		const double l_max_norm_edit_distance = (max_norm_edit_distance * 1.5);

        pair<size_t, const_UnitigMap<UnitigData>> um_solid1 = v_s[i_s]; // Last hit on the leftmost solid region
        pair<size_t, const_UnitigMap<UnitigData>> um_solid2 = has_end_pt ? v_s[i_s + 1] : make_pair(s.length() - opt.k, const_UnitigMap<UnitigData>()); // First hit on the rightmost solid region

        size_t len_weak_region = um_solid2.first - um_solid1.first + opt.k;
        size_t min_cov_vertex = 0xffffffffffffffffULL;

        const int64_t min_start = um_solid1.first - opt.insert_sz;
		const int64_t min_end = um_solid2.first + opt.insert_sz;

		const char* s_start = s.c_str() + um_solid1.first;
		const char* q_start = q.empty() ? nullptr : q.c_str() + um_solid1.first;

		bool pass_min_qual_score = true;

		ResultCorrection res_corrected(len_weak_region);

		string s_corrected;
		string q_corrected;

		PairID r_s, r_w, r_e;

		vector<pair<size_t, char>> v_ambiguity;
		vector<pair<size_t, const_UnitigMap<UnitigData>>> l_v_w;

		auto setUncorrected = [&](const size_t pos, const size_t len, const size_t qual){

			s_corrected = s.substr(pos, len);

			if (out_qual) q_corrected = long_read_correct ? q.substr(pos, len) : string(len_weak_region, qual);
			//q_corrected = q.substr(pos, len_weak_region);
		};

		auto addUncorrected = [&](const size_t pos, const size_t len, const size_t qual){

			s_corrected += s.substr(pos, len);

			if (out_qual) q_corrected += long_read_correct ? q.substr(pos, len) : string(len_weak_region, qual);
			//q_corrected += q.substr(pos, len_weak_region);
		};

		if (rc != nullptr) {
		
			r_s = rc->getTargetPairID();
			r_w = rc->getWeakPairID();
			r_e = rc->getSourcePairID();
		
			min_cov_vertex = rc->getMinCovVertex();
		}

		//if (!q.empty() && (len_weak_region > opt.insert_sz)){ // Read has quality scores and region to correct is more than 500 bp long
		//
		//	pass_min_qual_score = false;
		//
		//	int sum_char_window = 0;
		//	char prev_c = 0;
		//
		//	for (size_t i = um_solid1.first + opt.k; i < um_solid1.first + 2*opt.k - 1; ++i) sum_char_window += q[i];
		//
		//	const size_t min_qs_char = 33 + opt.min_qv;
		//
		//	for (size_t i = um_solid1.first + opt.k; (i < um_solid2.first - opt.k) && !pass_min_qual_score; ++i){
		//
		//		sum_char_window = sum_char_window - prev_c + q[i+opt.k-1];
		//		pass_min_qual_score = ((sum_char_window / opt.k) >= min_qs_char);
		//		prev_c = q[i];
		//	}
		//}
		
		if (pass_min_qual_score) {

			Roaring r_id_part_solid;
			Roaring r_id_neighbors;

			if (rc == nullptr){
				
				if (r_s.isEmpty()){

					size_t i_s_s = i_s;

					unordered_set<const PairID*> s_rs;

					const UnitigData* ud = v_s[i_s_s].second.getData();
					const PairID* pid_ptr = &(ud->get_readID());

					s_rs.insert(pid_ptr);
					r_id_part_solid.add(ud->getConnectedComp());

					min_cov_vertex = min(pid_ptr->cardinality() - 1, min_cov_vertex);

					while ((i_s_s > 0) && (v_s[i_s_s].first > min_start)){

						--i_s_s;

						ud = v_s[i_s_s].second.getData();
						pid_ptr = &(ud->get_readID());

						if (!v_s[i_s_s].second.isSameReferenceUnitig(v_s[i_s_s + 1].second) && s_rs.insert(pid_ptr).second){

							s_rs.insert(pid_ptr);
							r_id_part_solid.add(ud->getConnectedComp());

							min_cov_vertex = min(pid_ptr->cardinality() - 1, min_cov_vertex);
						}
					}

					for (const PairID* pid_ptr : s_rs){

						const PairID& pid = *pid_ptr;

						for (const uint32_t id : pid) r_s.add(id);
					}
				}

				if (has_end_pt && r_e.isEmpty()){

			    	size_t i_e_e = i_s + 2;

					unordered_set<const PairID*> s_re;

					const UnitigData* ud = v_s[i_s + 1].second.getData();
					const PairID* pid_ptr = &(ud->get_readID());

					s_re.insert(pid_ptr);
					r_id_part_solid.add(ud->getConnectedComp());

					min_cov_vertex = min(pid_ptr->cardinality() - 1, min_cov_vertex);

					while ((i_e_e < v_s.size()) && (v_s[i_e_e].first < min_end)){

						ud = v_s[i_e_e].second.getData();
						pid_ptr = &(ud->get_readID());

						if (!v_s[i_e_e].second.isSameReferenceUnitig(v_s[i_e_e - 1].second) && s_re.insert(pid_ptr).second){

							r_id_part_solid.add(ud->getConnectedComp());

							min_cov_vertex = min(pid_ptr->cardinality() - 1, min_cov_vertex);
						}

						++i_e_e;
					}

					for (const PairID* pid_ptr : s_re){

						const PairID& pid = *pid_ptr;

						for (const uint32_t id : pid) r_e.add(id);
					}
				}
			}
			else r_id_neighbors = rc->getPartitions();

			if (!v_w.empty()) {
				
				const size_t start_bound = um_solid1.first + opt.insert_sz;
				const size_t end_bound = (um_solid2.first <= opt.insert_sz) ? 0 : (um_solid2.first - opt.insert_sz);

				size_t i_w_s = 0;
				size_t i_w_e = i_w - static_cast<size_t>(i_w >= v_w.size());

				vector<pair<size_t, const_UnitigMap<UnitigData>>>::const_iterator it_front, it_back;

				// Establish boundaries of semi-solid k-mers we want to use
				while ((i_w_e > 0) && (v_w[i_w_e].first > min_start)) --i_w_e;

				i_w_e += static_cast<size_t>(((i_w_e < v_w.size()) && (v_w[i_w_e].first < min_start)));
				i_w_s = i_w_e;

				while ((i_w_e < v_w.size()) && (v_w[i_w_e].first < min_end)) ++i_w_e;

				it_front = v_w.begin() + i_w_s;
				it_back = v_w.begin() + i_w_e;

				if (all_partitions != nullptr){

					if (r_id_neighbors.cardinality() == 0){

						// Add partitions from unique semi-solid k-mers to to list of solid partitions
						while ((i_w_s < v_w.size()) && (v_w[i_w_s].first < min_end)){

					    	if (hasUniquePosition(v_w, i_w_s)) r_id_part_solid.add(v_w[i_w_s].second.getData()->getConnectedComp());

							++i_w_s;
						}

						// Get all partitions neighboring solid partitions
						for (const uint32_t id_part_s_solid : r_id_part_solid){
						
							for (const uint32_t id_part_neighbor : all_partitions[id_part_s_solid]) r_id_neighbors.add(id_part_neighbor);
						}
					}

					// Keep only semi-solid k-mers from list of neighboring or solid partitions
					for (; it_front != it_back; ++it_front){

						if ((it_front->first <= start_bound) || (it_front->first >= end_bound)){ // Semi solid k-mer is within range of partition radius

							if (r_id_neighbors.contains(it_front->second.getData()->getConnectedComp())) l_v_w.push_back(*it_front); // Semi solid k-mer is a neighbor
						}
						else l_v_w.push_back(*it_front);
					}
				}
				else {

					// Keep only semi-solid k-mers from list of neighboring or solid partitions
					for (; it_front != it_back; ++it_front) l_v_w.push_back(*it_front);
				}

				if (!l_v_w.empty()){ // Remove all positions having 2+ semi solid k-mers

					vector<pair<size_t, const_UnitigMap<UnitigData>>> l_v_w_tmp;

					const size_t l_v_w_sz = l_v_w.size();

					for (int64_t i = 0; i < l_v_w_sz; ++i){

						if ((i != 0) && (l_v_w[i].first == l_v_w[i-1].first)) l_v_w[i].second.isEmpty = true;
						else if ((i != l_v_w_sz-1) && (l_v_w[i].first == l_v_w[i+1].first)) l_v_w[i].second.isEmpty = true;
					}

					for (auto& p : l_v_w){

						if (!p.second.isEmpty) l_v_w_tmp.push_back(move(p));
					}

					l_v_w = move(l_v_w_tmp);
				}
			}

			if ((rc == nullptr) && r_w.isEmpty() && !l_v_w.empty()) {

				unordered_set<const PairID*> s_rw;

				for (const auto& p : l_v_w){

			    	const PairID* pid_ptr = &(p.second.getData()->get_readID());

			    	if (s_rw.insert(pid_ptr).second) min_cov_vertex = min(pid_ptr->cardinality() - 1, min_cov_vertex);
				}

				if (!s_rw.empty()){

					for (const PairID* pid_ptr : s_rw){

						const PairID& pid = *pid_ptr;

						for (const uint32_t id : pid) r_w.add(id);
					}
				}
			}

			if (rc == nullptr){

				res_corrected.setSourcePairID(r_s);
				res_corrected.setWeakPairID(r_w);
				res_corrected.setTargetPairID(r_e);
				res_corrected.setPartitions(r_id_neighbors);
				res_corrected.setMinCovVertex(min_cov_vertex);
			}

			const size_t cov_vertex = (long_read_correct || (min_cov_vertex == 0xffffffffffffffffULL)) ? opt.min_cov_vertices : max(static_cast<size_t>((min_cov_vertex + 1) * 0.1), opt.min_cov_vertices);

			PairID r;

			if (has_end_pt){

				if (r_w.isEmpty()) r = r_s & r_e;
				else if (long_read_correct) r = r_s & r_w & r_e;
				else r = (r_s & r_w) | (r_w & r_e);
			}
			else if (!r_w.isEmpty()) r = r_s & r_w;

			if (r.cardinality() >= cov_vertex) r_s = move(r);
			else {

				r_s |= r_w;
				r_s |= r_e;
			}

			r_w.clear();
			r_e.clear();

			TinyBloomFilter<uint32_t> bf_s(r_s.cardinality(), 14);

			for (const uint32_t id : r_s) bf_s.insert(id);

	    	pair<vector<Path<UnitigData>>, vector<Path<UnitigData>>> paths1 = extractSemiWeakPaths(opt, s, bf_s, r_s, um_solid1, um_solid2, l_v_w, 0, cov_vertex, long_read_correct);

	    	if (paths1.first.empty()){

	    		size_t i_w_s = 0;
	    		size_t i_w_e = 0;

	        	while (paths1.first.empty() && !paths1.second.empty()){

	    			const pair<int, int> align = selectBestPrefixAlignment(s.c_str() + um_solid1.first, len_weak_region, paths1.second, l_max_norm_edit_distance);

	    			if (align.first == -1) break;

	    			bool next_pos_valid = false;

	    			size_t next_pos = um_solid1.first + align.second + opt.k;
	    			size_t best_id = i_w_s;

	    			while (!next_pos_valid){

						while ((i_w_s < l_v_w.size()) && (l_v_w[i_w_s].first < next_pos)) ++i_w_s;

						const size_t len_to_next_seed = l_v_w[i_w_s].first - um_solid1.first;

						next_pos_valid = (i_w_s < l_v_w.size()) && (l_v_w[i_w_s].first < um_solid2.first - opt.k) && ((l_v_w[i_w_s].first - um_solid1.first) < max_len_weak_anchors);

						if (next_pos_valid){

							i_w_e = i_w_s + 1;

							while ((i_w_e < l_v_w.size()) && (l_v_w[i_w_e].first == l_v_w[i_w_s].first)) ++i_w_e;

			            	if (i_w_e - i_w_s != 1){

		            			next_pos_valid = false;
		            			next_pos = l_v_w[i_w_s].first + 1;
			            	}
			            	else best_id = i_w_s;
						}
						else break;
					}

					if (!next_pos_valid) break;

					const Path<UnitigData>::PathOut path_out = paths1.second[align.first].toStringVector();
					const vector<pair<size_t, char>> v_amb = getAmbiguityVector(path_out.toVector(), opt.k);

					for (const auto p : v_amb) v_ambiguity.push_back({s_corrected.length() + p.first, p.second});

					s_corrected += path_out.toString() + s.substr(um_solid1.first + align.second + 1, l_v_w[i_w_s].first - um_solid1.first - align.second - 1);

					if (out_qual){

						q_corrected += path_out.toQualityString();

						if (long_read_correct) q_corrected += q.substr(um_solid1.first + align.second + 1, l_v_w[i_w_s].first - um_solid1.first - align.second - 1);
						else q_corrected += string(l_v_w[i_w_s].first - um_solid1.first - align.second - 1, getQual(0.0));
						//q_corrected += q.substr(um_solid1.first + align.second + 1, l_v_w[i_w_s].first - um_solid1.first - align.second - 1);
					}

					res_corrected.addCorrectedPosOldSeq(um_solid1.first - v_s[i_s].first, um_solid1.first + align.second + 1 - v_s[i_s].first);

	            	um_solid1 = l_v_w[best_id];
	            	len_weak_region = um_solid2.first - um_solid1.first + opt.k;

	            	paths1 = extractSemiWeakPaths(opt, s, bf_s, r_s, um_solid1, um_solid2, l_v_w, i_w_e, cov_vertex, long_read_correct);
	    		}

	        	if (!paths1.first.empty()){

	        		const pair<int, int> p_align = selectBestAlignment(paths1.first, s.c_str() + um_solid1.first, len_weak_region);
	        		const Path<UnitigData>::PathOut path_out = paths1.first[p_align.first].toStringVector();
	        		const vector<pair<size_t, char>> v_amb = getAmbiguityVector(path_out.toVector(), opt.k);

	        		for (const auto p : v_amb) v_ambiguity.push_back({s_corrected.length() + p.first, p.second});

	        		s_corrected += path_out.toString();

	        		if (out_qual) q_corrected += path_out.toQualityString();

					res_corrected.addCorrectedPosOldSeq(um_solid1.first - v_s[i_s].first, um_solid2.first - v_s[i_s].first + opt.k);
	        	}
	        	else if (!paths1.second.empty()){

	        		const pair<int, int> align = selectBestPrefixAlignment(s.c_str() + um_solid1.first, len_weak_region, paths1.second, l_max_norm_edit_distance);

	        		if (align.first == -1) addUncorrected(um_solid1.first, len_weak_region, getQual(0.0));
	        		else {

	        			const Path<UnitigData>::PathOut path_out = paths1.second[align.first].toStringVector();
	        			const vector<pair<size_t, char>> v_amb = getAmbiguityVector(path_out.toVector(), opt.k);

			        	for (const auto p : v_amb) v_ambiguity.push_back({s_corrected.length() + p.first, p.second});

			        	s_corrected += path_out.toString() + s.substr(um_solid1.first + align.second + 1, len_weak_region - align.second - 1);

			        	if (out_qual){

			        		q_corrected += path_out.toQualityString();

			        		if (long_read_correct) q_corrected += q.substr(um_solid1.first + align.second + 1, len_weak_region - align.second - 1);
			        		else q_corrected += string(len_weak_region - align.second - 1, getQual(0.0));
			        		//q_corrected += q.substr(um_solid1.first + align.second + 1, len_weak_region - align.second - 1);
			        	}

						res_corrected.addCorrectedPosOldSeq(um_solid1.first - v_s[i_s].first, um_solid1.first + align.second + 1 - v_s[i_s].first);
					}
	    		}
	    		else if (!s_corrected.empty()) addUncorrected(um_solid1.first, len_weak_region, getQual(0.0));
	    		else setUncorrected(v_s[i_s].first, len_weak_region, getQual(0.0));
	    	}
	    	else {

	    		const pair<int, int> p_align = selectBestAlignment(paths1.first, s.c_str() + um_solid1.first, len_weak_region);
	    		const Path<UnitigData>::PathOut path_out = paths1.first[p_align.first].toStringVector();
	    		const vector<pair<size_t, char>> v_amb = getAmbiguityVector(path_out.toVector(), opt.k);

	    		for (const auto p : v_amb) v_ambiguity.push_back(p);

	    		s_corrected = path_out.toString();

	    		if (out_qual) q_corrected = path_out.toQualityString();

				res_corrected.addCorrectedPosOldSeq(0, len_weak_region);
	    	}
    	}
    	else setUncorrected(v_s[i_s].first, len_weak_region, static_cast<char>(126));

    	if (!q_corrected.empty()){

	    	if (long_read_correct){

	    		// Make sure that every q-score value which was set to "real bad" (126) during the first pass is replaced by the minimum qual score (33)
	    		std::replace(q_corrected.begin(), q_corrected.end(), static_cast<char>(126), getQual(0.0));

	    		// First k q-score values are copied back from input quality sequence as they correspond to a solid match
	    		if (!q.empty()){

	    			for (size_t i = 0; i < opt.k; ++i) q_corrected[i] = q[i];
	    		}
	    	}
	    	else {

	    		// First k q-score values are maximum as they correspond to a solid match
	    		for (size_t i = 0; i < opt.k; ++i) q_corrected[i] = getQual(1.0);
	    	}
    	}

    	fixAmbiguity(dbg, s_corrected, q_corrected, s_start, q_start, res_corrected.getLengthOldSequence(), v_ambiguity, opt.min_qv, true);

    	if (res_corrected.getNbCorrectedPosOldSeq() == res_corrected.getLengthOldSequence()){

	    	const Kmer end_s(s.c_str() + s.length() - opt.k);
	    	const Kmer end_s_corr(s_corrected.c_str() + s_corrected.length() - opt.k);

	    	// Region has all its bases corrected, first and last k-mers are solid -> region is corrected
	    	if (end_s == end_s_corr) res_corrected.setCorrected();
	    }
	    
	    if (!res_corrected.isCorrected()) {

    		// Corrected region might span the weak region too far
			EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0));
			EdlibAlignResult align = edlibAlign(s_start, um_solid2.first - v_s[i_s].first + opt.k, s_corrected.c_str(), s_corrected.length(), config);

			if (align.editDistance >= 0){

				size_t end_location = align.endLocations[0];

				for (size_t j = 1; j < align.numLocations; ++j){

					if (align.endLocations[j] > end_location) end_location = align.endLocations[j];
				}

				s_corrected = s_corrected.substr(0, end_location + 1);

				if (!q_corrected.empty()) q_corrected = q_corrected.substr(0, end_location + 1);
			}

			edlibFreeAlignResult(align);
    	}

    	res_corrected.setSequence(move(s_corrected));
    	res_corrected.setQuality(move(q_corrected));

		return res_corrected;
	};

    if (v_um_solid[0].first != 0){ // Read starts with a weak region (no left solid region)

    	const size_t i_solid_rev = v_um_solid_rev.size() - 1;

    	size_t i_weak_rev = v_um_weak_rev.size();

    	while ((i_weak_rev > 0) && (v_um_weak_rev[i_weak_rev - 1].first > v_um_solid_rev[i_solid_rev].first)) --i_weak_rev;

    	const ResultCorrection bw_corrected = correct(s_bw, q_fw, v_um_solid_rev, v_um_weak_rev, i_solid_rev, i_weak_rev, nullptr).reverseComplement();

    	corrected_s << bw_corrected.getSequence().substr(0, bw_corrected.getSequence().length() - opt.k);

    	if (out_qual) corrected_q << bw_corrected.getQuality().substr(0, bw_corrected.getQuality().length() - opt.k); // For now everything corrected is score 0
    }

	while (i_solid < v_um_solid.size() - 1){

	    if (v_um_solid[i_solid].first != (v_um_solid[i_solid + 1].first - 1)){ // Start of a weak region of the read

	    	bool sameUnitig = false;

	        const const_UnitigMap<UnitigData>& start_um = v_um_solid[i_solid].second;
	        const const_UnitigMap<UnitigData>& end_um = v_um_solid[i_solid + 1].second;

	    	if (long_read_correct || (v_um_solid[i_solid + 1].first < (v_um_solid[i_solid].first + opt.k))){

				// Check same unitig and same strand
	        	sameUnitig = (start_um.getUnitigHead() == end_um.getUnitigHead()) && (start_um.strand == end_um.strand);
	        	// Start and end positions are compatible
	        	sameUnitig = sameUnitig && ((start_um.strand && (start_um.dist < end_um.dist)) || (!start_um.strand && (start_um.dist > end_um.dist)));
        	}

        	if (sameUnitig) {
        		
        		const_UnitigMap<UnitigData> um_sub = v_um_solid[i_solid].second;

        		um_sub.dist = min(start_um.dist, end_um.dist);
        		um_sub.len = (start_um.strand ? (end_um.dist - start_um.dist) : (start_um.dist - end_um.dist)) + 1;

        		const string s_um_sub = um_sub.mappedSequenceToString();

        		corrected_s << s_fw.substr(prev_pos, v_um_solid[i_solid].first - prev_pos) + s_um_sub.substr(0, s_um_sub.length() - opt.k);

        		if (out_qual){

        			corrected_q << (long_read_correct ? q_fw.substr(prev_pos, v_um_solid[i_solid].first - prev_pos) : string(v_um_solid[i_solid].first - prev_pos, getQual(1.0))) +
        							string(s_um_sub.length() - opt.k, getQual(1.0));
        		}
	        }
			else if (v_um_solid[i_solid + 1].first >= (v_um_solid[i_solid].first + opt.k)){

	        	const ResultCorrection fw_corrected = correct(s_fw, q_fw, v_um_solid, v_um_weak, i_solid, i_weak, nullptr);

	        	if (fw_corrected.isCorrected()) {

	        		const size_t l_solid = v_um_solid[i_solid].first - prev_pos;
	        		const string corrected_subseq_fw = s_fw.substr(prev_pos, l_solid) + fw_corrected.getSequence();

	        		corrected_s << corrected_subseq_fw.substr(0, corrected_subseq_fw.length() - opt.k);

	        		if (out_qual){

	        			const string corrected_subqual_fw = (long_read_correct ? q_fw.substr(prev_pos, l_solid) : string(l_solid, getQual(1.0))) + fw_corrected.getQuality();

	        			corrected_q << corrected_subqual_fw.substr(0, corrected_subqual_fw.length() - opt.k);
	        		}
	        	}
	        	else {

		    		size_t i_solid_bw = v_um_solid_rev.size() - i_solid - 2;
		    		size_t i_weak_bw = v_um_weak_rev.size() - i_weak;

		    		while ((i_weak_bw > 0) && (v_um_weak_rev[i_weak_bw - 1].first > v_um_solid_rev[i_solid_bw].first)) --i_weak_bw;

		        	const ResultCorrection bw_corrected = correct(s_bw, q_bw, v_um_solid_rev, v_um_weak_rev, i_solid_bw, i_weak_bw, &fw_corrected).reverseComplement();

		        	if (bw_corrected.isCorrected()){

		        		const size_t l_solid = (s_bw.length() - v_um_solid_rev[i_solid_bw + 1].first - opt.k) - prev_pos;
		        		const string corrected_subseq_bw = s_fw.substr(prev_pos, l_solid) + bw_corrected.getSequence();

		        		corrected_s << corrected_subseq_bw.substr(0, corrected_subseq_bw.length() - opt.k);

		        		if (out_qual){

		        			const string corrected_subqual_bw = (long_read_correct ? q_fw.substr(prev_pos, l_solid) : string(l_solid, getQual(1.0))) + bw_corrected.getQuality();

		        			corrected_q << corrected_subqual_bw.substr(0, corrected_subqual_bw.length() - opt.k);
		        		}
		        	}
		        	else {

		        		string l_ref = s_fw.substr(v_um_solid[i_solid].first, v_um_solid[i_solid + 1].first - v_um_solid[i_solid].first + opt.k);
		        		
		        		pair<string, string> p_consensus = generateConsensus(&fw_corrected, &bw_corrected, l_ref, max_norm_edit_distance * 1.5);

		        		if (p_consensus.first.length() == 0) {

		        			p_consensus.first = move(l_ref);

		        			if (out_qual){

		        				if (long_read_correct) p_consensus.second = q_fw.substr(v_um_solid[i_solid].first, v_um_solid[i_solid + 1].first - v_um_solid[i_solid].first + opt.k);
		        				else p_consensus.second = string(opt.k, getQual(1.0)) + string(v_um_solid[i_solid + 1].first - v_um_solid[i_solid].first, getQual(0.0));
		        			}
		        		}

		        		const size_t l_solid = v_um_solid[i_solid].first - prev_pos;
				        const string corrected_subseq_consensus = s_fw.substr(prev_pos, l_solid) + p_consensus.first;

				        corrected_s << corrected_subseq_consensus.substr(0, corrected_subseq_consensus.length() - opt.k);

				        if (out_qual){

		        			const string corrected_subqual_consensus = (long_read_correct ? q_fw.substr(prev_pos, l_solid) : string(l_solid, getQual(1.0))) + p_consensus.second;

		        			corrected_q << corrected_subqual_consensus.substr(0, corrected_subqual_consensus.length() - opt.k);
		        		}
			        }
	        	}
	        }
			else { // No correction

	        	corrected_s << s_fw.substr(prev_pos, v_um_solid[i_solid + 1].first - prev_pos);

	        	if (out_qual){

	        		if (long_read_correct) corrected_q << q_fw.substr(prev_pos, v_um_solid[i_solid + 1].first - prev_pos);
	        		else corrected_q << string(v_um_solid[i_solid + 1].first - prev_pos, getQual(1.0));
	        	}
        	}

        	prev_pos = v_um_solid[i_solid + 1].first;
        }
        
        ++i_solid;
    }

    if (v_um_solid[v_um_solid.size() - 1].first < s_fw.length() - opt.k){ // There is a weak region at the end which has no

   		const ResultCorrection fw_corrected = correct(s_fw, q_fw, v_um_solid, v_um_weak, i_solid, i_weak, nullptr);
   		const size_t l_solid = v_um_solid[i_solid].first - prev_pos;

        corrected_s << s_fw.substr(prev_pos, l_solid) << fw_corrected.getSequence(); // Push right solid region

        if (out_qual) corrected_q << (long_read_correct ? q_fw.substr(prev_pos, l_solid) : string(l_solid, getQual(1.0))) << fw_corrected.getQuality();
    }
	else {

		corrected_s << s_fw.substr(prev_pos);

		if (out_qual) corrected_q << (long_read_correct ? q_fw.substr(prev_pos) : string(s_fw.length() - prev_pos, getQual(1.0)));
	}

    return {corrected_s.str(), corrected_q.str()};
}