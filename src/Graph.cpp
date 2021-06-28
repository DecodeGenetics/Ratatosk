#include "Graph.hpp"

pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> getSeeds(const Correct_Opt& opt, const CompactedDBG<UnitigData>& dbg,
																													const string& s, const string& q, 
																													const bool long_read_correct, const uint64_t hap_id,
																													unordered_map<Kmer, vector<const_UnitigMap<UnitigData>>, KmerHash>& m_km_um,
																													const bool fill_map) {

    auto comp_pair = [](const pair<size_t, const_UnitigMap<UnitigData>>& p1, const pair<size_t, const_UnitigMap<UnitigData>>& p2) {

        if (p1.first == p2.first) return (p1.second.mappedSequenceToString() < p2.second.mappedSequenceToString());

        return (p1.first < p2.first);
    };

    auto removeEmptyUnitigMap = [](vector<pair<size_t, const_UnitigMap<UnitigData>>>& v) {

    	vector<pair<size_t, const_UnitigMap<UnitigData>>> v2;

    	for (const auto& p : v){

    		if (!p.second.isEmpty) v2.push_back(p);
    	}

    	v = move(v2);
    };

    auto fillMap = [&](const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v) {

		for (const auto& p : v) {

			const Kmer km(s.c_str() + p.first);
			const pair<unordered_map<Kmer, vector<const_UnitigMap<UnitigData>>, KmerHash>::iterator, bool> p_it = m_km_um.insert({km, vector<const_UnitigMap<UnitigData>>()});

			vector<const_UnitigMap<UnitigData>>& l_v = p_it.first->second;

			if (p_it.second) l_v.push_back(p.second);
			else {

				bool same = false;

				for (const auto& um : l_v) same = (same || (um == p.second));

				if (!same) l_v.push_back(p.second);
			}
		}
    };

    if (s.length() <= opt.k) return {vector<pair<size_t, const_UnitigMap<UnitigData>>>(), vector<pair<size_t, const_UnitigMap<UnitigData>>>()};

    // Find all exact or inexact k-mer hits between the query and the graph
    vector<pair<size_t, const_UnitigMap<UnitigData>>> v_um, v_um_solid, v_um_weak;

    if (!m_km_um.empty()) { // Some k-mers have already been searched for on a previous round, no need to re-search those

    	string l_s = std::string(s.length(), 'N');

    	vector<pair<size_t, const_UnitigMap<UnitigData>>> l_v_um;

       	for (size_t i = 0; i < s.length() - opt.k + 1; ++i){

       		bool isValidKmer = true;

       		for (size_t j = i; (j < (i + opt.k)) && isValidKmer; ++j) isValidKmer = isDNA(s[j]);

       		if (isValidKmer) {

	    		const Kmer km(s.c_str() + i);
	    		const unordered_map<Kmer, vector<const_UnitigMap<UnitigData>>, KmerHash>::const_iterator it = m_km_um.find(km);

	    		if (it != m_km_um.end()) {

	    			for (const auto& um : it->second) l_v_um.push_back({i, um});
	    		}
	    		else l_s.replace(i, opt.k, s, i, opt.k);
	    	}
	    	else l_s.replace(i, opt.k, s, i, opt.k);
    	}

    	v_um = dbg.searchSequence(l_s, true, true, true, true, true);

    	if (fill_map) fillMap(v_um); 

    	v_um.insert(v_um.end(), l_v_um.begin(), l_v_um.end());
    }
    else {

    	v_um = dbg.searchSequence(s, true, true, true, true, true);

    	if (fill_map) fillMap(v_um); 
    }

    sort(v_um.begin(), v_um.end(), comp_pair); // Sort k-mer hits by query position and mapped k-mer similarity

    // The next piece of code separate query hits into 2 vectors:
    // - v_um_solid: contains solid and semi-solid k-mers. Solid k-mer is an exact match with the graph.
    //   Semi-solid k-mer is a non-solid k-mer with one inexact match in the graph (maximum one indel or substitution).
    // - v_um_weak: contains semi-weak k-mers. Semi-weak k-mer is a non-solid k-mer with more than one inexact match
    //   in the graph (maximum one indel or substitution).

    // Pull out solid matches
    {
    	const size_t v_um_sz = v_um.size();

	    int64_t exact_match = -1;

	    size_t i = 0;
    	size_t prev_idx = 0xffffffffffffffffULL;

	    while (i < v_um_sz){ // For every exact or inexact k-mer hit between the query and the graph

	        const size_t idx = v_um[i].first;
	        const string km_ref = s.substr(idx, opt.k);
	        const string km_mapped = v_um[i].second.mappedSequenceToString();

	        if (idx != prev_idx){ // New k-mer on the query

	            // For previous query k-mer, there was no exact match but exactly one inexact match,
	            // use that match as a seed;
	            if (exact_match != -1) v_um_solid.push_back(v_um[exact_match]);

	            prev_idx = idx;
	            exact_match = -1;
	        }

	        if ((km_ref == km_mapped) && (exact_match == -1)) exact_match = i;

	        ++i;

	        while ((i < v_um_sz) && (v_um[i].first == idx) && (v_um[i].second.mappedSequenceToString() == km_mapped)) ++i;
	    }

	    if (exact_match != -1) v_um_solid.push_back(v_um[exact_match]);
	}

	// Make sure solid matches overlap by k-1
	if (v_um_solid.size() >= 2) {

    	for (size_t i = 1; i < v_um_solid.size(); ++i){

    		if ((v_um_solid[i].first != (v_um_solid[i-1].first + 1)) && (v_um_solid[i].first < (v_um_solid[i-1].first + opt.k))) {

    			v_um_solid[i-1].second.isEmpty = true;

    			int64_t j = i-2;

    			while ((j >= 0) && (v_um_solid[j].first == (v_um_solid[j+1].first - 1)) && (v_um_solid[i].first < (v_um_solid[j].first + opt.k))){

    				v_um_solid[j--].second.isEmpty = true;
    			}
    		}
    	}

    	removeEmptyUnitigMap(v_um_solid);
	}

    // Build the list of semi-weak k-mers
    {
	    size_t i = 0, j = 0;

    	const size_t v_um_sz = v_um.size();
	    const size_t v_um_solid_sz = v_um_solid.size();

	    while ((i != v_um_sz) && (j != v_um_solid_sz)) {

	     	if (v_um[i].first < v_um_solid[j].first){

	     		if (v_um_weak.empty() || (v_um[i] != v_um_weak.back())) v_um_weak.push_back(v_um[i]);
		     	
		     	++i;
	     	}
	     	else if (v_um[i].first > v_um_solid[j].first) ++j;
	     	else ++i;
	    }

	    if (j == v_um_solid.size()){

	    	while (i != v_um_sz){

	     		v_um_weak.push_back(v_um[i]);
	     		++i;
	    	}
	    }
	}

	//Remove duplicated weak matches
	if (!v_um_weak.empty()){

		vector<pair<size_t, const_UnitigMap<UnitigData>>> v_um_weak_tmp(1, v_um_weak[0]);

		const size_t v_um_weak_sz = v_um_weak.size();

		for (size_t i = 1; i < v_um_weak_sz; ++i){

			if ((v_um_weak[i].first != v_um_weak[i-1].first) || (v_um_weak[i].second.mappedSequenceToString() != v_um_weak[i-1].second.mappedSequenceToString())) v_um_weak_tmp.push_back(v_um_weak[i]);
		}

		v_um_weak = move(v_um_weak_tmp);
	}

	// Runs of solid anchors must be consistent:
	// - Anchor i+1 must be on same unitig as anchor i or on a successor unitig
	// - In the latter case, both unitigs must share a minimum number of colors
	{
    	for (size_t i = 1; i < v_um_solid.size(); ++i){ // Make sure solid regions are color-consistent

    		if ((v_um_solid[i].first - v_um_solid[i-1].first) == 1) {

    			const_UnitigMap<UnitigData>& um_left = v_um_solid[i-1].second;
    			const_UnitigMap<UnitigData>& um_right = v_um_solid[i].second;

    			if (!um_left.isEmpty && !um_right.isEmpty && !um_left.isSameReferenceUnitig(um_right)){

    				const string s_tail_left = (um_left.strand ? um_left.getUnitigTail() : um_left.getUnitigHead().twin()).toString();
    				const string s_head_right = (um_right.strand ? um_right.getUnitigHead() : um_right.getUnitigTail().twin()).toString();

    				bool invalidAnchors = (s_tail_left.substr(1, opt.k-1) != s_head_right.substr(0, opt.k-1));

    				if (!invalidAnchors) invalidAnchors = (getNumberSharedPairID(um_left.getData()->getPairID(), um_right.getData()->getPairID(), opt.min_cov_vertices) < opt.min_cov_vertices);

    				if (invalidAnchors) {

		    			size_t i_l = i-1;
		    			size_t i_r = i+1;

		    			i_l -= (i_l != 0);

		    			while ((i_l > 0) && (v_um_solid[i_l].first == (v_um_solid[i_l+1].first - 1)) && v_um_solid[i_l].second.isSameReferenceUnitig(um_left)){

		    				v_um_solid[i_l].second.isEmpty = true;
		    				--i_l;
		    			}

		    			while ((i_r < v_um_solid.size()) && (v_um_solid[i_r].first == (v_um_solid[i_r-1].first + 1)) && v_um_solid[i_r].second.isSameReferenceUnitig(um_right)){

		    				v_um_solid[i_r].second.isEmpty = true;
		    				++i_r;
		    			}

		    			um_left.isEmpty = true;
		    			um_right.isEmpty = true;
    				}
    			}
    		}
    	}

    	removeEmptyUnitigMap(v_um_solid);
    }

    {
    	for (size_t i = 1; i < v_um_solid.size(); ++i){ // Make sure solid regions are color-consistent

    		const size_t diff_pos = v_um_solid[i].first - v_um_solid[i-1].first;

    		if ((diff_pos > 1) && (diff_pos < (opt.insert_sz/2))) {

    			size_t i_l = i-1;
    			size_t i_r = i+1;

    			SharedPairID pid_left = v_um_solid[i-1].second.getData()->getPairID();
    			SharedPairID pid_right = v_um_solid[i].second.getData()->getPairID();

    			const_UnitigMap<UnitigData> prev_um_left = v_um_solid[i-1].second;
    			const_UnitigMap<UnitigData> prev_um_right = v_um_solid[i].second;

    			i_l -= (i_l != 0);

    			while ((i_l > 0) && (v_um_solid[i_l].first == (v_um_solid[i_l+1].first - 1))){

    				if (!v_um_solid[i_l].second.isSameReferenceUnitig(prev_um_left)) {

	    				pid_left |= v_um_solid[i_l].second.getData()->getPairID();
	    				prev_um_left = v_um_solid[i_l].second;
    				}

    				--i_l;
    			}

    			while ((i_r < v_um_solid.size()) && (v_um_solid[i_r].first == (v_um_solid[i_r-1].first + 1))){

    				if (!v_um_solid[i_r].second.isSameReferenceUnitig(prev_um_right)) {

	    				pid_right |= v_um_solid[i_r].second.getData()->getPairID();
	    				prev_um_right = v_um_solid[i_r].second;
    				}

    				++i_r;
    			}

    			if (getNumberSharedPairID(pid_left, pid_right, opt.min_cov_vertices) < opt.min_cov_vertices){

    				for (++i_l; i_l < i_r; ++i_l) v_um_solid[i_l].second.isEmpty = true;

    				i = i_r;
    			}
    		}
    	}

    	removeEmptyUnitigMap(v_um_solid);
    }

    pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = {move(v_um_solid), move(v_um_weak)};

    return p;
}

void detectSNPs(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt){

	if (opt.verbose) cout << "Ratatosk::detectSNPs(): Scanning graph for k-mers with substitutions (SNPs)" << endl;

	const size_t k = dbg.getK();

    auto comp_pair = [](const pair<size_t, const_UnitigMap<UnitigData>>& p1, const pair<size_t, const_UnitigMap<UnitigData>>& p2) {

        return (p1.first < p2.first);
    };

	if (opt.nb_threads == 1){

		size_t count_vertices = 0;

	    for (auto& um : dbg){

			const string seq_ref = um.referenceUnitigToString();

			vector<pair<size_t, UnitigMap<UnitigData>>> v_um = dbg.searchSequence(seq_ref, false, false, false, true, false);

			if (!v_um.empty()){

				UnitigData* ud = um.getData();

				string seq_final = seq_ref;
				string seq_tried = seq_ref;

				unordered_set<Kmer, KmerHash> s_valid_unitigs, s_invalid_unitigs;

				local_graph_traversal lgt_fw, lgt_bw;

				sort(v_um.begin(), v_um.end(), comp_pair);

		    	for (size_t i = 0; i < v_um.size(); ++i){

		    		const pair<size_t, UnitigMap<UnitigData>>& p = v_um[i];

	    			const Kmer head = p.second.getUnitigHead();
		    		const string km_snp = p.second.mappedSequenceToString();
		    		const size_t pos_snp_km = cstrMatch(km_snp.c_str(), seq_ref.c_str() + p.first);

		    		if (pos_snp_km < k){

		    			bool nuc_a_f = false, nuc_c_f = false, nuc_g_f = false, nuc_t_f = false;
			    		bool nuc_a_t = false, nuc_c_t = false, nuc_g_t = false, nuc_t_t = false;
			    		bool nuc_a_k = false, nuc_c_k = false, nuc_g_k = false, nuc_t_k = false;

			    		getAmbiguityRev(seq_final[p.first + pos_snp_km], nuc_a_f, nuc_c_f, nuc_g_f, nuc_t_f);
			    		getAmbiguityRev(seq_tried[p.first + pos_snp_km], nuc_a_t, nuc_c_t, nuc_g_t, nuc_t_t);
			    		getAmbiguityRev(km_snp[pos_snp_km], nuc_a_k, nuc_c_k, nuc_g_k, nuc_t_k);

			    		const char cf = getAmbiguity(nuc_a_f || nuc_a_k, nuc_c_f || nuc_c_k, nuc_g_f || nuc_g_k, nuc_t_f || nuc_t_k);
			    		const char ct = getAmbiguity(nuc_a_t || nuc_a_k, nuc_c_t || nuc_c_k, nuc_g_t || nuc_g_k, nuc_t_t || nuc_t_k);

			    		if (seq_tried[p.first + pos_snp_km] != ct) { // We didn't try that SNP candidate before

			    			seq_tried[p.first + pos_snp_km] = ct;

			    			if (s_valid_unitigs.find(head) != s_valid_unitigs.end()) seq_final[p.first + pos_snp_km] = cf; // Unitig was already analyzed before and was valid
			    			else if (s_invalid_unitigs.find(head) == s_invalid_unitigs.end()) { // Unitig was not already analyzed before or it was invalid

			    				if (isValidSNPcandidate(lgt_fw, lgt_bw, um, p.second, opt.insert_sz, opt.min_cov_vertices)) {

			    					seq_final[p.first + pos_snp_km] = cf;

			    					s_valid_unitigs.insert(head);
			    				}
			    				else s_invalid_unitigs.insert(head);
			    			}
			    		}
			    	}
			    }

			    for (size_t i = 0; i < seq_final.length(); ++i){

			    	if (!isDNA(seq_final[i])) ud->add_ambiguity_char(i, seq_final[i]);
			    }

			    ud->runOptimizeAmbiguityChar();
			}

			++count_vertices;

			if (opt.verbose && ((count_vertices%1000000)==0)) cout << "Ratatosk::detectSNPs(): Processed " << count_vertices << " vertices." << endl;
		}
	}
	else {

		auto worker_function = [&](CompactedDBG<UnitigData>::iterator it_c, CompactedDBG<UnitigData>::iterator it_e) {

			size_t l_nb_fp_avoided = 0;

			while (it_c != it_e){

				const UnitigMap<UnitigData>& um = *it_c;
				const string seq_ref = um.referenceUnitigToString();

				vector<pair<size_t, UnitigMap<UnitigData>>> v_um = dbg.searchSequence(seq_ref, false, false, false, true, false);

				if (!v_um.empty()){

					UnitigData* ud = um.getData();

					unordered_set<Kmer, KmerHash> s_valid_unitigs, s_invalid_unitigs;

					local_graph_traversal lgt_fw, lgt_bw;

					string seq_final = seq_ref;
					string seq_tried = seq_ref;

					sort(v_um.begin(), v_um.end(), comp_pair);

			    	for (size_t i = 0; i < v_um.size(); ++i){

			    		const pair<size_t, UnitigMap<UnitigData>>& p = v_um[i];

		    			const Kmer head = p.second.getUnitigHead();
			    		const string km_snp = p.second.mappedSequenceToString();
			    		const size_t pos_snp_km = cstrMatch(km_snp.c_str(), seq_ref.c_str() + p.first);

			    		if (pos_snp_km < k){

			    			bool nuc_a_f = false, nuc_c_f = false, nuc_g_f = false, nuc_t_f = false;
				    		bool nuc_a_t = false, nuc_c_t = false, nuc_g_t = false, nuc_t_t = false;
				    		bool nuc_a_k = false, nuc_c_k = false, nuc_g_k = false, nuc_t_k = false;

				    		getAmbiguityRev(seq_final[p.first + pos_snp_km], nuc_a_f, nuc_c_f, nuc_g_f, nuc_t_f);
				    		getAmbiguityRev(seq_tried[p.first + pos_snp_km], nuc_a_t, nuc_c_t, nuc_g_t, nuc_t_t);
				    		getAmbiguityRev(km_snp[pos_snp_km], nuc_a_k, nuc_c_k, nuc_g_k, nuc_t_k);

				    		const char cf = getAmbiguity(nuc_a_f || nuc_a_k, nuc_c_f || nuc_c_k, nuc_g_f || nuc_g_k, nuc_t_f || nuc_t_k);
				    		const char ct = getAmbiguity(nuc_a_t || nuc_a_k, nuc_c_t || nuc_c_k, nuc_g_t || nuc_g_k, nuc_t_t || nuc_t_k);

				    		if (seq_tried[p.first + pos_snp_km] != ct) { // We didn't try that SNP candidate before

				    			seq_tried[p.first + pos_snp_km] = ct;

				    			if (s_valid_unitigs.find(head) != s_valid_unitigs.end()) seq_final[p.first + pos_snp_km] = cf; // Unitig was already analyzed before and was valid
				    			else if (s_invalid_unitigs.find(head) == s_invalid_unitigs.end()) { // Unitig was not already analyzed before or it was invalid

				    				if (isValidSNPcandidate(lgt_fw, lgt_bw, um, p.second, opt.insert_sz, opt.min_cov_vertices)) {

				    					seq_final[p.first + pos_snp_km] = cf;

				    					s_valid_unitigs.insert(head);
				    				}
				    				else s_invalid_unitigs.insert(head);
				    			}
				    		}
				    	}
				    }

				    for (size_t i = 0; i < seq_final.length(); ++i){

				    	if (!isDNA(seq_final[i])) ud->add_ambiguity_char(i, seq_final[i]);
				    }

				    ud->runOptimizeAmbiguityChar();
				}

				++it_c;
			}
		};

		{
			const size_t slice = 1024;

	        bool stop = false;

	        size_t count_before = 0, count_after = 0;

	        vector<thread> workers; // need to keep track of threads so we can join them

	        mutex mutex_graph;

	        CompactedDBG<UnitigData>::iterator it_c = dbg.begin();
	        CompactedDBG<UnitigData>::iterator it_e = dbg.end();

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	            workers.emplace_back(

	                [&, t]{

						CompactedDBG<UnitigData>::iterator l_it_c, l_it_e;

	                    while (true) {

	                        {
	                            unique_lock<mutex> lock(mutex_graph);

	                            if (stop) return;

	                            l_it_c = it_c;

						        size_t i = 0;

						        count_before = count_after;

						        while ((i < slice) && (it_c != it_e)){

						        	++i;
						        	++it_c;
						        }

	                            if (it_c == it_e) stop = true;

	                            count_after += i;
								l_it_e = it_c;

								if (opt.verbose && ((count_before/1000000) != (count_after/1000000))) cout << "Ratatosk::detectSNPs(): Processed " << count_after << " vertices." << endl;
	                        }

	                        worker_function(l_it_c, l_it_e);
	                    }
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
	    }	
	}
}

/*unordered_set<uint64_t> preprocessGraphQueries(const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt) {

	const size_t k = dbg.getK();

	unordered_set<uint64_t> us_mini_mini;

    auto comp_pair = [](const pair<size_t, const_UnitigMap<UnitigData>>& p1, const pair<size_t, const_UnitigMap<UnitigData>>& p2) {

        return (p1.first < p2.first);
    };

	auto worker_function = [&](CompactedDBG<UnitigData>::const_iterator it_c, CompactedDBG<UnitigData>::const_iterator it_e) {

		unordered_set<uint64_t> us_mini_mini;

		while (it_c != it_e){

			const string unitig = it_c->referenceUnitigToString();

			vector<pair<size_t, const_UnitigMap<UnitigData>>> v_um = dbg.searchSequence(unitig, false, false, false, true, false);

			if (!v_um.empty()){

				size_t pos_v_um = 0;

				minHashIterator<RepHash> it_min_s(unitig.c_str(), unitig.length(), dbg.getK(), dbg.getG(), RepHash(), true), it_min_e;

				sort(v_um.begin(), v_um.end(), comp_pair);

				while (it_min_s != it_min_e) {

					const size_t unitig_pos = it_min_s.getKmerPosition();

					uint64_t min_h_min = it_min_s.getHash();

					while ((pos_v_um < v_um.size()) && (v_um[pos_v_um].first < unitig_pos)) ++pos_v_um;

					while ((pos_v_um < v_um.size()) && (v_um[pos_v_um].first == unitig_pos)) {

						const string kmer = v_um[pos_v_um].second.mappedSequenceToString();

						if (cstrMatch(kmer.c_str(), unitig.c_str() + unitig_pos) != k) {

							const minHashIterator<RepHash> l_it_min_s(kmer.c_str(), kmer.length(), dbg.getK(), dbg.getG(), RepHash(), true);

							min_h_min = min(min_h_min, l_it_min_s.getHash());
						}

						++pos_v_um;
					}

					us_mini_mini.insert(min_h_min);

					++it_min_s; // Increment k-mer position, does not jump to next minimizer in sequence
				}
			}

			++it_c;
		}

		return us_mini_mini;
	};

	{
		const size_t slice = 1024;

        bool stop = false;

        size_t count_before = 0, count_after = 0;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_graph, mutex_mini;

        CompactedDBG<UnitigData>::const_iterator it_c = dbg.begin();
        CompactedDBG<UnitigData>::const_iterator it_e = dbg.end();

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

					CompactedDBG<UnitigData>::const_iterator l_it_c, l_it_e;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_graph);

                            if (stop) return;

                            l_it_c = it_c;

					        size_t i = 0;

					        count_before = count_after;

					        while ((i < slice) && (it_c != it_e)){

					        	++i;
					        	++it_c;
					        }

                            if (it_c == it_e) stop = true;

                            count_after += i;
							l_it_e = it_c;

							if (opt.verbose && ((count_before/1000000) != (count_after/1000000))) cout << "Processed " << count_after << " vertices." << endl;
                        }

                        const unordered_set<uint64_t> l_us_mini_mini = worker_function(l_it_c, l_it_e);

                        {
                            unique_lock<mutex> lock(mutex_mini);

                            for (const uint64_t h : l_us_mini_mini) us_mini_mini.insert(h);
                        }
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
    }

	return us_mini_mini;
}*/

bool readGraphData(const string& input_filename, CompactedDBG<UnitigData>& dbg, const bool pid_only, const bool verbose){

    ifstream infile;
    istream in(0);

    Kmer head;

    if (verbose) cout << "Ratatosk::readGraphData(): Reading data associated to unitigs" << endl;

    infile.open(input_filename.c_str(), std::ofstream::in | std::ofstream::binary);
    in.rdbuf(infile.rdbuf());

    unordered_map<uint64_t, const PairID*, CustomHashUint64_t> m_sz;

    m_sz.reserve(dbg.size());

    while (in.good()){

    	head.read(in);

    	const UnitigMap<UnitigData> um = dbg.find(head, true);

    	if (!um.isEmpty) {

    		um.getData()->read(in, pid_only);

    		// See if we can switched our global owned PairID to a global shared PairID 
    		{
	    		SharedPairID& spid = um.getData()->getPairID();

	    		const pair<const PairID*, const PairID*> pids = spid.getPairIDs(); // Get pointers to <global, local> PairIDs

	    		if (pids.first != nullptr) { // IF there is a global PairID

	    			const vector<uint32_t> v_id = pids.first->toVector();
	    			const uint64_t h = XXH64(&v_id[0], v_id.size() * sizeof(uint32_t), 0); // Compute a hash from the global ids

	    			// Insert <hash, global PairID pointer>
	    			pair<unordered_map<uint64_t, const PairID*, CustomHashUint64_t>::iterator, bool> p = m_sz.insert(pair<uint64_t, const PairID*>(h, pids.first));

	    			if (!p.second) { // If the hash already exist, it means another SharedPairID might have the same Global PairID

	    				const PairID* g_pid = p.first->second;

	    				if (*g_pid == *(pids.first)) spid.setGlobalPairID(g_pid); // If the two global PairID have the same elements, 
	    			}
	    		}
    		}
    	}
    	else {

    		if (verbose) cerr << "Ratatosk::readGraphData(): Cannot associate data to unitig as the matching unitig is not found in the graph. Abort.";

    		infile.close();

    		return false;
    	}
    }

    return true;
}

void writeGraphData(const string& output_filename, const CompactedDBG<UnitigData>& dbg, const bool verbose){

    ofstream outfile;
    ostream out(0);

    if (verbose) cout << "Ratatosk::writeGraphData(): Writing data associated to unitigs" << endl;

    outfile.open(output_filename.c_str(), std::ofstream::out | std::ofstream::binary);
    out.rdbuf(outfile.rdbuf());

    for (const auto& um : dbg){

    	um.getUnitigHead().write(out);
    	um.getData()->write(out);
    }
}

void writeGraphPairID(CompactedDBG<UnitigData>& dbg, const string& out_fn, const Correct_Opt& opt) {

    ofstream outfile;
    ostream out(0);

    outfile.open(out_fn.c_str(), std::ofstream::out | std::ofstream::binary | std::ofstream::app);
    out.rdbuf(outfile.rdbuf());

    for (auto& um : dbg){

    	UnitigData* ud = um.getData();

    	if ((ud->getKmerCoverage(um) >= opt.max_cov_vertices) && (ud->getPairID().cardinality() >= opt.max_cov_vertices)) {

			um.getUnitigHead().write(out);

	    	ud->write(out, true);
	    	ud->clear(false, false, false, false, true, false, false);
    	}
    }
}

unordered_map<Kmer, vector<streampos>, KmerHash> mergeDiskPairIDs(const string& fn) {

	unordered_map<Kmer, vector<streampos>, KmerHash> umap_pos;

	Kmer head;

    ifstream infile;
    istream in(0);

    infile.open(fn.c_str(), std::ofstream::in | std::ofstream::binary);
    in.rdbuf(infile.rdbuf());

    while (in.good()){

  		SharedPairID l_spid;

   		if (in.good()) {

   			head.read(in);

   			umap_pos.insert({head, vector<streampos>()}).first->second.push_back(in.tellg());

    		if (in.good()) l_spid.read(in);
    	}
	}

    return umap_pos;	
}

void addCoverage(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, HapReads& hap_reads, const bool long_read_correct){

	if (dbg.isInvalid() || (dbg.length() == 0)) vector<Kmer>();

	vector<string> v_seq_in(opt.filename_seq_in);

	const string fn_tmp = opt.filename_long_out + ".tmp";

	const size_t k = dbg.getK();
	const size_t nb_pairs = countRecords(v_seq_in, !long_read_correct, opt.h_seed) + 1;

    const size_t thread_seq_buf_sz = opt.buffer_sz;
    const size_t thread_name_buf_sz = (thread_seq_buf_sz / (k + 1)) + 1;

    const char c_min_confidence_2nd_pass = getQual(opt.min_confidence_2nd_pass, opt.out_qual);

	const uint64_t off_limit_id = 0xffffffffffffffffULL;

	size_t nextID = 1;
	size_t nb_reads = 0;

	uint64_t last_id = 0;

	bool is_unallocated = true;

    double sampling_rate = 0.0;

    size_t len_read;
    size_t pos_read;

    string s, q;

    unordered_map<Kmer, vector<streampos>, KmerHash> umap_pids;

	vector<pair<uint64_t, uint64_t>> read_map_ids_sr(nb_pairs, pair<uint64_t, uint64_t>(off_limit_id, 0x0ULL));

    unordered_map<uint64_t, uint64_t, CustomHashUint64_t> name_hmap;

	FileParser fp(v_seq_in);

	LockGraph lck_g(opt.nb_threads * 1024);
	LockGraph lck_h(opt.nb_threads * 1024);

	name_hmap.reserve(nb_pairs);

	remove(fn_tmp.c_str());

    auto cmpPair = [](const pair<uint64_t, uint64_t>& a, const pair<uint64_t, uint64_t>& b) {

    	return ((a.second < b.second) || ((a.second == b.second) && (a.first < b.first)));
    };

    auto cmpKmer = [](const pair<Kmer, uint64_t>& a, const pair<Kmer, uint64_t>& b) {

    	return (a.second < b.second);
    };

    auto cmpKmerCovReversed = [](const pair<Kmer, uint64_t>& a, const pair<Kmer, uint64_t>& b) {

    	return (a.second > b.second);
    };

    auto cmpUnitigCovReversed = [](const UnitigMap<UnitigData>& a, const UnitigMap<UnitigData>& b) {

    	return (a.getData()->getKmerCoverage(a) > b.getData()->getKmerCoverage(b));
    };

    auto worker_function1 = [&](char* seq_buf, const size_t seq_buf_sz, const uint64_t* id_buf) {

    	size_t it_id_buf = 0;

        const char* str_end = seq_buf + seq_buf_sz;

        while (seq_buf < str_end) { // for each input

	        const int len = strlen(seq_buf);
	        const uint64_t id = id_buf[it_id_buf];

	        toUpperCase(seq_buf, len);

            KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;
            UnitigMap<UnitigData> um_canonical;

            unordered_set<Kmer, KmerHash> s_km;

            size_t len_centroid = 0;

            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

                const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
	            const UnitigMap<UnitigData> um = dbg.findUnitig(seq_buf, p_.second, len);

	            if (!um.isEmpty) { // Read maps to a Unitig

	            	const uint64_t h = um.getUnitigHead().hash();
	            	
	            	UnitigData* ud = um.getData();

	            	lck_g.lock_unitig(h);

					ud->increaseUnphasedCoverage(um.len); // No haplotype considered for now, we're interested purely in kmer coverage per unitig

	            	const size_t l_nb_reads_centroid = ud->getPairID().size();

	            	lck_g.unlock_unitig(h);

	            	if (um_canonical.isEmpty || ((l_nb_reads_centroid != 0) && (um.size > len_centroid)) || ((len_centroid == 0) && (um.size > um_canonical.size))) {

	            		um_canonical = um;

	            		if ((l_nb_reads_centroid != 0) && (um.size > len_centroid)) len_centroid = um.size;
	            	}

	            	if (!long_read_correct) s_km.insert(um.getUnitigHead().rep());

	            	it_kmer_h += um.len - 1;
	            }
	        }

	        if (!um_canonical.isEmpty){

            	const uint64_t h = um_canonical.getUnitigHead().hash();

            	lck_g.lock_unitig(h);

            	um_canonical.getData()->getPairID().add(id);

            	lck_g.unlock_unitig(h);
            }

            if (!long_read_correct) {

	        	lck_h.lock_unitig(id);

	        	for (const auto& km : s_km) read_map_ids_sr[id].second += km.hash();

	        	lck_h.unlock_unitig(id);
        	}

            seq_buf += len + 1;
            ++it_id_buf;
        }
    };

    auto worker_function2 = [&](char* seq_buf, const size_t seq_buf_sz, const uint64_t* id_buf) {

	    size_t it_id_buf = 0;

        const char* str_end = seq_buf + seq_buf_sz;

        while (seq_buf < str_end) { // for each input

	        const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_name_hmap = name_hmap.find(id_buf[it_id_buf]);

	    	const uint64_t read_name_id = (it_name_hmap == name_hmap.end()) ? 0 : it_name_hmap->second;

	        const size_t len = strlen(seq_buf);

            if ((read_name_id != 0) && (read_name_id != off_limit_id)) { // Read was selected by subsampling, is not duplicate nor invalid

	        	const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_read2hap = hap_reads.read2hap.find(id_buf[it_id_buf]);

		    	const uint64_t hap_name_id = (it_read2hap == hap_reads.read2hap.end()) ? off_limit_id : it_read2hap->second;

		    	toUpperCase(seq_buf, len);

	            KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;

	            {
		            lck_h.lock_unitig(hap_name_id + 1);

		            if (hap_name_id == off_limit_id) hap_reads.hap2unphasedReads.add(read_name_id - 1);
		            else hap_reads.hap2phasedReads[hap_name_id].add(read_name_id - 1);

		            lck_h.unlock_unitig(hap_name_id + 1);
	        	}

	            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

	                const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
		            const UnitigMap<UnitigData> um = dbg.findUnitig(seq_buf, p_.second, len);

		            if (!um.isEmpty){

		            	const uint64_t h = um.getUnitigHead().hash();

		            	UnitigData* ud = um.getData();

		            	lck_g.lock_unitig(h);

		            	ud->getPairID().add(read_name_id - 1);

		            	if (hap_name_id == off_limit_id) ud->increaseUnphasedCoverage(um.len);
		            	else {

		            		ud->increasePhasedCoverage(um.len);
		            		ud->get_hapID().add(hap_name_id);
		            	}

		            	lck_g.unlock_unitig(h);

		            	it_kmer_h += um.len - 1;
		            }
		        }
	    	}

            seq_buf += len + 1;
            ++it_id_buf;
        }
    };

    auto reading_function1 = [&](char* seq_buf, size_t& seq_buf_sz, uint64_t* name_hash_buf) {

        size_t file_id = 0;
        size_t it_name_hash_buf = 0;

        const size_t sz_seq_buf = thread_seq_buf_sz - k;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_seq_buf){

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

            	pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if ((len_read >= k) && (!long_read_correct || (len_read >= opt.min_len_2nd_pass))) {

                	if (pos_read == 0){ // New read

		            	const uint64_t h_name = XXH64(fp.getNameString(), strlen(fp.getNameString()), opt.h_seed);
			        	const pair<unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator, bool> p = name_hmap.insert({h_name, nextID});

			        	last_id = p.second ? nextID : p.first->second;
			        	nextID += static_cast<size_t>(p.second);

			        	if (long_read_correct) {

			        		const char* qual = fp.getQualityScoreString();

			        		for (size_t j = 0; j != len_read; ++j) {

			        			if (qual[j] < c_min_confidence_2nd_pass) s[j] = 'N';
			        		}
			        	}
                	}
                	
                	name_hash_buf[it_name_hash_buf++] = last_id;

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';

                        pos_read += sz_seq_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    auto reading_function2 = [&](char* seq_buf, size_t& seq_buf_sz, uint64_t* name_hash_buf) {

        size_t file_id = 0;
        size_t it_name_hash_buf = 0;

        const size_t sz_seq_buf = thread_seq_buf_sz - k;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_seq_buf){

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

            	pos_read &= static_cast<size_t>(new_reading) - 1;

            	if (opt.verbose && new_reading && ((++nb_reads % 1000000) == 0)) cout << "Ratatosk::addCoverage(): Processed " << nb_reads << " reads." << endl;

                len_read = s.length();
                s_str = s.c_str();

                if ((len_read >= k) && (!long_read_correct || (len_read >= opt.min_len_2nd_pass))) {

		        	if (long_read_correct && new_reading) {

		        		const char* qual = fp.getQualityScoreString();

		        		for (size_t j = 0; j != len_read; ++j) {

		        			if (qual[j] < c_min_confidence_2nd_pass) s[j] = 'N';
		        		}
		        	}

                	name_hash_buf[it_name_hash_buf++] = XXH64(fp.getNameString(), strlen(fp.getNameString()), opt.h_seed);

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';

                        pos_read += sz_seq_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    auto getDiskPairID = [&](istream& in, const const_UnitigMap<UnitigData>& um) {

    	SharedPairID spid = um.getData()->getPairID();

		const unordered_map<Kmer, vector<streampos>, KmerHash>::const_iterator it_umap_um = umap_pids.find(um.getUnitigHead());

		if (it_umap_um != umap_pids.end()) {

			for (const auto fpos : it_umap_um->second) {

				SharedPairID spid_tmp;

				in.seekg(fpos);
				spid_tmp.read(in);

				spid |= spid_tmp.toPairID();
			}
		}

		return spid;
    };

    auto getDiskPairIDs = [&](istream& in, const const_UnitigMap<UnitigData>& um) {

		const unordered_map<Kmer, vector<streampos>, KmerHash>::const_iterator it_umap_um = umap_pids.find(um.getUnitigHead());

    	SharedPairID spid = um.getData()->getPairID();

		vector<SharedPairID> v_spids;

		if (it_umap_um != umap_pids.end()) {

			v_spids.reserve(it_umap_um->second.size() + 1);
			v_spids.push_back(move(spid));

			for (const auto fpos : it_umap_um->second) {

				in.seekg(fpos);
				spid.read(in);

				v_spids.push_back(move(spid));
			}
		}
		else v_spids.push_back(move(spid));

		return v_spids;
    };

	auto postProcessUnitigs = [&](istream& in, CompactedDBG<UnitigData>::iterator it_s, CompactedDBG<UnitigData>::iterator it_e) {

		while (it_s != it_e) {

			UnitigData* ud = it_s->getData();

			const SharedPairID pids = getDiskPairID(in, *it_s);

			if (ud->getPhasedKmerCoverage(*it_s) < opt.min_cov_vertices) ud->get_hapID().clear(); // Not enough k-mer coverage for the phased reads -> delete
			else ud->get_hapID().runOptimize();

			ud->setBranching((it_s->getPredecessors().cardinality() > 1) || (it_s->getSuccessors().cardinality() > 1));

			for (const auto& um_succ : it_s->getSuccessors()) {

				const SharedPairID pids_succ = getDiskPairID(in, um_succ);

				if (getNumberSharedPairID(pids, pids_succ, opt.min_cov_vertices) >= opt.min_cov_vertices) ud->setSharedPids(it_s->strand, um_succ.getMappedHead().toString()[k-1]);
			}

			{
				UnitigMap<UnitigData> um_rev(*it_s);

				um_rev.strand = !it_s->strand;

				for (const auto& um_succ : um_rev.getSuccessors()) {

					const SharedPairID pids_succ = getDiskPairID(in, um_succ);

					if (getNumberSharedPairID(pids, pids_succ, opt.min_cov_vertices) >= opt.min_cov_vertices) ud->setSharedPids(um_rev.strand, um_succ.getMappedHead().toString()[k-1]);
				}
			}

			++it_s;
		}
	};

    {
    	if (opt.verbose) cout << "Ratatosk::addCoverage(): Anchoring reads on graph." << endl;

        char** buffer_seq = new char*[opt.nb_threads];
        size_t* buffer_seq_sz = new size_t[opt.nb_threads];
        uint64_t** buffer_name = new uint64_t*[opt.nb_threads];

	    len_read = 0;
	    pos_read = 0;
	    nextID = 1;

	    if (opt.nb_threads == 1){

        	buffer_seq[0] = new char[thread_seq_buf_sz];
        	buffer_name[0] = new uint64_t[thread_name_buf_sz];

            while (!reading_function1(buffer_seq[0], buffer_seq_sz[0], buffer_name[0])) worker_function1(buffer_seq[0], buffer_seq_sz[0], buffer_name[0]);
	    }
	    else {

        	bool eof_fastx = false;

	        vector<thread> workers; // need to keep track of threads so we can join them

	        mutex mutex_file;

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	        	buffer_seq[t] = new char[thread_seq_buf_sz];
	        	buffer_name[t] = new uint64_t[thread_name_buf_sz];

	            workers.emplace_back(

	                [&, t]{

	                    while (true) {

	                        {
	                            unique_lock<mutex> lock(mutex_file);

	                            if (eof_fastx) return;

	                            eof_fastx = reading_function1(buffer_seq[t], buffer_seq_sz[t], buffer_name[t]);
	                        }

	                        worker_function1(buffer_seq[t], buffer_seq_sz[t], buffer_name[t]);
	                    }
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
	    }

        for (size_t t = 0; t < opt.nb_threads; ++t){

        	delete[] buffer_seq[t];
        	delete[] buffer_name[t];
        }

        delete[] buffer_seq;
        delete[] buffer_name;
        delete[] buffer_seq_sz;

	    fp.close();

	    if (opt.verbose) cout << "Ratatosk::addCoverage(): Detecting and removing duplicated reads." << endl;

		nextID = 1;

	    std::random_device rd; //Seed
	    std::default_random_engine generator(rd()); //Random number generator
	    std::uniform_real_distribution<> distribution(0.0, 1.0); //Distribution on which to apply the generator

	    for (const auto& um : dbg) {

			const SharedPairID& pid = um.getData()->getPairID();

			vector<pair<uint64_t, uint64_t>> vp;

			vp.reserve(pid.size());

			for (const auto id : pid) vp.push_back(pair<uint64_t, uint64_t>(id, read_map_ids_sr[id].second));

			sort(vp.begin(), vp.end(), cmpPair);

			size_t prev_id = nextID;

			for (size_t i = 0; i < vp.size(); ++i){

				pair<uint64_t, uint64_t>& p = read_map_ids_sr[vp[i].first];

				if (long_read_correct) {

					if (p.first == off_limit_id) p.first = ((opt.sampling_rate == 1.0) || (distribution(generator) <= opt.sampling_rate)) ? nextID++ : 0;
				}
				else if (p.first == off_limit_id){

					if ((i == 0) || (vp[i].second != vp[i-1].second)) {

						prev_id = ((opt.sampling_rate == 1.0) || (distribution(generator) <= opt.sampling_rate)) ? nextID++ : 0;
					}
					
					p.first = prev_id;
				}
				else prev_id = p.first;
			}
		}

		for (auto& h : name_hmap) h.second = read_map_ids_sr[h.second].first;

		read_map_ids_sr.clear();

		for (auto& um : dbg) um.getData()->clear(); // Clear all unitig data
	}

    {
    	if (opt.verbose) cout << "Ratatosk::addCoverage(): Coloring graph with reads." << endl;

    	const size_t limit_sz_read = 0x00000000ffffffffULL;

        char** buffer_seq = new char*[opt.nb_threads];
        size_t* buffer_seq_sz = new size_t[opt.nb_threads];
        uint64_t** buffer_name = new uint64_t*[opt.nb_threads];

    	size_t curr_sz_read = 0;

	    vector<string> v_fn_out_part; // vector of temporary file(name)s

	    len_read = 0;
	    pos_read = 0;

	    fp = FileParser(v_seq_in);

	    if (opt.nb_threads == 1){

        	buffer_seq[0] = new char[thread_seq_buf_sz];
        	buffer_name[0] = new uint64_t[thread_name_buf_sz];

            while (!reading_function2(buffer_seq[0], buffer_seq_sz[0], buffer_name[0])) {

            	worker_function2(buffer_seq[0], buffer_seq_sz[0], buffer_name[0]);

            	curr_sz_read += buffer_seq_sz[0];

            	if (curr_sz_read >= limit_sz_read) {

		        	writeGraphPairID(dbg, fn_tmp, opt);

		        	curr_sz_read = 0;
            	}
            }
	    }
	    else {

        	bool eof_fastx = false;

	        vector<thread> workers; // need to keep track of threads so we can join them

	        mutex mutex_file;

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	        	buffer_seq[t] = new char[thread_seq_buf_sz];
	        	buffer_name[t] = new uint64_t[thread_name_buf_sz];
	        }

	        while (!eof_fastx) {

		        for (size_t t = 0; t < opt.nb_threads; ++t){

		            workers.emplace_back(

		                [&, t]{

		                    while (true) {

		                        {
		                            unique_lock<mutex> lock(mutex_file);

		                            if (eof_fastx || (curr_sz_read >= limit_sz_read)) return;

		                            eof_fastx = reading_function2(buffer_seq[t], buffer_seq_sz[t], buffer_name[t]);
					            	curr_sz_read += buffer_seq_sz[t];
		                        }

		                        worker_function2(buffer_seq[t], buffer_seq_sz[t], buffer_name[t]);
		                    }
		                }
		            );
		        }

		        for (auto& t : workers) t.join();

		        if (curr_sz_read >= limit_sz_read) {

		        	writeGraphPairID(dbg, fn_tmp, opt);

		        	curr_sz_read = 0;

		        	workers.clear();
		        }
	    	}
    	}

        for (size_t t = 0; t < opt.nb_threads; ++t){

        	delete[] buffer_seq[t];
        	delete[] buffer_name[t];
        }

        delete[] buffer_seq;
        delete[] buffer_name;
        delete[] buffer_seq_sz;

        name_hmap.clear();
        fp.close();

    	if (check_files(fn_tmp, false, false)) umap_pids = mergeDiskPairIDs(fn_tmp); // Reassemble data from disk
	}

	{
		if (opt.verbose) cout << "Ratatosk::addCoverage(): Collecting branching information." << endl;

	    const size_t slice = 1024;

	    CompactedDBG<UnitigData>::iterator it_s = dbg.begin();
	    CompactedDBG<UnitigData>::iterator it_e = dbg.end();

	    if ((opt.nb_threads == 1) || (dbg.size() < slice)){

	    	ifstream infile;
	    	istream in(0);

	    	if (check_files(fn_tmp, false, false)) {

		    	infile.open(fn_tmp.c_str(), std::ofstream::in | std::ofstream::binary);
		    	in.rdbuf(infile.rdbuf());
	    	}

	    	postProcessUnitigs(in, it_s, it_e);
	    }
	    else {

	        vector<thread> workers;

	        mutex mutex_graph;

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	            workers.emplace_back(

	                [&]{

	                	CompactedDBG<UnitigData>::iterator it_s_l;
	                	CompactedDBG<UnitigData>::iterator it_e_l;

				    	ifstream infile;
				    	istream in(0);

				    	if (check_files(fn_tmp, false, false)) {

					    	infile.open(fn_tmp.c_str(), std::ofstream::in | std::ofstream::binary);
					    	in.rdbuf(infile.rdbuf());
				    	}

	                	while (true){

	                		{
	                			unique_lock<mutex> lock(mutex_graph);

	                			if (it_s == it_e) return;
	                			else {

	                				it_s_l = it_s;
	                				it_e_l = it_s;

	                				for (size_t i = 0; (i < slice) && (it_e_l != it_e); ++i) ++it_e_l;

	                				it_s = it_e_l;
	                			}
	                		}

	                		postProcessUnitigs(in, it_s_l, it_e_l);
	                	}
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
	    }
	}

	const size_t estimated_hap_cov = estimateHaplotypeCoverage(dbg);

	if (estimated_hap_cov >= 10) {

		if (opt.verbose) cout << "Ratatosk::addCoverage(): Subsampling reads." << endl;

		fstream in(0);

		vector<pair<Kmer, size_t>> v_km;

		const double step = 0.05;

		const size_t sz_bitset = ((nextID + 63) / 64) + 1;

		double sampling_rate = 5.0 / static_cast<double>(estimated_hap_cov);

	    std::random_device rd; // Seed
	    std::default_random_engine generator(rd()); // Random number generator
	    std::uniform_real_distribution<> distribution(0.0, 1.0); // Distribution on which to apply the generator

		uint64_t* all_ids = new uint64_t[sz_bitset]();
		uint64_t* all_subsampled_ids = new uint64_t[sz_bitset]();
		uint64_t* sampled_ids = new uint64_t[sz_bitset]();

		if (check_files(fn_tmp, false, false)) in.open(fn_tmp.c_str(), ios::in | ios::out | ios::binary | ios::app);

		v_km.reserve(dbg.size());

		for (const auto& um : dbg) v_km.push_back(pair<Kmer, size_t>(um.getUnitigHead(), um.getData()->getKmerCoverage(um)));

		sort(v_km.begin(), v_km.end(), cmpKmerCovReversed); // Sort unitigs by coverage

		size_t prev_max_km_cov = 0;
		size_t prev_min_km_cov = 0;

		for (const auto& um : dbg) {

			const UnitigData* ud = um.getData();

			if (!ud->isBranching()) {

				const SharedPairID spid = getDiskPairID(in, um);
				const PairID pid = subsample(spid, opt.min_cov_vertices);

				for (const uint32_t id : pid) all_subsampled_ids[id >> 6] |= (0x1ULL << (id & 0x3fULL));
			}
		}

		for (double i = 1.0; i > 0.0; i -= step) {

			const size_t pos_min = v_km.size() * i - static_cast<size_t>((v_km.size() * i) == v_km.size());
			const size_t pos_max = v_km.size() * (i - step);

			const size_t min_km_cov = v_km[pos_min].second;
			const size_t max_km_cov = v_km[pos_max].second;

			if ((max_km_cov != min_km_cov) && (max_km_cov != prev_max_km_cov)) {

				memset(sampled_ids, 0, sz_bitset * sizeof(uint64_t));

				if (opt.nb_threads == 1) {
				
					for (auto& um : dbg) {

						const size_t um_km_cov = um.getData()->getKmerCoverage(um);

						if ((um_km_cov >= min_km_cov) && (um_km_cov < max_km_cov)) {

							const vector<SharedPairID> v_spids = getDiskPairIDs(in, um);

							for (const auto& spid : v_spids) {

								for (const uint32_t id : spid) sampled_ids[id >> 6] |= (0x1ULL << (id & 0x3fULL));
							}
						}
					}
				}
				else {

					const size_t slice = 1024;

			        vector<thread> workers;

			        mutex mtx_g, mtx_f;

			        CompactedDBG<UnitigData>::iterator it_s = dbg.begin();
			        CompactedDBG<UnitigData>::iterator it_e = dbg.end();

			        uint64_t** l_sampled_ids = new uint64_t*[opt.nb_threads];

			        for (size_t t = 0; t < opt.nb_threads; ++t){

			        	l_sampled_ids[t] = new uint64_t[sz_bitset]();

			            workers.emplace_back(

			                [&, t]{

			                	CompactedDBG<UnitigData>::iterator it_s_l;
			                	CompactedDBG<UnitigData>::iterator it_e_l;

			                	while (true){

			                		{
			                			unique_lock<mutex> lock(mtx_g);

			                			if (it_s == it_e) return;
			                			else {

			                				it_s_l = it_s;
			                				it_e_l = it_s;

			                				for (size_t j = 0; (j < slice) && (it_e_l != it_e); ++j) ++it_e_l;

			                				it_s = it_e_l;
			                			}
			                		}

			                		while (it_s_l != it_e_l) {

										const size_t um_km_cov = it_s_l->getData()->getKmerCoverage(*it_s_l);

										if ((um_km_cov >= min_km_cov) && (um_km_cov < max_km_cov)) {

											vector<SharedPairID> v_spids;

					                		{
					                			unique_lock<mutex> lock(mtx_f);

												v_spids = getDiskPairIDs(in, *it_s_l);
											}

											for (const auto& spid : v_spids) {

												for (const uint32_t id : spid) l_sampled_ids[t][id >> 6] |= (0x1ULL << (id & 0x3fULL));
											}
										}

										++it_s_l;
			                		}
			                	}
			                }
			            );
			        }

			        for (auto& t : workers) t.join();

			        for (size_t t = 0; t < opt.nb_threads; ++t){

			        	for (size_t j = 0; j < sz_bitset; ++j) sampled_ids[j] |= l_sampled_ids[t][j];

			        	delete[] l_sampled_ids[t];
			        }

			        delete[] l_sampled_ids;
				}

				if (min_km_cov >= 5) { // No subsampling for very low coverage vertices

					// Compute sampling rate
					if (opt.verbose) cout << "Ratatosk::addCoverage(): Subsampling reads (" << sampling_rate << ") coloring vertices with k-mer coverage [" << min_km_cov << ", " << max_km_cov << "[" << endl;

					for (size_t i = 0; i < sz_bitset; ++i) {

						sampled_ids[i] &= ~all_ids[i]; // Do not select ids that were already selected previously

						if (sampled_ids[i] != 0) {

							for (size_t j = 0; j < 64; ++j) {

								if (static_cast<bool>(sampled_ids[i] & (0x1ULL << j)) && (distribution(generator) <= sampling_rate)) all_subsampled_ids[i] |= (0x1ULL << j);
							}
						}
					}
				}
				else {

					for (size_t i = 0; i < sz_bitset; ++i) all_subsampled_ids[i] |= sampled_ids[i];
				}

				for (size_t i = 0; i < sz_bitset; ++i) all_ids[i] |= sampled_ids[i];
			}

			prev_max_km_cov = max_km_cov;
			prev_min_km_cov = min_km_cov;
		}

		if (opt.verbose) cout << "Ratatosk::addCoverage(): Coloring graph with subsampled reads." << endl;

		if (opt.nb_threads == 1) {

			const bool hasDiskPairID = !umap_pids.empty();

			for (auto& um : dbg) {

				const Kmer head = um.getUnitigHead();

				const SharedPairID disk_spids = getDiskPairID(in, um);

				SharedPairID& spids = um.getData()->getPairID();

				spids.clear();

				for (const uint32_t id : disk_spids) {

					if (static_cast<bool>(all_subsampled_ids[id >> 6] & (0x1ULL << (id & 0x3fULL)))) spids.add(id);
				}

				spids.runOptimize();

				if (hasDiskPairID) {

					const unordered_map<Kmer, vector<streampos>, KmerHash>::iterator it_umap_um = umap_pids.find(head);

					if (it_umap_um != umap_pids.end()) {

						in.seekg(0, ios::end); // Place file cursor at the end of the file

						it_umap_um->second.clear(); // Clear all disk PairID associated with that unitig
						it_umap_um->second.push_back(in.tellg()); // Insert position of next disk PairID associated with that unitig

						head.write(in); // Write to disk position of unitig head k-mer

	    				um.getData()->write(in, true); // Write to disk the new PairID of that unitig
	    				um.getData()->clear(false, false, false, false, true, false, false); //
					}
				}
			}
		}
		else {

			const size_t slice = (dbg.size() / opt.nb_threads) + 1;

			const bool hasDiskPairID = !umap_pids.empty();

	        vector<thread> workers;

	        mutex mtx_g, mtx_f;

	        CompactedDBG<UnitigData>::iterator it_s = dbg.begin();
	        CompactedDBG<UnitigData>::iterator it_e = dbg.end();

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	            workers.emplace_back(

	                [&]{

	                	CompactedDBG<UnitigData>::iterator it_s_l;
	                	CompactedDBG<UnitigData>::iterator it_e_l;

	                	while (true){

	                		{
	                			unique_lock<mutex> lock(mtx_g);

	                			if (it_s == it_e) return;
	                			else {

	                				it_s_l = it_s;
	                				it_e_l = it_s;

	                				for (size_t i = 0; (i < slice) && (it_e_l != it_e); ++i) ++it_e_l;

	                				it_s = it_e_l;
	                			}
	                		}

	                		while (it_s_l != it_e_l) {

	                			if (it_s_l->getData()->isBranching()) {

									const Kmer head = it_s_l->getUnitigHead();

									SharedPairID& spids = it_s_l->getData()->getPairID();

									SharedPairID spid_tmp;

			                		if (hasDiskPairID) {

			                			unique_lock<mutex> lock(mtx_f);

			                			spid_tmp = getDiskPairID(in, *it_s_l);
			                		}
			                		else spid_tmp = spids;

									spids.clear();

									for (const uint32_t id : spid_tmp) {

										if (static_cast<bool>(all_subsampled_ids[id >> 6] & (0x1ULL << (id & 0x3fULL)))) spids.add(id);
									}

									spids.runOptimize();

			                		if (hasDiskPairID) {

										const unordered_map<Kmer, vector<streampos>, KmerHash>::iterator it_umap_um = umap_pids.find(head);

										if (it_umap_um != umap_pids.end()) {

			                				unique_lock<mutex> lock(mtx_f);

											in.seekg(0, ios::end); // Place file cursor at the end of the file

											it_umap_um->second.clear(); // Clear all disk PairID associated with that unitig
											it_umap_um->second.push_back(in.tellg()); // Insert position of next disk PairID associated with that unitig

											head.write(in); // Write to disk position of unitig head k-mer

						    				it_s_l->getData()->write(in, true); // Write to disk the new PairID of that unitig
						    				it_s_l->getData()->clear(false, false, false, false, true, false, false); // Clear data
										}
									}
								}

								++it_s_l;
	                		}
	                	}
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
		}

		delete[] all_ids;
		delete[] all_subsampled_ids;
		delete[] sampled_ids;
	}

	{
		if (opt.verbose) cout << "Ratatosk::addCoverage(): Compacting colors." << endl;

    	ifstream infile;
    	istream in(0);

   		vector<pair<Kmer, size_t>> v_km;

    	if (check_files(fn_tmp, false, false)) {

	    	infile.open(fn_tmp.c_str(), std::ofstream::in | std::ofstream::binary);
	    	in.rdbuf(infile.rdbuf());
    	}

		v_km.reserve(dbg.size());

		for (const auto& um : dbg) {

			// Only use branching vertices
			if (um.getData()->isBranching() && (um.getData()->getKmerCoverage(um) >= 3 * estimated_hap_cov)) {

				v_km.push_back(pair<Kmer, size_t>(um.getUnitigHead(), um.getData()->getKmerCoverage(um)));
			}
		}

		sort(v_km.begin(), v_km.end(), cmpKmerCovReversed); // Sort unitigs by coverage

		for (const auto& p_km : v_km) {

			const UnitigMap<UnitigData> um = dbg.find(p_km.first, true);

			if (!um.isEmpty) {

				UnitigData* um_ud = um.getData();

				if (!um_ud->getVisitStatus()) {

					const SharedPairID spids = getDiskPairID(in, um);

					PairID inter = spids.toPairID();

					size_t max_card_inter = spids.cardinality() * opt.min_color_sharing;

					unordered_set<Kmer, KmerHash> s_km_visited, s_km_valid; // To remember which nodes have been visited locally
					queue<UnitigMap<UnitigData>> q;

					s_km_visited.insert(p_km.first.rep());
					q.push(um);

					while (!q.empty()) {

						const UnitigMap<UnitigData> l_um = q.front();

						vector<UnitigMap<UnitigData>> v_um_neigh;

						v_um_neigh.reserve(8); // Only 8 possible neighbors (4 predecessors + 4 successors)
						q.pop();

						for (const auto& l_um_neigh : l_um.getSuccessors()) {

							if (s_km_visited.insert(l_um_neigh.getUnitigHead().rep()).second && !l_um_neigh.getData()->getVisitStatus()) v_um_neigh.push_back(l_um_neigh);
						}

						for (const auto& l_um_neigh : l_um.getPredecessors()) {

							if (s_km_visited.insert(l_um_neigh.getUnitigHead().rep()).second && !l_um_neigh.getData()->getVisitStatus()) v_um_neigh.push_back(l_um_neigh);
						}

						sort(v_um_neigh.begin(), v_um_neigh.end(), cmpUnitigCovReversed); // Sort unitigs by coverage

						for (const auto& l_um_neigh : v_um_neigh) {

							const SharedPairID l_spid_neigh = getDiskPairID(in, l_um_neigh);

							PairID l_inter = inter & l_spid_neigh.toPairID();

							// New intersection shares at least (opt.min_color_sharing * 100)% of its colors with previously "accepted" unitigs
							if ((l_inter.cardinality() >= (l_spid_neigh.cardinality() * opt.min_color_sharing)) && (l_inter.cardinality() >= max_card_inter)) {

								inter = move(l_inter);
								max_card_inter = max(max_card_inter, static_cast<size_t>(l_spid_neigh.cardinality() * opt.min_color_sharing));

								s_km_valid.insert(l_um_neigh.getUnitigHead().rep());
								q.push(l_um_neigh);
							}
						}
					}

					if (!s_km_valid.empty()) {

						SharedPairID new_spid;

						inter.runOptimize();

						new_spid.setGlobalPairID(inter);
						new_spid.setLocalPairID(spids.toPairID());

						new_spid.runOptimize();

						um_ud->getPairID() = move(new_spid);
						um_ud->setVisitStatus(true);

						const pair<const PairID*, const PairID*> p = um.getData()->getPairID().getPairIDs();

						for (const auto& km : s_km_valid) {

							const UnitigMap<UnitigData> um_neigh = dbg.find(km, true);

							if (!um_neigh.isEmpty) {

								const SharedPairID l_spid = getDiskPairID(in, um_neigh);

								SharedPairID l_new_spid;

								l_new_spid.setGlobalPairID(p.first);
								l_new_spid.setLocalPairID(l_spid.toPairID());

								l_new_spid.runOptimize();

								um_neigh.getData()->getPairID() = move(l_new_spid);
								um_neigh.getData()->setVisitStatus(true);
							}
						}
					}
				}
			}
		}

		if (opt.nb_threads == 1) {

			for (const auto& um : dbg){

				UnitigData* ud = um.getData();

				if (!ud->getVisitStatus()) {

					SharedPairID spids = getDiskPairID(in, um);

					ud->getPairID() = move(spids);
					ud->getPairID().runOptimize();
				}
				else ud->setVisitStatus(false);
			}
		}
		else {

			const size_t slice = (dbg.size() / opt.nb_threads) + 1;

	        vector<thread> workers;

	        mutex mtx_g, mtx_f;

	        CompactedDBG<UnitigData>::iterator it_s = dbg.begin();
	        CompactedDBG<UnitigData>::iterator it_e = dbg.end();

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	            workers.emplace_back(

	                [&]{

	                	CompactedDBG<UnitigData>::iterator it_s_l;
	                	CompactedDBG<UnitigData>::iterator it_e_l;

	                	while (true){

	                		{
	                			unique_lock<mutex> lock(mtx_g);

	                			if (it_s == it_e) return;
	                			else {

	                				it_s_l = it_s;
	                				it_e_l = it_s;

	                				for (size_t i = 0; (i < slice) && (it_e_l != it_e); ++i) ++it_e_l;

	                				it_s = it_e_l;
	                			}
	                		}

	                		while (it_s_l != it_e_l) {

								UnitigData* ud = it_s_l->getData();

								if (!ud->getVisitStatus()) {

									SharedPairID spids;

									{
										unique_lock<mutex> lock(mtx_f);

										spids = getDiskPairID(in, *it_s_l);
									}

									ud->getPairID() = move(spids);
									ud->getPairID().runOptimize();
								}
								else ud->setVisitStatus(false);

								++it_s_l;
	                		}
	                	}
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
		}
	}

	if (!umap_pids.empty() && (remove(fn_tmp.c_str()) != 0)) cerr << "Ratatosk::addCoverage(): Could not remove tmp file " << fn_tmp << endl;
}

void addPhasing(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, HapReads& hap_r, const bool long_read_phasing, const bool mapHapReads) {

	const size_t k = dbg.getK();
	const uint64_t off_limit_id = 0xffffffffffffffffULL;

    string s;

	LockGraph lck_g(opt.nb_threads * 1024);

    size_t len_read;
    size_t pos_read;

    size_t nb_reads = 0;

    const size_t thread_seq_buf_sz = opt.buffer_sz;
    const size_t thread_name_buf_sz = (thread_seq_buf_sz / (k + 1)) + 1;

    const vector<string>& v_phase_filenames = long_read_phasing ? opt.filenames_long_phase : opt.filenames_short_phase;

    vector<string> v_read_filenames = long_read_phasing ? opt.filenames_long_in : opt.filename_seq_in;

    if (long_read_phasing) v_read_filenames.insert(v_read_filenames.end(), opt.filenames_helper_long_in.begin(), opt.filenames_helper_long_in.end());

    auto reading_function1 = [&](FileParser& fp, char* seq_buf, size_t& seq_buf_sz, uint64_t* hap_hash_buf) {

        size_t file_id = 0;
        size_t it_name_hash_buf = 0;

        const size_t sz_seq_buf = thread_seq_buf_sz - k;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_seq_buf){

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

            	pos_read &= static_cast<size_t>(new_reading) - 1;

            	if (opt.verbose && new_reading && ((++nb_reads % 1000000) == 0)) cout << "Ratatosk::addPhasing(): Processed " << nb_reads << " reads." << endl;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k){

                	const uint64_t h_name = XXH64(fp.getNameString(), strlen(fp.getNameString()), opt.h_seed);
                	const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_read2hap = hap_r.read2hap.find(h_name);

                	if ((it_read2hap != hap_r.read2hap.end()) && (it_read2hap->second != off_limit_id)){

                		hap_hash_buf[it_name_hash_buf++] = it_read2hap->second;

	                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

	                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

	                        seq_buf[thread_seq_buf_sz - 1] = '\0';

	                        pos_read += sz_seq_buf - seq_buf_sz;
	                        seq_buf_sz = thread_seq_buf_sz;

	                        break;
	                    }
	                    else {

	                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

	                        seq_buf_sz += (len_read - pos_read) + 1;
	                        pos_read = len_read;
	                    }
                	}
                	else pos_read = len_read;
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    auto worker_function1 = [&](char* seq_buf, const size_t seq_buf_sz, const uint64_t* hap_buf) {

    	uint64_t it_min_h;

	    RepHash rep;

	    size_t it_hap_buf = 0;

        const char* str_end = seq_buf + seq_buf_sz;

        while (seq_buf < str_end) { // for each input

	        const int len = strlen(seq_buf);

        	toUpperCase(seq_buf, len);

            KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;

            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

                const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
	            const UnitigMap<UnitigData> um = dbg.findUnitig(seq_buf, p_.second, len);

	            if (!um.isEmpty) { // Read maps to a Unitig

	            	const uint64_t h = um.getUnitigHead().hash();

	            	UnitigData* ud = um.getData();

	            	lck_g.lock_unitig(h);

	            	ud->get_hapID().add(hap_buf[it_hap_buf]);
		            //ud->increasePhasedCoverage(um.len);

	            	lck_g.unlock_unitig(h);

	            	it_kmer_h += um.len - 1;
	            }
	        }

            seq_buf += len + 1;
            ++it_hap_buf;
        }
    };

	auto loadPhasing = [&](const vector<string>& v_filenames_phase) {

		const size_t sz_buffer = 16384;

    	char buffer[sz_buffer];

    	int buffer_occupancy = 0;

    	uint64_t hap2read_sz = 0;

    	if (opt.verbose) cout << "Ratatosk::addPhasing(): Associating " << (long_read_phasing ? "long" : "short") << " reads to their haploblocks" << endl;

		for (const auto& filename : v_filenames_phase){

            gzFile gfp = gzopen(filename.c_str(), "r");

            if (gfp != Z_NULL) {

            	int l_read = gzread(gfp, buffer, sz_buffer);

	            while (l_read != 0){

	            	const char* l_buffer = buffer;
	            	const char* l_next_buffer = nullptr;

	            	size_t l_buffer_occupancy = buffer_occupancy + l_read;

	            	l_next_buffer = static_cast<const char*>(memchr(l_buffer, '\n', l_buffer_occupancy));

	            	while (l_next_buffer != nullptr) {

		            	if (l_buffer[0] != '#'){ // This is a comment or a header

						    const char* haplotype_cstr = strchr(l_buffer, '\t');
						    const char* haploblock_cstr = (haplotype_cstr == nullptr) ? nullptr : strchr(haplotype_cstr + 1, '\t');

						    if (haploblock_cstr != nullptr) {

						    	const string haplotype_str(haplotype_cstr + 1, haploblock_cstr - (haplotype_cstr + 1));
						    	const string haploblock_str(haploblock_cstr + 1, l_next_buffer - (haploblock_cstr + 1));

						    	if ((haploblock_str != "none") && (haplotype_str != "none")) {

							    	const pair<unordered_map<string, uint64_t, CustomHashString>::const_iterator, bool> p_block = hap_r.hapBlock2id.insert({haploblock_str, hap_r.block_id});
							    	const pair<unordered_map<string, uint64_t, CustomHashString>::const_iterator, bool> p_type = hap_r.hapType2id.insert({haplotype_str, hap_r.type_id});

							    	const uint64_t hap_block_type = (p_block.first->second << 1) + p_type.first->second;
							    	const uint64_t h_record = XXH64(l_buffer, haplotype_cstr - l_buffer, opt.h_seed);

							    	hap_r.block_id += static_cast<size_t>(p_block.second);
							    	hap_r.type_id += static_cast<size_t>(p_type.second);

							    	hap2read_sz = max(hap2read_sz, hap_block_type);

					        		pair<unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::iterator, bool> p_read2hap = hap_r.read2hap.insert({h_record, hap_block_type});

					        		// Insertion did not take place:
					        		// - Long read phasing: ID collision -> Set read ID to unphased
					        		// - Short read phasing: Collision is due to mates from same pair. If Haploblock is incompatible, set read ID to unphased
					        		if (!p_read2hap.second && (long_read_phasing || (!long_read_phasing && (p_read2hap.first->second != hap_block_type)))) p_read2hap.first->second = off_limit_id;
				        		}
						    }
							else cerr << "Ratatosk::addPhasing(): Malformed phasing record. Record is discarded." << endl;
						}

						l_buffer_occupancy -= (l_next_buffer - l_buffer) + 1;
						l_buffer = l_next_buffer + 1;
						l_next_buffer = static_cast<const char*>(memchr(l_buffer, '\n', l_buffer_occupancy));
					}

					if (l_buffer_occupancy != 0) memmove(buffer, l_buffer, l_buffer_occupancy);
					
					buffer_occupancy = l_buffer_occupancy;
					l_read = gzread(gfp, (buffer + buffer_occupancy), (sz_buffer - buffer_occupancy));
	            }

	            gzclose(gfp);

	            hap_r.hap2phasedReads = vector<PairID>(hap2read_sz + 1);
        	}
        	else cerr << "Ratatosk::addPhasing(): Could not open file " << filename << " for phasing. File is discarded." << endl;
    	}

    	if (opt.verbose) cout << "Ratatosk::addPhasing(): " << hap_r.read2hap.size() << " phased reads added." << endl;
    	if (opt.verbose) cout << "Ratatosk::addPhasing(): " << hap_r.hapBlock2id.size() << " haploblocks with " << hap_r.hapType2id.size() << " haplotypes compacted." << endl;
	};

    auto runHapMapping = [&](const vector<string>& filenames) {

    	if (opt.verbose) cout << "Ratatosk::addPhasing(): Associating haploblocks to unitigs." << endl;

        char** buffer_seq = new char*[opt.nb_threads];
        size_t* buffer_seq_sz = new size_t[opt.nb_threads];
        uint64_t** buffer_name = new uint64_t*[opt.nb_threads];

        FileParser fp(filenames);

	    len_read = 0;
	    pos_read = 0;

	    if (opt.nb_threads == 1){

        	buffer_seq[0] = new char[thread_seq_buf_sz];
        	buffer_name[0] = new uint64_t[thread_name_buf_sz];

            while (!reading_function1(fp, buffer_seq[0], buffer_seq_sz[0], buffer_name[0])) worker_function1(buffer_seq[0], buffer_seq_sz[0], buffer_name[0]);
	    }
	    else {

        	bool stop = false;

	        vector<thread> workers; // need to keep track of threads so we can join them

	        mutex mutex_file;

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	        	buffer_seq[t] = new char[thread_seq_buf_sz];
	        	buffer_name[t] = new uint64_t[thread_name_buf_sz];

	            workers.emplace_back(

	                [&, t]{

	                    while (true) {

	                        {
	                            unique_lock<mutex> lock(mutex_file);

	                            if (stop) return;

	                            stop = reading_function1(fp, buffer_seq[t], buffer_seq_sz[t], buffer_name[t]);
	                        }

	                        worker_function1(buffer_seq[t], buffer_seq_sz[t], buffer_name[t]);
	                    }
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
    	}

    	for (auto& um : dbg){

    		UnitigData* ud = um.getData();
    		PairID& hap_ids = ud->get_hapID();

			/*if (ud->getPhasedKmerCoverage(um) < opt.min_cov_vertices) hap_ids.clear(); // Not enough k-mer coverage for the phased reads -> delete
			else*/ hap_ids.runOptimize();

			//ud->resetCoverage();
		}

        for (size_t t = 0; t < opt.nb_threads; ++t){

        	delete[] buffer_seq[t];
        	delete[] buffer_name[t];
        }

        delete[] buffer_seq;
        delete[] buffer_name;
        delete[] buffer_seq_sz;

        fp.close();
	};

	if (!v_phase_filenames.empty() && !v_read_filenames.empty()){

		loadPhasing(v_phase_filenames);

		if (mapHapReads) runHapMapping(v_read_filenames);
	}
}

pair<BlockedBloomFilter, unordered_set<uint64_t, CustomHashUint64_t>> buildBBF(const vector<string>& v_filenames_in, const Correct_Opt& opt, const size_t k, const size_t g, const bool long_reads) {

	size_t nb_unique_kmers = 0;
	size_t nb_non_unique_kmers = 0;

	// Estimate parameters of the Blocked Bloom Filter
	{
	    KmerStream_Build_opt kms_opt;

	    kms_opt.threads = opt.nb_threads;
	    kms_opt.verbose = opt.verbose;
	    kms_opt.files = v_filenames_in;
	    kms_opt.k = k;
	    kms_opt.g = g;
	    kms_opt.q = 0;

	    KmerStream kms(kms_opt);

	    nb_unique_kmers = max(1UL, kms.KmerF0());
	    nb_non_unique_kmers = max(1UL, nb_unique_kmers - min(nb_unique_kmers, kms.Kmerf1()));

	    if (opt.verbose) cout << "Estimated " << nb_unique_kmers << " unique k-mers and " << nb_non_unique_kmers << " non-unique k-mers." << endl;
	}

	BlockedBloomFilter bf_uniq(nb_unique_kmers, opt.nb_bits_unique_kmers_bf);
	BlockedBloomFilter bf_non_uniq(nb_non_unique_kmers, opt.nb_bits_non_unique_kmers_bf);

	unordered_set<uint64_t, CustomHashUint64_t> name_hset;

	FileParser fp(v_filenames_in);

    string s;

    size_t len_read = 0;
    size_t pos_read = 0;

    const bool multi_threaded = (opt.nb_threads != 1);

    // Main worker thread
    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz) {

        const char* str_end = seq_buf + seq_buf_sz;

        while (seq_buf < str_end) { // for each input

            const int len = strlen(seq_buf);

            toUpperCase(seq_buf, len);

            KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;
        	minHashIterator<RepHash> it_min(seq_buf, len, k, g, RepHash(), true);

            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

                const pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

                it_min += (p_.second - it_min.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

                const uint64_t min_hr = it_min.getHash();

                if (!bf_uniq.insert(p_.first, min_hr, multi_threaded)) bf_non_uniq.insert(p_.first, min_hr, multi_threaded);
            }

            seq_buf += len + 1;
        }
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz) {

        size_t file_id = 0;

        const size_t sz_buf = opt.buffer_sz - k;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_buf) {

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k){

                	//if (pos_read == 0) name_hset.insert(XXH64(fp.getNameString(), strlen(fp.getNameString()), opt.h_seed));

                    if ((opt.buffer_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(seq_buf + seq_buf_sz, s_str + pos_read, opt.buffer_sz - seq_buf_sz - 1);

                        seq_buf[opt.buffer_sz - 1] = '\0';

                        pos_read += sz_buf - seq_buf_sz;
                        seq_buf_sz = opt.buffer_sz;

                        break;
                    }
                    else {

                        strcpy(seq_buf + seq_buf_sz, s_str + pos_read);

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else return true;
        }

        return false;
    };

    {
        bool stop = false;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_file;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    char* buffer_seq = new char[opt.buffer_sz]();

                    size_t buffer_seq_sz = 0;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_file);

                            if (stop) break;

                            stop = reading_function(buffer_seq, buffer_seq_sz);
                        }

                        worker_function(buffer_seq, buffer_seq_sz);
                    }

                    delete[] buffer_seq;
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    fp.close();

	return pair<BlockedBloomFilter, unordered_set<uint64_t, CustomHashUint64_t>>(bf_non_uniq, name_hset);
}

BlockedBloomFilter buildBBF(const CompactedDBG<UnitigData>& dbg) {

	BlockedBloomFilter bbf(dbg.nbKmers(), 14);

	for (const auto& um : dbg) {

		const string unitig = um.referenceUnitigToString();

        KmerHashIterator<RepHash> it_kmer_h(unitig.c_str(), unitig.length(), dbg.getK()), it_kmer_h_end;
        minHashIterator<RepHash> it_min(unitig.c_str(), unitig.length(), dbg.getK(), dbg.getG(), RepHash(), true);

        for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

            const pair<uint64_t, int> p = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

            it_min += (p.second - it_min.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

            bbf.insert(p.first, it_min.getHash(), false);
        }
	}

	return bbf;
}

string retrieveMissingReads(const Correct_Opt& opt){

    Correct_Opt opt_pass_lr(opt);

    opt_pass_lr.k = opt_pass_lr.small_k;

    opt_pass_lr.filename_seq_in.clear();
    opt_pass_lr.filename_seq_in.insert(opt_pass_lr.filename_seq_in.end(), opt_pass_lr.filenames_long_in.begin(), opt_pass_lr.filenames_long_in.end());
    opt_pass_lr.filename_seq_in.insert(opt_pass_lr.filename_seq_in.end(), opt_pass_lr.filenames_helper_long_in.begin(), opt_pass_lr.filenames_helper_long_in.end());

    CompactedDBG<UnitigData> dbg_lr(opt_pass_lr.k);

	const size_t k = dbg_lr.getK();
	const size_t g = dbg_lr.getG();

    std::atomic<size_t> nb_reads_added;

    nb_reads_added = 0;

	string filename_out_extra_sr = opt.filename_long_out + "_extra_sr.fasta";

	if (opt.verbose) cout << "Ratatosk::retrieveMissingReads(): Creating index of short reads" << endl;

	const pair<BlockedBloomFilter, unordered_set<uint64_t, CustomHashUint64_t>> p_bf_sr = buildBBF(opt.filename_seq_in, opt, k, g, false);

	if (opt.verbose) cout << "Ratatosk::retrieveMissingReads(): Creating index of long reads" << endl;

    dbg_lr.build(opt_pass_lr);

    if (opt.verbose) cout << "Ratatosk::retrieveMissingReads(): Adding coverage to long read graph" << endl;

    const size_t nb_km_lr = dbg_lr.nbKmers();

    if (nb_km_lr != 0){

	    BlockedBloomFilter bf_lr(max(nb_km_lr, static_cast<size_t>(1)), opt.nb_bits_unique_kmers_bf);

	    auto addKmersToBBF = [&](CompactedDBG<UnitigData>::iterator it_s, CompactedDBG<UnitigData>::iterator it_e, const bool multi_threaded) {

	        while (it_s != it_e){

	        	const string unitig = it_s->referenceUnitigToString();

	            KmerHashIterator<RepHash> it_kmer_h(unitig.c_str(), unitig.length(), k), it_kmer_h_end;
	        	minHashIterator<RepHash> it_min(unitig.c_str(), unitig.length(), k, g, RepHash(), true);

	            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

	                const pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

	                it_min += (p_.second - it_min.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

	                const uint64_t min_hr = it_min.getHash();

	                bf_lr.insert(p_.first, min_hr, multi_threaded);
	            }

	            ++it_s;
	        }
	    };

	    const size_t slice = 1024;

	    CompactedDBG<UnitigData>::iterator it_s = dbg_lr.begin();
	    CompactedDBG<UnitigData>::iterator it_e = dbg_lr.end();

	    if ((opt.nb_threads == 1) || (dbg_lr.size() < slice)) addKmersToBBF(it_s, it_e, false);
	    else {

	        vector<thread> workers;

	        mutex mutex_graph;

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	            workers.emplace_back(

	                [&]{

	                	CompactedDBG<UnitigData>::iterator it_s_l;
	                	CompactedDBG<UnitigData>::iterator it_e_l;

	                	while (true){

	                		{
	                			unique_lock<mutex> lock(mutex_graph);

	                			if (it_s == it_e) return;
	                			else {

	                				it_s_l = it_s;
	                				it_e_l = it_s;

	                				for (size_t i = 0; (i < slice) && (it_e_l != it_e); ++i) ++it_e_l;

	                				it_s = it_e_l;
	                			}
	                		}

	                		addKmersToBBF(it_s_l, it_e_l, true);
	                	}
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
	    }

	    dbg_lr.clear();

		if (opt.verbose) cout << "Ratatosk::retrieveMissingReads(): Querying full short read set for missing reads" << endl;

	    ofstream outfile;
	    ostream out(0);

	    string s;
	    string n;

	    mutex mutex_file_in, mutex_file_out;

		FileParser fp(opt.filenames_short_all);

	    outfile.open(filename_out_extra_sr.c_str());
	    out.rdbuf(outfile.rdbuf());

	    // Main worker thread
	    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz, const vector<string>& v_read_names) {

	    	size_t i = 0;
	    	size_t l_nb_reads_added = 0;

	        const char* str_end = seq_buf + seq_buf_sz;

	        while (seq_buf < str_end) { // for each input

	            const int len = strlen(seq_buf);
	            const int len_km = len - k + 1;

	            //const uint64_t name_h = XXH64(v_read_names[i].c_str(), v_read_names[i].length(), opt.h_seed);

	            size_t nb_km = 0;
	            size_t remaining_len_km = len_km;

	            toUpperCase(seq_buf, len); // Put characters in upper case

				/*if (p_bf_sr.second.find(name_h) != p_bf_sr.second.end()) {

	        		++l_nb_reads_added;

	        		unique_lock<mutex> lock(mutex_file_out);

	        		out << '>' << v_read_names[i] << '\n' << string(seq_buf) << '\n';
				}
				else {*/

		            KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;
		        	minHashIterator<RepHash> it_min(seq_buf, len, k, g, RepHash(), true);

		            for (; (it_kmer_h != it_kmer_h_end) && (nb_km < opt.min_nb_km_unmapped) && (nb_km + remaining_len_km >= opt.min_nb_km_unmapped); ++it_kmer_h) {

		                const pair<uint64_t, int> p = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

		                it_min += (p.second - it_min.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

		                const uint64_t min_hr = it_min.getHash();

		                nb_km += static_cast<size_t>(!p_bf_sr.first.contains(p.first, min_hr) && bf_lr.contains(p.first, min_hr));
		                remaining_len_km = len_km - p.second;
		            }

		            if (nb_km >= opt.min_nb_km_unmapped) {

		        		++l_nb_reads_added;

		        		unique_lock<mutex> lock(mutex_file_out);

		        		out << '>' << v_read_names[i] << '\n' << string(seq_buf) << '\n';
		            }

					/*KmerHashIterator<RepHash> it_kmer_h(seq_buf, len, k), it_kmer_h_end;
		        	minHashIterator<RepHash> it_min(seq_buf, len, k, g, RepHash(), true);

		        	size_t max_len_run = 0;

		            for (; (it_kmer_h != it_kmer_h_end); ++it_kmer_h) {

		                const pair<uint64_t, int> p = *it_kmer_h; // <k-mer hash, k-mer position in sequence>

		                it_min += (p.second - it_min.getKmerPosition()); //If one or more k-mer were jumped because contained non-ACGT char.

		                const uint64_t min_hr = it_min.getHash();

	                	if (p_bf_sr.first.contains(p.first, min_hr) || bf_lr.contains(p.first, min_hr)) ++nb_km;
	                	else {

	                		max_len_run = max(max_len_run, nb_km);
	                		nb_km = 0;
	                	}
		            }

		            if (max_len_run >= opt.min_nb_km_unmapped) {

		        		++l_nb_reads_added;

		        		unique_lock<mutex> lock(mutex_file_out);

		        		out << '>' << v_read_names[i] << '\n' << string(seq_buf) << '\n';
		            }*/
	        	//}

	            seq_buf += len + 1;
	            ++i;
	        }

	        return l_nb_reads_added;
	    };

	    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz, vector<string>& v_read_names) {

	        size_t file_id = 0;

	        seq_buf_sz = 0;

	        while (seq_buf_sz < opt.buffer_sz) {

	        	const bool isNewRead = (s.length() == 0);

	            if (!isNewRead || fp.read(s, file_id)) {

	                const size_t len_s = s.length();

	                if (len_s >= k){

	                	if (isNewRead) n = fp.getNameString();

	                    if ((opt.buffer_sz - seq_buf_sz) < (len_s + 1)) break;
	                    else {

	                    	v_read_names.push_back(move(n));

	                        strcpy(seq_buf + seq_buf_sz, s.c_str());

	                        seq_buf_sz += len_s + 1;

	                        s.clear();
	                    }
	                }
	                else s.clear();
	            }
	            else return true;
	        }

	        return false;
	    };

	    {
	        bool stop = false;

	        size_t nb_processed = 0;

	        vector<thread> workers; // need to keep track of threads so we can join them

	        for (size_t t = 0; t < opt.nb_threads; ++t){

	            workers.emplace_back(

	                [&]{

	                    char* buffer_seq = new char[opt.buffer_sz]();

	                    size_t buffer_seq_sz = 0;

	                    vector<string> v_read_names;

	                    while (true) {

	                        {
	                            unique_lock<mutex> lock(mutex_file_in);

	                            if (stop) break;

	                            const size_t nb_processed_before = nb_processed;

	                            stop = reading_function(buffer_seq, buffer_seq_sz, v_read_names);
	                            nb_processed += v_read_names.size();

	                            if ((nb_processed / 1000000) != (nb_processed_before / 1000000)) {

	                            	cout << "Ratatosk::retrieveMissingReads(): " << nb_processed << " short reads queried." << endl;
	                            }
	                        }

	                        nb_reads_added += worker_function(buffer_seq, buffer_seq_sz, v_read_names);
	                        v_read_names.clear();
	                    }

	                    delete[] buffer_seq;
	                }
	            );
	        }

	        for (auto& t : workers) t.join();
	    }

		if (opt.verbose) cout << "Ratatosk::retrieveMissingReads(): Added " << nb_reads_added << " short reads to dataset." << endl;

	    fp.close();
	    outfile.close();
	}

	if (nb_reads_added == 0){

		remove(filename_out_extra_sr.c_str());
		filename_out_extra_sr = "";
	} 

    return filename_out_extra_sr;
}

/*size_t estimateCoverage(const CompactedDBG<UnitigData>& dbg) {

	size_t tot_cov = 0, nb_km = 0;

	const size_t k = dbg.getK();

	for (const auto& um_start : dbg) {

		if (um_start.getData()->isBranching()) {

			vector<const_UnitigMap<UnitigData>> v_um;

			bool branches = false, is_snp = true;

			for (const auto& um_succ : um_start.getSuccessors()){

				v_um.push_back(um_succ);

				is_snp = (is_snp && (um_succ.len == k));
				branches = (branches || um_succ.getData()->isBranching());
			}

			if ((v_um.size() >= 2) && is_snp && !branches) { // At least two non-branching successors of with k k-mers

				const_UnitigMap<UnitigData> um_end;

				bool is_simple_bubble = true;

				for (const auto& um_succ : v_um){

					for (const auto& um_succ_succ : um_succ.getSuccessors()){

						if (um_end.isEmpty) um_end = um_succ_succ;
						else is_simple_bubble = is_simple_bubble && um_succ_succ.isSameReferenceUnitig(um_end) && (um_succ_succ.strand == um_end.strand);
					}
				}

				if (is_simple_bubble) {

					for (const auto& um_succ : v_um){

						nb_km += um_succ.len;
						tot_cov += um_succ.getData()->getCoverage();
					}
				}
			}
		}
	}

	return ((nb_km == 0) ? 0 : (tot_cov / nb_km));
}*/

size_t estimateHaplotypeCoverage(const CompactedDBG<UnitigData>& dbg) {

	size_t tot_cov = 0, nb_km = 0;

	const size_t k = dbg.getK();

	for (const auto& um_start : dbg) {

		if (um_start.getSuccessors().cardinality() > 1) {

			vector<const_UnitigMap<UnitigData>> v_um;

			bool branches = false;

			for (const auto& um_succ : um_start.getSuccessors()){

				v_um.push_back(um_succ);

				branches = (branches || (um_succ.getSuccessors().cardinality() > 1) || (um_succ.getPredecessors().cardinality() > 1));
			}

			if ((v_um.size() >= 2) && !branches) { // At least two non-branching successors of with k k-mers

				const_UnitigMap<UnitigData> um_end;

				bool is_simple_bubble = true;

				for (const auto& um_succ : v_um){

					for (const auto& um_succ_succ : um_succ.getSuccessors()){

						if (um_end.isEmpty) um_end = um_succ_succ;
						else is_simple_bubble = is_simple_bubble && um_succ_succ.isSameReferenceUnitig(um_end) && (um_succ_succ.strand == um_end.strand);
					}
				}

				if (is_simple_bubble) {

					for (const auto& um_succ : v_um){

						nb_km += um_succ.len;
						tot_cov += um_succ.getData()->getCoverage();
					}
				}
			}
		}
	}

	return ((nb_km == 0) ? 0 : (tot_cov / nb_km));
}

/*void expandConnectedComponent(CompactedDBG<UnitigData>& dbg, const size_t id_mapped, const size_t id_visited) { 

	for (const auto& um : dbg){

		// If this is a mapped read that has not been visited
		if (um.getData()->getConnectedComp() == id_mapped) {

			queue<UnitigMap<UnitigData>> q;

			q.push(um);

			while (!q.empty()){

				const UnitigMap<UnitigData> traversed_um = q.front();

				q.pop();

				if (traversed_um.getData()->getConnectedComp() == id_mapped){

					traversed_um.getData()->setConnectedComp(id_visited);

					for (const auto& traversed_um_succ : traversed_um.getSuccessors()) q.push(traversed_um_succ);
					for (const auto& traversed_um_pred : traversed_um.getPredecessors()) q.push(traversed_um_pred);
				}
			}
		}
	}

	return;
}*/

/*Roaring* createPartitions(CompactedDBG<UnitigData>& dbg, const vector<Kmer>& v_centroids, const double max_path_len, const size_t nb_threads, const bool verbose) {

	auto expandPartition = [&](const UnitigMap<UnitigData>& src, const size_t id_partition){

		size_t sum = 0;

		KmerHashTable<uint8_t> km_h;
		queue<UnitigMap<UnitigData>> q1;

		q1.push(src);
		km_h.insert(src.getUnitigHead(), 0);

		while (!q1.empty()){

			const UnitigMap<UnitigData> um = q1.front();

			q1.pop();
			um.getData()->setConnectedComp(id_partition);

			sum += um.len;

			if (sum >= max_path_len) break;

			for (const auto& um_succ : um.getSuccessors()){

				if ((um_succ.getData()->getConnectedComp() == 0) && km_h.insert(um_succ.getUnitigHead(), 0).second) q1.push(um_succ);
			}

			for (const auto& um_pred : um.getPredecessors()){

				if ((um_pred.getData()->getConnectedComp() == 0) && km_h.insert(um_pred.getUnitigHead(), 0).second) q1.push(um_pred);
			}
		}

		return sum;
	};

	auto mergeNeighbors = [&](const UnitigMap<UnitigData>& src, vector<size_t>& v_part_sz){

		const size_t id_part_src = src.getData()->getConnectedComp();

		size_t min_bp_part = 0xffffffffffffffffULL;
		size_t min_id_part = id_part_src;

		Roaring part_neighbors;

		KmerHashTable<uint8_t> km_h;
		queue<UnitigMap<UnitigData>> q1;

		q1.push(src);
		km_h.insert(src.getUnitigHead(), 0);

		while (!q1.empty()){

			const UnitigMap<UnitigData> um = q1.front();

			q1.pop();

			for (const auto& um_succ : um.getSuccessors()){

				if (km_h.insert(um_succ.getUnitigHead(), 0).second){

					const size_t id_part_succ = um_succ.getData()->getConnectedComp();

					if (id_part_succ == id_part_src) q1.push(um_succ);
					else part_neighbors.add(id_part_succ);
				}
			}

			for (const auto& um_pred : um.getPredecessors()){

				if (km_h.insert(um_pred.getUnitigHead(), 0).second){

					const size_t id_part_pred = um_pred.getData()->getConnectedComp();

					if (id_part_pred == id_part_src) q1.push(um_pred);
					else part_neighbors.add(id_part_pred);
				}
			}
		}

		for (const uint32_t part_id : part_neighbors){

			if (v_part_sz[part_id] < min_bp_part) {

				min_bp_part = v_part_sz[part_id];
				min_id_part = part_id;
			}
		}

		if (min_bp_part != 0xffffffffffffffffULL){

			v_part_sz[min_id_part] += v_part_sz[id_part_src];
			v_part_sz[id_part_src] = 0xffffffffffffffffULL;

			q1.push(src);
			km_h.clear();
			km_h.insert(src.getUnitigHead(), 0);

			while (!q1.empty()){

				const UnitigMap<UnitigData> um = q1.front();

				q1.pop();
				um.getData()->setConnectedComp(min_id_part);

				for (const auto& um_succ : um.getSuccessors()){

					if ((um_succ.getData()->getConnectedComp() == id_part_src) && (km_h.insert(um_succ.getUnitigHead(), 0).second)) q1.push(um_succ);
				}

				for (const auto& um_pred : um.getPredecessors()){

					if ((um_pred.getData()->getConnectedComp() == id_part_src) && (km_h.insert(um_pred.getUnitigHead(), 0).second)) q1.push(um_pred);
				}
			}

			return true;
		}

		return false;
	};

	auto getNeighbors = [&](const UnitigMap<UnitigData>& src){

		PairID p_id_part;
		KmerHashTable<uint8_t> km_h;
		queue<UnitigMap<UnitigData>> q1, q2;

		Roaring neighbors;

		const size_t id_part_src = src.getData()->getConnectedComp();

		q1.push(src);
		km_h.insert(src.getUnitigHead(), 0);
		neighbors.add(id_part_src);

		while (!q1.empty()){

			const UnitigMap<UnitigData> um = q1.front();

			q1.pop();

			p_id_part |= um.getData()->getPairID();

			for (const auto& um_succ : um.getSuccessors()){

				if (km_h.insert(um_succ.getUnitigHead(), 0).second){

					if (um_succ.getData()->getConnectedComp() == id_part_src) q1.push(um_succ);
					else q2.push(um_succ);
				}
			}

			for (const auto& um_pred : um.getPredecessors()){

				if (km_h.insert(um_pred.getUnitigHead(), 0).second){

					if (um_pred.getData()->getConnectedComp() == id_part_src) q1.push(um_pred);
					else q2.push(um_pred);
				}
			}
		}

		TinyBloomFilter<uint32_t> tbf_id_part(p_id_part.cardinality(), 14);

		for (const uint32_t id_part : p_id_part) tbf_id_part.insert(id_part);

		while (!q2.empty()){

			const UnitigMap<UnitigData> um = q2.front();

			const PairID& p_id_um = um.getData()->getPairID();
			const size_t part_id_um = um.getData()->getConnectedComp();

			const bool is_already_neighbor = neighbors.contains(part_id_um);

			q2.pop();

			if (p_id_um.isEmpty() || is_already_neighbor || hasEnoughSharedPairID(tbf_id_part, p_id_part, p_id_um, 1)){

				if (!is_already_neighbor) neighbors.add(part_id_um);

				for (const auto& um_succ : um.getSuccessors()){

					if (km_h.insert(um_succ.getUnitigHead(), 0).second) q2.push(um_succ);
				}

				for (const auto& um_pred : um.getPredecessors()){

					if (km_h.insert(um_pred.getUnitigHead(), 0).second) q2.push(um_pred);
				}
			}
		}

		return neighbors;
	};

	vector<Kmer> v_centroids_tmp;
	vector<size_t> v_part_sz;

	// Create partitions from centroids that expand 500 bp in each direction
	size_t id_part = 1;

	const size_t k = dbg.getK();

	v_part_sz.push_back(0);

	if (v_centroids.empty()){

		for (const auto& um : dbg){

			if (um.getData()->getConnectedComp() == 0){

				const size_t sz = expandPartition(um, id_part);

				++id_part;

				v_part_sz.push_back(sz);
				v_centroids_tmp.push_back(um.getUnitigHead());
			}
		}
	}
	else {

		for (const auto km_centroid : v_centroids){

			const UnitigMap<UnitigData> um = dbg.find(km_centroid);

			if (!um.isEmpty && (um.getData()->getConnectedComp() == 0)){

				UnitigMap<UnitigData> um_tmp = um;

				um_tmp.dist = 0;
				um_tmp.len = um.size - k + 1;

				const size_t sz = expandPartition(um_tmp, id_part);

				++id_part;

				v_part_sz.push_back(sz);
			}
		}
	}

	if (verbose) cout << "Ratatosk::createPartitions(): Created " << (id_part - 1) << " partitions." << endl;

	const vector<Kmer>& v_centroids_ref = v_centroids.empty() ? v_centroids_tmp : v_centroids;

	size_t nb_too_small = 0;
	size_t prev_nb_too_small = 0xffffffffffffffffULL;

	for (const auto sz_part : v_part_sz) nb_too_small += static_cast<size_t>(sz_part < max_path_len);

	while ((nb_too_small > 1) && (nb_too_small != prev_nb_too_small)) {

		if (verbose) cout << "Ratatosk::createPartitions(): Merging " << nb_too_small << " partitions." << endl;

		prev_nb_too_small = nb_too_small;

		for (const auto km_centroid : v_centroids_ref){

			const UnitigMap<UnitigData> um = dbg.find(km_centroid);

			if (!um.isEmpty){

				const size_t id_part_um = um.getData()->getConnectedComp();

				if ((id_part_um != 0) && (v_part_sz[id_part_um] < max_path_len)) {

					UnitigMap<UnitigData> um_tmp = um;

					um_tmp.dist = 0;
					um_tmp.len = um.size - k + 1;

					mergeNeighbors(um_tmp, v_part_sz);
				}
			}
		}

		nb_too_small = 0;

		for (const auto sz_part : v_part_sz) nb_too_small += static_cast<size_t>(sz_part < max_path_len);
	}

	Roaring* adjacency = new Roaring[id_part + 1];

	size_t processed_centroid = 0;

	if (nb_threads == 1){

		for (const auto km_centroid : v_centroids_ref){

			const UnitigMap<UnitigData> um = dbg.find(km_centroid);

			if (!um.isEmpty){

				const size_t id_part_um = um.getData()->getConnectedComp();

				if ((id_part_um != 0) && adjacency[id_part_um].isEmpty()){

					UnitigMap<UnitigData> um_tmp = um;

					um_tmp.dist = 0;
					um_tmp.len = um.size - k + 1;

					adjacency[id_part_um] = getNeighbors(um_tmp);

					++processed_centroid;

					if (verbose) cout << "Ratatosk::createPartitions(): Processed " << processed_centroid << " centroids." << endl;
				}
			}
		}
	}
	else {

		const size_t id_unreachable = id_part + 1;

        vector<thread> workers; // need to keep track of threads so we can join them

        LockGraph lck_g(nb_threads * 1024);

        std::atomic<size_t> i;
        std::atomic<size_t> processed_centroid;

        i = 0;
        processed_centroid = 0;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                	while (true){

	                	const size_t l_i = i++;

	                	if (l_i >= v_centroids_ref.size()) return;

						const UnitigMap<UnitigData> um = dbg.find(v_centroids_ref[l_i]);

						if (!um.isEmpty){

							const size_t id_part_um = um.getData()->getConnectedComp();

							if (id_part_um != 0){

								lck_g.lock_unitig(id_part_um);

								if (adjacency[id_part_um].isEmpty()) {

									adjacency[id_part_um].add(id_unreachable); // Add fake partition ID that cannot be reached, no other thread cannot process this partition anymore

									lck_g.unlock_unitig(id_part_um);

									UnitigMap<UnitigData> um_tmp = um;

									um_tmp.dist = 0;
									um_tmp.len = um.size - k + 1;

									Roaring r = getNeighbors(um_tmp);

									lck_g.lock_unitig(id_part_um);

									adjacency[id_part_um] = move(r);

									lck_g.unlock_unitig(id_part_um);

									const size_t l_processed_centroid = ++processed_centroid;

									if (verbose && (l_processed_centroid % 100 == 0)) cout << "Ratatosk::createPartitions(): Processed " << l_processed_centroid << " centroids." << endl;
								}
								else lck_g.unlock_unitig(id_part_um);
							}
						}
                	}
                }
			);
        }

        for (auto& t : workers) t.join();

        for (size_t i = 0; i < id_part + 1; ++i) adjacency[i].remove(id_unreachable); // Remove fake partition ID that cannot be reached
	}

	for (size_t i = 0; i < id_part + 1; ++i) adjacency[i].add(0); // Unvisited vertices are always considered to be nearby

	return adjacency;
}*/

/*size_t detectSTRs(CompactedDBG<UnitigData>& dbg, const size_t max_len_cycle, const size_t max_cov_km, const size_t nb_threads, const bool verbose){

	const size_t k = dbg.getK();

	if (verbose) cout << "Ratatosk::detectSTRs(): Scanning graph for unitigs in short cycles (STRs)" << endl;

	auto comparePathsDist = [](const pair<size_t, UnitigMap<UnitigData>>& a, const pair<size_t, UnitigMap<UnitigData>>& b){

		return a.first > b.first;
	};

	auto BFS = [&](const UnitigMap<UnitigData>& um, const bool strand, const size_t max_len_cycle){

		if (um.size >= max_len_cycle) return false;

		KmerHashTable<uint8_t> km_h;
		queue<pair<size_t, UnitigMap<UnitigData>>> q;

		UnitigMap<UnitigData> start_um(um);

		start_um.dist = 0;
		start_um.len = start_um.size - k + 1;
		start_um.strand = strand;

		const Kmer start_head = um.getMappedHead();
		const size_t max_len_path = max_len_cycle - k;

		q.push({0, start_um});
		km_h.insert(start_head, 0);

		while (!q.empty()){

			const pair<size_t, UnitigMap<UnitigData>> traversed_path = q.front();

			q.pop();

			for (const auto& um_succ : traversed_path.second.getSuccessors()){

				const Kmer head_succ = um_succ.getMappedHead();
				const pair<KmerHashTable<uint8_t>::iterator, bool> p_it_succ = km_h.insert(head_succ, 0);

				if (head_succ == start_head) return true;
				if (p_it_succ.second && (traversed_path.first + 1 <= max_len_path)) q.push({traversed_path.first + 1, um_succ});
			}
		}

		return false;
	};

	auto detectCycle = [&](const UnitigMap<UnitigData>& um, const bool strand, const size_t max_len_cycle) {

		if (um.size >= max_len_cycle) return false;
		//if (!BFS(um, strand, max_len_cycle)) return false;

		KmerHashTable<size_t> km_h;
		priority_queue<pair<size_t, UnitigMap<UnitigData>>, vector<pair<size_t, UnitigMap<UnitigData>>>, decltype(comparePathsDist)> q_traversal(comparePathsDist);

		UnitigMap<UnitigData> start_um(um);

		start_um.dist = 0;
		start_um.len = start_um.size - k + 1;
		start_um.strand = strand;

		const Kmer start_head = um.getMappedHead();

		q_traversal.push({0, start_um});

		while (!q_traversal.empty()){

			const pair<size_t, UnitigMap<UnitigData>> traversed_path = q_traversal.top();

			q_traversal.pop();

			for (const auto& um_succ : traversed_path.second.getSuccessors()){

				//if (um_succ.getData()->getCoverage() != 0){

					const Kmer head_succ = um_succ.getMappedHead();
					const KmerHashTable<size_t>::iterator it_succ = km_h.find(head_succ);
					const size_t dist_succ = (it_succ == km_h.end()) ? 0xffffffffffffffffULL : *it_succ; 

					if (dist_succ > (traversed_path.first + (um_succ.size - k + 1))){

						if (dist_succ == 0xffffffffffffffffULL) km_h.insert(head_succ, traversed_path.first + (um_succ.size - k + 1));
						else *it_succ = traversed_path.first + (um_succ.size - k + 1);

						if ((traversed_path.first + k) <= max_len_cycle) q_traversal.push({traversed_path.first + (um_succ.size - k + 1), um_succ});
					}

					if (head_succ == start_head){

						const KmerHashTable<size_t>::const_iterator it_succ = km_h.find(head_succ);

						if ((it_succ != km_h.end()) && (*it_succ <= max_len_cycle)) return true;
					}
				//}
			}
		}

		return false;
	};

	if (nb_threads == 1){

		size_t nb_repeats = 0;

		for (const auto& um : dbg){

			if (um.size <= max_len_cycle){ // Motif is short enough to be an STR or poly-A/C/G/T sequence

				const double km_frequency = um.getData()->getKmerCoverage(um);
				const double seq_entropy = getEntropy(um.referenceUnitigToString().c_str(), um.size);

				if (((km_frequency >= max_cov_km) && (seq_entropy <= 1)) || (seq_entropy <= 0.5)) { // Sequence has very high coverage or low sequence entropy

					if (detectCycle(um, true, max_len_cycle)){

						um.getData()->setCycle();
						++nb_repeats;
					}
					else if (detectCycle(um, false, max_len_cycle)) {

						um.getData()->setCycle();
						++nb_repeats;
					}
				}
			}
		}

		return nb_repeats;
	}
	else {

		atomic<size_t> nb_repeats;

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_graph;

        CompactedDBG<UnitigData>::iterator it_dbg = dbg.begin();

        nb_repeats = 0;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                	UnitigMap<UnitigData> um;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_graph);

                            if (it_dbg == dbg.end()) return;

                            um = *it_dbg;

                            ++it_dbg;
                        }

                        bool cycle = false;

						if (um.size <= max_len_cycle){ // Motif is short enough to be an STR or poly-A/C/G/T sequence

							const double km_frequency = um.getData()->getKmerCoverage(um);
							const double seq_entropy = getEntropy(um.referenceUnitigToString().c_str(), um.size);

							// Sequence has very high coverage and low-ish sequence entropy or just low sequence entropy
							if (((km_frequency >= max_cov_km) && (seq_entropy <= 1)) || (seq_entropy <= 0.5)) {

								if (detectCycle(um, true, max_len_cycle)){

									um.getData()->setCycle();
									++nb_repeats;
								}
								else if (detectCycle(um, false, max_len_cycle)) {

									um.getData()->setCycle();
									++nb_repeats;
								}
							}
						}
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        return nb_repeats;
	}
}*/

// Detect unitigs which are in a short cycle
void detectShortCycles(CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt) {

	auto detectShortCycle = [&](const UnitigMap<UnitigData>& um, const bool fw) {

		queue<pair<UnitigMap<UnitigData>, PairID>> q_um;

		unordered_set<Kmer, KmerHash> us_km;

		q_um.push({um, um.getData()->getPairID().toPairID()});

		while (!q_um.empty()) {

			const pair<UnitigMap<UnitigData>, PairID> p_um = q_um.front();

			q_um.pop();

			if (fw) {

				for (const auto& um_succ : p_um.first.getSuccessors()) {

					if (um.isSameReferenceUnitig(um_succ) && (um_succ.strand == um.strand)) {

						um.getData()->setVisitStatus(true);

						return true;
					}
					else if (us_km.insert(um_succ.getMappedHead()).second) { // Insertion took place, vertex was never visited before

						PairID p_int = p_um.second & um_succ.getData()->getPairID().toPairID();

						if (p_int.cardinality() >= opt.min_cov_vertices) q_um.push({um_succ, move(p_int)});
					}
				}
			}
			else {

				for (const auto& um_pred : p_um.first.getPredecessors()) {

					if (um.isSameReferenceUnitig(um_pred) && (um_pred.strand == um.strand)) {

						um.getData()->setVisitStatus(true);

						return true;
					}
					else if (us_km.insert(um_pred.getMappedHead()).second) { // Insertion took place, vertex was never visited before

						PairID p_int = p_um.second & um_pred.getData()->getPairID().toPairID();

						if (p_int.cardinality() >= opt.min_cov_vertices) q_um.push({um_pred, move(p_int)});
					}
				}
			}
		}

		return false;
	};

	if (opt.verbose) cout << "Ratatosk::detectShortCycles(): Detecting unitigs in short cycles" << endl;

	if (opt.nb_threads == 1){

		size_t count_short_cycles = 0;

	    for (const auto& um : dbg) {

	    	um.getData()->setVisitStatus(false);

	    	count_short_cycles += static_cast<size_t>(detectShortCycle(um, true) || detectShortCycle(um, false));
	    }

	    if (opt.verbose) cout << "Ratatosk::detectShortCycles(): Found " << count_short_cycles << " in short cycles. " << endl;
	}
	else {

		const size_t slice = 1024;

        bool stop = false;

        atomic<size_t> count_short_cycles(0);

        vector<thread> workers; // need to keep track of threads so we can join them

        mutex mutex_graph;

        CompactedDBG<UnitigData>::iterator it_c = dbg.begin();
        CompactedDBG<UnitigData>::iterator it_e = dbg.end();

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

					CompactedDBG<UnitigData>::iterator l_it_c, l_it_e;

                    while (true) {

                        {
                            unique_lock<mutex> lock(mutex_graph);

                            if (stop) return;

                            l_it_c = it_c;

					        size_t i = 0;

					        while ((i < slice) && (it_c != it_e)){

					        	++i;
					        	++it_c;
					        }

                            stop = (it_c == it_e);
							l_it_e = it_c;
                        }

                        {
							size_t l_count_short_cycles = 0;

							while (l_it_c != l_it_e){

								l_it_c->getData()->setVisitStatus(false);

								l_count_short_cycles += static_cast<size_t>(detectShortCycle(*l_it_c, true) || detectShortCycle(*l_it_c, false));
								++l_it_c;
							}

							count_short_cycles += l_count_short_cycles;
                        }
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        if (opt.verbose) cout << "Ratatosk::detectShortCycles(): Found " << count_short_cycles << " in short cycles. " << endl;
	}
}
