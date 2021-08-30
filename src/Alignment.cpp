#include "Alignment.hpp"

pair<int, int> selectBestAlignment(const vector<Path<UnitigData>>& candidates, const char* ref, const size_t len_ref){

	const size_t k = Kmer::k;

	double best_edit_dist;
	int best_end_loc = -1;
	int best_cand_id = -1;

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, len_ref, config);

	best_edit_dist = static_cast<double>(align.editDistance) / max(candidate.length(), len_ref);
	best_end_loc = align.endLocations[0];
	best_cand_id = 0;

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();

		const size_t max_len_norm = max(candidate.length(), len_ref);

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, len_ref, config);

		if ((align.editDistance >= 0) && ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist)) {

			best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
			best_end_loc = align.endLocations[0];
			best_cand_id = i;
		}

		edlibFreeAlignResult(align);

		p_cand_before = move(p_cand_curr);
	}

	return {best_cand_id, best_end_loc};
}

pair<int, int> selectBestPrefixAlignment(const char* ref, const size_t ref_len, const vector<Path<UnitigData>>& candidates, const double cut_threshold_norm_edit){

	const size_t k = Kmer::k;

	double best_edit_dist;

	int best_cand_id = -1;
	int best_end_loc = -1;

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / candidate.length();
	best_cand_id = 0;
	best_end_loc = align.endLocations[0];

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = candidate.length();

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if ((align.editDistance >= 0) && ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist)) {

			best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
			best_cand_id = i;
			best_end_loc = align.endLocations[0];
		}

		edlibFreeAlignResult(align);

		p_cand_before = move(p_cand_curr);
	}

	if ((cut_threshold_norm_edit > 0.0) && (best_edit_dist > cut_threshold_norm_edit)) return {-1, -1};

	return {best_cand_id, best_end_loc};
}

pair<int, int> selectBestPrefixAlignment(const char* ref, const size_t ref_len, const vector<const Path<UnitigData>*>& candidates, const double cut_threshold_norm_edit){

	const size_t k = Kmer::k;

	double best_edit_dist;

	int best_cand_id = -1;
	int best_end_loc = -1;

	Path<UnitigData>::PathOut p_cand_before = candidates[0]->toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));
	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / candidate.length();
	best_cand_id = 0;
	best_end_loc = align.endLocations[0];

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i]->toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = candidate.length();

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if ((align.editDistance >= 0) && ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist)) {

			best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
			best_cand_id = i;
			best_end_loc = align.endLocations[0];
		}

		edlibFreeAlignResult(align);

		p_cand_before = move(p_cand_curr);
	}

	if ((cut_threshold_norm_edit > 0.0) && (best_edit_dist > cut_threshold_norm_edit)) return {-1, -1};

	return {best_cand_id, best_end_loc};
}

pair<int, pair<int, int>> selectBestPrefixAlignment2(const char* ref, const size_t ref_len, const vector<Path<UnitigData>>& candidates, const bool get_end_location){

	const size_t k = Kmer::k;

	double best_edit_dist;

	int best_cand_id = -1;
	int best_end_loc_target = -1;
	int best_end_loc_src = -1;

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / max(candidate.length(), ref_len);
	best_cand_id = 0;
	best_end_loc_target = align.endLocations[0];

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = max(candidate.length(), ref_len);

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if ((align.editDistance >= 0) && ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist)) {

			best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
			best_cand_id = i;
			best_end_loc_target = align.endLocations[0];
		}

		edlibFreeAlignResult(align);

		p_cand_before = move(p_cand_curr);
	}

	if (get_end_location){

		char* cigar = nullptr;

		size_t cigar_len = 0;

		const string best_candidate = candidates[best_cand_id].toString();
		const size_t max_len_norm = max(best_candidate.length(), ref_len);

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		align = edlibAlign(best_candidate.c_str(), best_candidate.length(), ref, ref_len, config);
		cigar = edlibAlignmentToCigar(align.alignment, align.alignmentLength, EDLIB_CIGAR_STANDARD);
		cigar_len = strlen(cigar);

		size_t cigar_pos = 0;
		size_t prev_cigar_pos = 0;
		size_t target_pos = 0;
		size_t query_pos = 0;

		while ((cigar_pos != cigar_len) && (target_pos <= best_end_loc_target)){

			if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

				if (cigar[cigar_pos] == 'M') { //match

					const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

					query_pos += cigar_l;
					target_pos += cigar_l;
				}
				else if ((cigar[cigar_pos] == 'I') || (cigar[cigar_pos] == 'S')) query_pos += atoi(&cigar[prev_cigar_pos]); //insertion or soft-clipping
				else if (cigar[cigar_pos] == 'D') target_pos += atoi(&cigar[prev_cigar_pos]);  //deletion

				prev_cigar_pos = cigar_pos + 1;
			}

			best_end_loc_src = query_pos - 1;

			++cigar_pos;
		}

		free(cigar);

		edlibFreeAlignResult(align);
	}

	return {best_cand_id, {best_end_loc_target, best_end_loc_src}};
}

bool test(const string& query_seq, const size_t query_len, const string& qual_seq, const size_t qual_len, const string& ref_seq, const size_t ref_len) {

	EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, edlib_iupac_alpha, sz_edlib_iupac_alpha);
	EdlibAlignResult align_fw = edlibAlign(query_seq.c_str(), query_len, ref_seq.c_str(), ref_len, config);

	char* cigar  = edlibAlignmentToCigar(align_fw.alignment, align_fw.alignmentLength, EDLIB_CIGAR_STANDARD);

	const size_t cigar_len = strlen(cigar);

	edlibFreeAlignResult(align_fw);

	size_t prev_cigar_pos = 0;
	size_t cigar_pos = 0;
	size_t query_pos = 0;
	size_t ref_pos = 0;

	bool ret = false;

	cout << "===" << endl;

	while ((cigar_pos != cigar_len) && (ref_pos < ref_len)){

		if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

			if (cigar[cigar_pos] == 'M') { //match

				const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

				query_pos += cigar_l;
				ref_pos += cigar_l;
			}
			else if ((cigar[cigar_pos] == 'I') || (cigar[cigar_pos] == 'S')) {

				if (atoi(&cigar[prev_cigar_pos]) >= 5) cout << "I" << atoi(&cigar[prev_cigar_pos]) << endl;

				query_pos += atoi(&cigar[prev_cigar_pos]); //insertion or soft-clipping
			}
			else if (cigar[cigar_pos] == 'D') {

				if (atoi(&cigar[prev_cigar_pos]) >= 2) {

					cout << "D" << atoi(&cigar[prev_cigar_pos]) << " " << query_pos << " " << ref_pos << " " << query_len << " " << ref_len << endl;
					ret = true;
				}

				ref_pos += atoi(&cigar[prev_cigar_pos]);  //deletion
			}

			prev_cigar_pos = cigar_pos + 1;
		}

		++cigar_pos;
	}

	free(cigar);

	if (ret) {

		cout << query_seq << endl;
		cout << qual_seq << endl;
	}

	return ret;
}

pair<string, string> generateConsensus(const ResultCorrection* fw_s, const ResultCorrection* bw_s, const string& ref_seq, const double max_norm_edit_distance){

	if ((bw_s->getNbCorrectedPosOldSeq() == 0) && (fw_s->getNbCorrectedPosOldSeq() != 0)) return {fw_s->getSequence(), fw_s->getQuality()};
	else if ((fw_s->getNbCorrectedPosOldSeq() == 0) && (bw_s->getNbCorrectedPosOldSeq() != 0)) return {bw_s->getSequence(), bw_s->getQuality()};
	else if (fw_s->getNbCorrectedPosOldSeq() + bw_s->getNbCorrectedPosOldSeq() == 0) return {string(), string()};

	if (bw_s->getNbCorrectedPosOldSeq() > fw_s->getNbCorrectedPosOldSeq()) swap(fw_s, bw_s);

	EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, edlib_iupac_alpha, sz_edlib_iupac_alpha);

	EdlibAlignResult align_fw = edlibAlign(fw_s->getSequence().c_str(), fw_s->getSequence().length(), ref_seq.c_str(), ref_seq.length(), config);
	EdlibAlignResult align_bw = edlibAlign(bw_s->getSequence().c_str(), bw_s->getSequence().length(), ref_seq.c_str(), ref_seq.length(), config);

	const double align_dist_norm_fw = static_cast<double>(align_fw.editDistance) / max(fw_s->getSequence().length(), ref_seq.length());
	const double align_dist_norm_bw = static_cast<double>(align_bw.editDistance) / max(bw_s->getSequence().length(), ref_seq.length());

	if ((max_norm_edit_distance > 0.0) && ((align_dist_norm_fw > max_norm_edit_distance) || (align_dist_norm_bw > max_norm_edit_distance))) {

		edlibFreeAlignResult(align_fw);
		edlibFreeAlignResult(align_bw);

		if ((align_dist_norm_fw > max_norm_edit_distance) && (align_dist_norm_bw > max_norm_edit_distance)) return {string(), string()};
		if (align_dist_norm_fw > max_norm_edit_distance) return {bw_s->getSequence(), bw_s->getQuality()};
		
		return {fw_s->getSequence(), fw_s->getQuality()};
	}

	char* cigar_fw  = edlibAlignmentToCigar(align_fw.alignment, align_fw.alignmentLength, EDLIB_CIGAR_STANDARD);
	char* cigar_bw  = edlibAlignmentToCigar(align_bw.alignment, align_bw.alignmentLength, EDLIB_CIGAR_STANDARD);

	const size_t cigar_fw_len = strlen(cigar_fw), cigar_bw_len = strlen(cigar_bw);

	size_t cigar_fw_pos = 0, cigar_bw_pos = 0;
	size_t prev_cigar_fw_pos = 0, prev_cigar_bw_pos = 0;
	size_t fw_s_pos = 0, bw_s_pos = 0;
	size_t ref_pos_fw = 0, ref_pos_bw = 0;

	size_t i = 0;

	edlibFreeAlignResult(align_fw);
	edlibFreeAlignResult(align_bw);

	stringstream ss;
	stringstream sq;

	auto moveIntoCIGAR = [](const size_t start, const size_t end, size_t& cigar_pos, size_t& prev_cigar_pos, size_t& query_pos, size_t& ref_pos, const char* cigar, const size_t cigar_len){

		size_t read_pos_start = query_pos, read_pos_end = query_pos;

		while ((cigar_pos != cigar_len) && (ref_pos < start)){

			if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

				if (cigar[cigar_pos] == 'M') { //match

					const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

					if (ref_pos + cigar_l > start){ // We span the boundary of the start position on reference

						read_pos_start = query_pos + (start - ref_pos);
						break;
					}

					query_pos += cigar_l;
					ref_pos += cigar_l;
				}
				else if ((cigar[cigar_pos] == 'I') || (cigar[cigar_pos] == 'S')) query_pos += atoi(&cigar[prev_cigar_pos]); //insertion or soft-clipping
				else if (cigar[cigar_pos] == 'D') ref_pos += atoi(&cigar[prev_cigar_pos]);  //deletion

				prev_cigar_pos = cigar_pos + 1;
				read_pos_start = query_pos;
			}

			++cigar_pos;
		}

		read_pos_end = read_pos_start;

		while ((cigar_pos != cigar_len) && (ref_pos < end)){

			if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

				if (cigar[cigar_pos] == 'M') { //match

					const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

					if (ref_pos + cigar_l > end) return make_pair(make_pair(read_pos_start, query_pos + (end - ref_pos)), end); // We span the boundary of the start position on reference

					query_pos += cigar_l;
					ref_pos += cigar_l;
				}
				else if ((cigar[cigar_pos] == 'I') || (cigar[cigar_pos] == 'S')) query_pos += atoi(&cigar[prev_cigar_pos]); //insertion or soft-clipping
				else if (cigar[cigar_pos] == 'D') ref_pos += atoi(&cigar[prev_cigar_pos]);  //deletion

				prev_cigar_pos = cigar_pos + 1;
				read_pos_end = query_pos;
			}

			++cigar_pos;
		}

		return make_pair(make_pair(read_pos_start, read_pos_end), ref_pos);
	};

	while (i < ref_seq.length()){

		int64_t len_fw = fw_s->getLengthCorrectedRegion(i);
		int64_t len_bw = bw_s->getLengthCorrectedRegion(i);

		pair<pair<size_t, size_t>, size_t> pos_read;

		if ((len_fw + len_bw) <= 0){

			len_fw = fw_s->getLengthUncorrectedRegion(i);
			len_bw = bw_s->getLengthUncorrectedRegion(i);

			if ((len_fw > len_bw) || (len_fw <= 0)) len_fw = -1;
			else len_bw = -1;
		}
		
		if (len_fw >= len_bw){

			pos_read = moveIntoCIGAR(i, i + len_fw, cigar_fw_pos, prev_cigar_fw_pos, fw_s_pos, ref_pos_fw, cigar_fw, cigar_fw_len);

			if (pos_read.first.second > pos_read.first.first){

				ss << fw_s->getSequence().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);
				sq << fw_s->getQuality().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);
			}
		}
		else {

			pos_read = moveIntoCIGAR(i, i + len_bw, cigar_bw_pos, prev_cigar_bw_pos, bw_s_pos, ref_pos_bw, cigar_bw, cigar_bw_len);

			if (pos_read.first.second > pos_read.first.first){

				ss << bw_s->getSequence().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);
				sq << bw_s->getQuality().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);
			}
		}

		i = pos_read.second;
	}

	free(cigar_fw);
	free(cigar_bw);

	if (max_norm_edit_distance > 0.0){

		const string ss_str = ss.str();

		EdlibAlignResult align = edlibAlign(ss_str.c_str(), ss_str.length(), ref_seq.c_str(), ref_seq.length(), edlibDefaultAlignConfig());

		const double align_dist_norm = static_cast<double>(align.editDistance) / max(ss_str.length(), ref_seq.length());

		edlibFreeAlignResult(align);

		if (align_dist_norm > max_norm_edit_distance) return {fw_s->getSequence(), fw_s->getQuality()};
	}
	
	return {ss.str(), sq.str()};
}

/*vector<pair<size_t, char>> filterAmbiguity(	const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const WeightsPairID& w_pid,
											const string& query, const vector<pair<size_t, char>>& v_ambiguity) {

	auto isValidAmb = [&](const string& q_sub, const char c, const size_t pos) {

		bool valid = false;
		bool found = false;

		string q = q_sub;

		q[pos] = c;

		for (KmerIterator it_km(q.c_str()), it_km_end; (it_km != it_km_end) && !valid; ++it_km) {

			const const_UnitigMap<UnitigData> um = dbg.findUnitig(q.c_str(), it_km->second, q.length());

			if (!um.isEmpty){

				const UnitigData* ud = um.getData();
				const PairID& p_ids = ud->getPairID();

				valid = hasEnoughSharedPairID(p_ids, w_pid.all_pids, opt.min_cov_vertices);
				found = true;
				it_km += um.len - 1;
			}
		}

		return (valid || !found);
	};

	vector<pair<size_t, char>> v_amb_out;

	for (const auto& p_amb : v_ambiguity) {

		const size_t pos_buff = (p_amb.first < (opt.k - 1)) ? 0 : (p_amb.first - opt.k + 1);
		const size_t len_buff = min(p_amb.first + opt.k, query.length()) - pos_buff;
		const size_t pos_snp_buff = p_amb.first - pos_buff;

		const string q_sub = query.substr(pos_buff, len_buff);

		bool is_a = false, is_c = false, is_g = false, is_t = false;

		getAmbiguityRev(p_amb.second, is_a, is_c, is_g, is_t);

		if (is_a) is_a = isValidAmb(q_sub, 'A', pos_snp_buff);
		if (is_c) is_c = isValidAmb(q_sub, 'C', pos_snp_buff);
		if (is_g) is_g = isValidAmb(q_sub, 'G', pos_snp_buff);
		if (is_t) is_t = isValidAmb(q_sub, 'T', pos_snp_buff);

		if (static_cast<size_t>(is_a) + static_cast<size_t>(is_c) + static_cast<size_t>(is_g) + static_cast<size_t>(is_t) > 1) v_amb_out.push_back({p_amb.first, getAmbiguity(is_a, is_c, is_g, is_t)});
	}

	return v_amb_out;
}*/

void fixAmbiguity(	const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const WeightsPairID& w_pid,
					string& query, string& quality, 
					const char* ref_seq, const char* ref_qual, const size_t ref_len,
					const uint64_t hap_id, const vector<pair<size_t, char>>& v_ambiguity){

	if (!v_ambiguity.empty()) {

		const size_t quality_len = quality.length();
		const size_t query_len = query.length();

		const char q_max_corr = getQual(1.0, opt.out_qual);
		const char q_min_corr = getQual(0.0, opt.out_qual);

		const char q_min_conf_corr = getQual(opt.min_confidence_snp_corr);

		const size_t k = dbg.getK();

		const char c_noCorrect = 'X';

		const uint64_t undetermined_hap_id = 0xffffffffffffffffULL;

		auto isPairIDconsistent = [&](const size_t pos_km, const char snp) {

			const size_t pos_buff = (pos_km < (k - 1)) ? 0 : (pos_km - k + 1);
			const size_t len_buff = min(pos_km + k, query_len) - pos_buff;
			const size_t pos_snp_buff = pos_km - pos_buff;

			string q_sub = query.substr(pos_buff, len_buff);

			q_sub[pos_snp_buff] = snp;

			bool isConsistent = false;

			for (KmerIterator it_km(q_sub.c_str()), it_km_end; (it_km != it_km_end) && !isConsistent; ++it_km) {

				const const_UnitigMap<UnitigData> um = dbg.findUnitig(q_sub.c_str(), it_km->second, q_sub.length());

				if (!um.isEmpty){

					const UnitigData* ud = um.getData();

					isConsistent = (getNumberSharedPairID(ud->getPairID(), w_pid.all_pids, opt.min_cov_vertices) >= opt.min_cov_vertices);

					it_km += um.len - 1;
				}
			}

			return isConsistent;
		};

		struct compare_ambiguity {

			bool operator()(const pair<size_t, char>& a, const pair<size_t, char>& b) const {

				return ((a.first < b.first) || ((a.first == b.first) && (a.second < b.second)));
			}
		};

		string query_tmp = query;

		unordered_map<size_t, char> m_ambiguity_safe;
		unordered_map<size_t, char> m_ambiguity_all;

		for (const auto p : v_ambiguity){

			if (quality[p.first] < q_min_conf_corr) {

				m_ambiguity_safe.insert(p);
				query_tmp[p.first] = p.second;
			}
		}

		m_ambiguity_all = m_ambiguity_safe;

		const EdlibAlignConfig config_ambiguous = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		EdlibAlignResult align = edlibAlign(query_tmp.c_str(), query_len, ref_seq, ref_len, config_ambiguous);

		char* cigar = edlibAlignmentToCigar(align.alignment, align.alignmentLength, EDLIB_CIGAR_STANDARD);

		const size_t cigar_len = strlen(cigar);

		size_t cigar_pos = 0;
		size_t prev_cigar_pos = 0;
		size_t query_pos = 0;
		size_t target_pos = align.startLocations[0];

		while (cigar_pos != cigar_len){

			if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

				const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

				if (cigar[cigar_pos] == 'M') { //match

					for (size_t q_pos = query_pos, t_pos = target_pos; q_pos < query_pos + cigar_l; ++q_pos, ++t_pos){

						if (!isDNA(query_tmp[q_pos])) { // SNP candidate detected from corrected sequence

							if (!isDNA(ref_seq[t_pos])) { // Ref already has ambiguity char. at that location -> SNP from prev. correction pass.

								unordered_map<size_t, char>::iterator it = m_ambiguity_safe.find(q_pos);

								if (it != m_ambiguity_safe.end()) it->second = c_noCorrect;
							}
							else if (quality[q_pos] >= q_min_corr)  {

								bool alpha_q[4] = {false, false, false, false};
								bool alpha_t[4] = {false, false, false, false};

								getAmbiguityRev(ref_seq[t_pos], alpha_t[0], alpha_t[1], alpha_t[2], alpha_t[3]);
								getAmbiguityRev(query_tmp[q_pos], alpha_q[0], alpha_q[1], alpha_q[2], alpha_q[3]);

								if ((alpha_q[0] && alpha_t[0]) || (alpha_q[1] && alpha_t[1]) || (alpha_q[2] && alpha_t[2]) || (alpha_q[3] && alpha_t[3])) {

									unordered_map<size_t, char>::iterator it = m_ambiguity_safe.find(q_pos);

									if (it != m_ambiguity_safe.end()) it->second = ref_seq[t_pos];
								}
							}

							unordered_map<size_t, char>::iterator it = m_ambiguity_all.find(q_pos);

							if (it != m_ambiguity_all.end()) it->second = ref_seq[t_pos];
						}
						else if (!isDNA(ref_seq[t_pos])) {

							if (quality[q_pos] < q_min_conf_corr) {

								m_ambiguity_safe.insert({q_pos, c_noCorrect});
								m_ambiguity_all.insert({q_pos, ref_seq[t_pos]});
							}
							else {

								bool alpha_q[4] = {false, false, false, false};
								bool alpha_t[4] = {false, false, false, false};

								getAmbiguityRev(ref_seq[t_pos], alpha_t[0], alpha_t[1], alpha_t[2], alpha_t[3]);
								getAmbiguityRev(query_tmp[q_pos], alpha_q[0], alpha_q[1], alpha_q[2], alpha_q[3]);

								if (!(alpha_q[0] && alpha_t[0]) && !(alpha_q[1] && alpha_t[1]) && !(alpha_q[2] && alpha_t[2]) && !(alpha_q[3] && alpha_t[3])) {

									m_ambiguity_safe.insert({q_pos, c_noCorrect});
									m_ambiguity_all.insert({q_pos, ref_seq[t_pos]});
								}
							}
						}
					}

					query_pos += cigar_l;
					target_pos += cigar_l;
				}
				else if (cigar[cigar_pos] == 'I') { //insertion

					for (size_t q_pos = query_pos; q_pos < query_pos + cigar_l; ++q_pos){

						if (!isDNA(query_tmp[q_pos])) {

							unordered_map<size_t, char>::iterator it_s = m_ambiguity_safe.find(q_pos);
							unordered_map<size_t, char>::iterator it_a = m_ambiguity_all.find(q_pos);

							if ((it_s != m_ambiguity_safe.end()) && (it_a != m_ambiguity_all.end())){

								it_a->second = it_s->second;
								it_s->second = c_noCorrect;
							}
						}
					}

					query_pos += cigar_l;
				}
				else if (cigar[cigar_pos] == 'S') query_pos += cigar_l; //soft-clipping
				else if (cigar[cigar_pos] == 'D') target_pos += cigar_l; //deletion

				prev_cigar_pos = cigar_pos + 1;
			}

			++cigar_pos;
		}

		free(cigar);
		edlibFreeAlignResult(align);

		set<pair<size_t, char>, compare_ambiguity> s_ambiguity;

		for (auto& p : m_ambiguity_safe){

			if (isDNA(p.second)){

				//if (isPairIDconsistent(p.first, p.second)){

					const size_t pos_buff = (p.first < (k - 1)) ? 0 : (p.first - k + 1);
					const size_t len_buff = min(p.first + k, query_len) - pos_buff;
					const size_t pos_snp_buff = p.first - pos_buff;

					string q_sub = query.substr(pos_buff, len_buff);

					q_sub[pos_snp_buff] = p.second;

					for (KmerIterator it_km(q_sub.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

						const pair<Kmer, int>& p_km = *it_km;
						const const_UnitigMap<UnitigData> um = dbg.findUnitig(q_sub.c_str(), p_km.second, q_sub.length());

						if (!um.isEmpty) {

							const PairID& hap_um = um.getData()->get_hapID();

							if ((hap_id == undetermined_hap_id) || isValidHap(hap_um, hap_id)) {

								const_UnitigMap<UnitigData> um_tmp = um;

								um_tmp.dist = 0;
								um_tmp.len = um_tmp.size - k + 1;

								const string unitig_seq = um_tmp.mappedSequenceToString();
								const vector<pair<size_t, char>> v_amb = um_tmp.getData()->get_ambiguity_char(um_tmp);

								size_t pos_snp_unitig = (pos_snp_buff - p_km.second) + um.dist;

								if (!um.strand) pos_snp_unitig = um.size - pos_snp_unitig - 1;

								for (const auto p_amb : v_amb){

									int64_t pos = static_cast<int64_t>(p_amb.first);

									if (pos <= pos_snp_unitig) pos = static_cast<int64_t>(p.first) - static_cast<int64_t>(pos_snp_unitig - pos);
									else pos = static_cast<int64_t>(p.first) + static_cast<int64_t>(pos - pos_snp_unitig);

									if ((pos >= 0) && (pos < query_len) && (pos != p.first)) {

										const unordered_map<size_t, char>::const_iterator it = m_ambiguity_safe.find(pos);

										if ((it != m_ambiguity_safe.end()) && !isDNA(it->second)) s_ambiguity.insert(pair<size_t, char>(pos, unitig_seq[p_amb.first]));
									}
								}
		                	}

			                it_km += um.len - 1;
						}
					}
				//}
				//else p.second = query[p.first];
			}
		}

		{
			const vector<pair<size_t, char>> v(s_ambiguity.begin(), s_ambiguity.end());

			for (int64_t i = 0; i < v.size(); ++i){

				if (((i == 0) || (v[i].first != v[i-1].first)) && ((i == v.size()-1) || (v[i].first != v[i+1].first))){

					unordered_map<size_t, char>::iterator it_s = m_ambiguity_safe.find(v[i].first);

					if (it_s != m_ambiguity_safe.end()){

						char c_amb = it_s->second;

						bool is_a_q = false, is_c_q = false, is_g_q = false, is_t_q = false;
						bool is_a_t = false, is_c_t = false, is_g_t = false, is_t_t = false;

						getAmbiguityRev(c_amb, is_a_t, is_c_t, is_g_t, is_t_t);
						getAmbiguityRev(v[i].second, is_a_q, is_c_q, is_g_q, is_t_q);

						if ((is_a_q && is_a_t) || (is_c_q && is_c_t) || (is_g_q && is_g_t) || (is_t_q && is_t_t)) it_s->second = v[i].second;
					}
				}
			}
		}

		for (const auto& p : m_ambiguity_safe){

			if ((p.second == c_noCorrect) || (quality[p.first] < q_min_corr)) {

				const unordered_map<size_t, char>::const_iterator it_a = m_ambiguity_all.find(p.first);

				if (it_a != m_ambiguity_all.end()){

					bool validHap = ((hap_id == undetermined_hap_id) || (p.second == c_noCorrect) || !isDNA(it_a->second));

					if (!validHap) {

						const size_t pos_buff = (p.first < (k - 1)) ? 0 : (p.first - k + 1);
						const size_t len_buff = min(p.first + k, query_len) - pos_buff;
						const size_t pos_snp_buff = p.first - pos_buff;

						string q_sub = query.substr(pos_buff, len_buff);

						q_sub[pos_snp_buff] = it_a->second;

						for (KmerIterator it_km(q_sub.c_str()), it_km_end; (it_km != it_km_end) && !validHap; ++it_km) {

							const const_UnitigMap<UnitigData> um = dbg.findUnitig(q_sub.c_str(), it_km->second, q_sub.length());

							if (!um.isEmpty){

								validHap = isValidHap(um.getData()->get_hapID(), hap_id);
								it_km += um.len - 1;
							}
						}
					}

					if (validHap) {

						query_tmp[p.first] = it_a->second;
						quality[p.first] = q_max_corr;
					}
					else query_tmp[p.first] = query[p.first]; // Replacing the SNP with base from raw sequence lead to incorrect phasing, keep the correction
					
				}
			}
			else if (!isDNA(p.second)) query_tmp[p.first] = query[p.first];
			else query_tmp[p.first] = p.second; // HERE -> Why not putting raw sequence base?
		}

		query = move(query_tmp);
	}
}

string fixSNPs(const Correct_Opt& opt, const CompactedDBG<UnitigData>& dbg, const string& s) {

	struct amb_t {

		bool dna_aplha[4];

		amb_t() {

			clear();
		}

		amb_t(const char amb_c) {

			clear();

			getAmbiguityRev(amb_c, dna_aplha[0], dna_aplha[1], dna_aplha[2], dna_aplha[3]);
		}

		void clear() {

			dna_aplha[0] = false;
			dna_aplha[1] = false;
			dna_aplha[2] = false;
			dna_aplha[3] = false;
		}
	};

	const size_t limit_nb_km_cand = 64;

	string out_s = s;

	if (s.length() >= opt.k) {

		map<size_t, amb_t> m_amb;

		for (size_t i = 0; i < s.length(); ++i){

			if (!isDNA(s[i])) m_amb.insert({i, amb_t(s[i])});
		}

		for (size_t i = 0; i < s.length(); ++i){

			if (!isDNA(out_s[i])){

				const size_t min_pos_amb = (i < (opt.k - 1)) ? 0 : (i - opt.k + 1);
				const size_t len_buff_amb = min(i + opt.k, s.length()) - min_pos_amb;
				const size_t pos_amb_buff = i - min_pos_amb;

				const string s_sub = out_s.substr(min_pos_amb, len_buff_amb);

				vector<pair<size_t, amb_t>> v_amb;

				map<size_t, amb_t>::const_iterator it_s = m_amb.lower_bound(min_pos_amb);
				map<size_t, amb_t>::const_iterator it_e = m_amb.upper_bound(min_pos_amb + len_buff_amb);

				v_amb.insert(v_amb.end(), it_s, it_e);

				if (!v_amb.empty()) {

					size_t nb_km_cand = 1;

					for (const auto& p_amb : v_amb) { // Compute number of k-mers to query

						nb_km_cand *= static_cast<size_t>(p_amb.second.dna_aplha[0]) + static_cast<size_t>(p_amb.second.dna_aplha[1]) + static_cast<size_t>(p_amb.second.dna_aplha[2]) + static_cast<size_t>(p_amb.second.dna_aplha[3]);

						if (nb_km_cand >= limit_nb_km_cand) break;
					}

					if (nb_km_cand < limit_nb_km_cand) { // Not too many possibilities

						unordered_set<char> s_amb_cand;

						for (auto& p_amb : v_amb) p_amb.first -= min_pos_amb; //Set ambiguity position coordinates within local substring rather than whole read

						for (size_t j = 0; (j < v_amb.size() * 4) && (s_amb_cand.size() <= 1); ++j) {

							string l_s_sub = s_sub;

							bool valid = true;

							for (size_t pos_v_amb = 0; (pos_v_amb < v_amb.size()) && valid; ++pos_v_amb) {

								const size_t subpos_v_amb = (j >> (pos_v_amb << 1)) & 0x3ULL;

								if (v_amb[pos_v_amb].second.dna_aplha[subpos_v_amb]) l_s_sub[v_amb[pos_v_amb].first] = alpha[subpos_v_amb];
								else valid = false;
							}

							if (valid && (s_amb_cand.find(l_s_sub[pos_amb_buff]) == s_amb_cand.end())) {

					            KmerHashIterator<RepHash> it_kmer_h(l_s_sub.c_str(), len_buff_amb, opt.k), it_kmer_h_end;

					            for (; it_kmer_h != it_kmer_h_end; ++it_kmer_h) {

					                const std::pair<uint64_t, int> p_ = *it_kmer_h; // <k-mer hash, k-mer position in sequence>
						            const const_UnitigMap<UnitigData> um = dbg.findUnitig(l_s_sub.c_str(), p_.second, len_buff_amb);

						            if (!um.isEmpty) {

						            	s_amb_cand.insert(l_s_sub[pos_amb_buff]);
						            	break;
						            }
						        }
							}
						}

						if (s_amb_cand.size() == 1){

							out_s[i] = *(s_amb_cand.begin());

							m_amb.erase(i); // These ambiguity characters are now invalid
						}
					}
				}
			}
		}
	}

	return out_s;
}

pair<int, int> selectBestSubstringAlignment(const char* ref, const size_t ref_len, const vector<Path<UnitigData>>& candidates, const double cut_threshold_norm_edit){

	const size_t k = Kmer::k;

	double best_edit_dist;

	int best_cand_id = -1;
	int best_end_loc = -1;

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha));

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / candidate.length();
	best_cand_id = 0;
	best_end_loc = align.endLocations[0];

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = candidate.length();

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, edlib_iupac_alpha, sz_edlib_iupac_alpha);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if  ((align.editDistance >= 0) && ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist)) {

			best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
			best_cand_id = i;
			best_end_loc = align.endLocations[0];
		}

		edlibFreeAlignResult(align);

		p_cand_before = move(p_cand_curr);
	}

	if ((cut_threshold_norm_edit > 0.0) && (best_edit_dist > cut_threshold_norm_edit)) return {-1, -1};

	return {best_cand_id, best_end_loc};
}

vector<pair<size_t, const_UnitigMap<UnitigData>>> keep_non_overlap(const Correct_Opt& opt/*, const WeightsPairID& w_pid*/, const char* ref, const vector<pair<size_t, const_UnitigMap<UnitigData>>>& v){

	struct var_info_t {

		size_t pos_s;
		size_t pos_e;

		bool keep;

		PairID pos_v;
		
		//PairID pid;

		set<Kmer> s_km;

		//const const_UnitigMap<UnitigData>* prev_um;

		var_info_t() : pos_s(0), pos_e(0), keep(true)/*, prev_um(nullptr)*/ {}
	};

	const size_t ref_len = strlen(ref);

	map<size_t, var_info_t> m_var;

	PairID pos_out;

	vector<pair<size_t, const_UnitigMap<UnitigData>>> v_out;

	for (size_t i = 0; i < v.size(); ++i) {

		const pair<size_t, const_UnitigMap<UnitigData>>& p_um = v[i];

		const string km_ref = string(ref + p_um.first, opt.k);
		const string km_query = p_um.second.mappedSequenceToString();

		const size_t l_match_pref = cstrMatch(km_ref.c_str(), km_query.c_str());

		uint8_t type_var = 0;
		uint8_t mis_ins = getAmbiguityIndex('.');

		if (l_match_pref + cstrMatch(km_ref.c_str() + l_match_pref + 1, km_query.c_str() + l_match_pref + 1) == (opt.k-1)) {

			type_var = 1; // mismatch
			mis_ins = getAmbiguityIndex(km_query[l_match_pref]); // mismatching character
		}
		else if (l_match_pref + cstrMatch(km_ref.c_str() + l_match_pref, km_query.c_str() + l_match_pref + 1) == (opt.k-1)){

			type_var = 2; // insertion
			mis_ins = getAmbiguityIndex(km_query[l_match_pref]); // inserted character
		}
		else if (l_match_pref + cstrMatch(km_ref.c_str() + l_match_pref + 1, km_query.c_str() + l_match_pref) == (opt.k-1)) type_var = 3; // deletion

		if ((type_var != 0) && (l_match_pref != 0) && (l_match_pref != opt.k-1)){

			const size_t pos = p_um.first + l_match_pref;
			const size_t key = (pos << 16) | ((static_cast<size_t>(mis_ins) & 0x00000000000000ffULL) << 8) | (static_cast<size_t>(type_var) & 0x00000000000000ffULL);

			const const_UnitigMap<UnitigData>& um = p_um.second;
			const SharedPairID& spid = um.getData()->getPairID();

			pair<map<size_t, var_info_t>::iterator, bool> p_it = m_var.insert(pair<size_t, var_info_t>(key, var_info_t()));

			var_info_t& var_info = p_it.first->second;

			if (p_it.second) {

				var_info.pos_s = p_um.first;
				var_info.pos_e = p_um.first + opt.k;

				//var_info.pid = spid.toPairID();
			}
			else {

				var_info.pos_s = min(var_info.pos_s, p_um.first);
				var_info.pos_e = max(var_info.pos_e, p_um.first + opt.k);

				//if (!v[i].second.isSameReferenceUnitig(*var_info.prev_um)) var_info.pid |= spid.toPairID();
			}

			/*if (*/var_info.s_km.insert(um.strand ? um.getUnitigHead() : um.getUnitigHead().twin())/*.second) var_info.pid |= spid.toPairID()*/;
			//var_info.prev_um = &(v[i].second);

			var_info.pos_v.add(i);
		}
	}

	for (auto& it1 : m_var) {

		if (it1.second.keep) {

			const size_t it1_pos = (it1.first >> 16);

			//const size_t lower_bound_pos = (it1.second.pos_s < opt.k-1) ? 0 : (it1.second.pos_s - opt.k + 1);
			//const size_t upper_bound_pos = ((it1.second.pos_e + opt.k) >= ref_len) ? ref_len : (it1.second.pos_e + opt.k);
			const size_t lower_bound_pos = (it1_pos < opt.k-1) ? 0 : (it1_pos - opt.k + 1);
			const size_t upper_bound_pos = ((it1_pos + opt.k) >= ref_len) ? ref_len : (it1_pos + opt.k);

			const size_t lower_bound_pos_key = lower_bound_pos << 16;
			const size_t upper_bound_pos_key = (upper_bound_pos << 16) + 0x000000000000ffffULL;

			map<size_t, var_info_t>::iterator it2 = m_var.lower_bound(lower_bound_pos_key);

			while (it1.second.keep && (it2 != m_var.end()) && (it2->first <= upper_bound_pos_key)) {

				const size_t it2_pos = (it2->first >> 16);

				const bool overlap1 = (it1_pos >= it2->second.pos_s) && (it1_pos < it2->second.pos_e);
				const bool overlap2 = (it2_pos >= it1.second.pos_s) && (it2_pos < it1.second.pos_e);

				// it2 kmer variant overlaps it1 variant
				if ((it1.first != it2->first) && (overlap1 || overlap2)) {

					bool sameUnitig = false;

					set<Kmer>::const_iterator it1_s = it1.second.s_km.begin();
					set<Kmer>::const_iterator it2_s = it2->second.s_km.begin();

					set<Kmer>::const_iterator it1_e = it1.second.s_km.end();
					set<Kmer>::const_iterator it2_e = it2->second.s_km.end();

					while ((it1_s != it1_e) && (it2_s != it2_e)) {

						if (*it1_s < *it2_s) ++it1_s;
						else if (*it2_s < *it1_s) ++it2_s;
						else {
						
							sameUnitig = true;
							break;
						}
					}

					if (!sameUnitig){

						/*const bool it1_safe = it1.second.keep && (getNumberSharedPairID(it1.second.pid, w_pid.weighted_pids, opt.min_cov_vertices) >= opt.min_cov_vertices);
						const bool it2_safe = it2->second.keep && (getNumberSharedPairID(it2->second.pid, w_pid.weighted_pids, opt.min_cov_vertices) >= opt.min_cov_vertices);

						it1.second.keep = (it1_safe && !it2_safe);
						it2->second.keep = (it2_safe && !it1_safe);*/

						it1.second.keep = false;
						it2->second.keep = false;

						/*{
							const size_t lower_bound_pos2 = it2->second.pos_s;
							const size_t upper_bound_pos2 = it2->second.pos_e;

							const char mis1 = getAmbiguityIndexRev(((it1.first >> 8) & 0x00000000000000ffULL));
							const char mis2 = getAmbiguityIndexRev(((it2->first >> 8) & 0x00000000000000ffULL));

							const string seq1 = string(ref + lower_bound_pos, upper_bound_pos - lower_bound_pos);
							const string seq2 = string(ref + lower_bound_pos2, upper_bound_pos2 - lower_bound_pos2);

							string seq_q1 = seq1;
							string seq_q2 = seq2;

							if ((it1.first & 0x00000000000000ffULL) == 1) seq_q1[it1_pos - lower_bound_pos] = mis1;
							else if ((it1.first & 0x00000000000000ffULL) == 2) seq_q1.insert(it1_pos - lower_bound_pos, 1, mis1);
							else seq_q1.erase(it1_pos - lower_bound_pos, 1);

							if ((it2->first & 0x00000000000000ffULL) == 1) seq_q2[it2_pos - lower_bound_pos2] = mis2;
							else if ((it2->first & 0x00000000000000ffULL) == 2) seq_q2.insert(it2_pos - lower_bound_pos2, 1, mis2);
							else seq_q2.erase(it2_pos - lower_bound_pos2, 1);

							cout << "1: " << seq1 << " " << seq_q1 << " " << lower_bound_pos << " " << upper_bound_pos << " " << it1_pos << " " << (it1.first & 0x00000000000000ffULL) << endl;
							cout << "2: " << seq2 << " " << seq_q2 << " " << lower_bound_pos2 << " " << upper_bound_pos2 << " " << it2_pos << " " << (it2->first & 0x00000000000000ffULL) << endl;
						}*/
					}
				}

				++it2;
			}

			if (it1.second.keep) {

				for (const uint32_t pos : it1.second.pos_v) pos_out.add(pos);
			}
		}
	}

	for (const uint32_t pos : pos_out) v_out.push_back(v[pos]);

	return v_out;
}