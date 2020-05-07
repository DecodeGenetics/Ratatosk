#include "Alignment.hpp"

pair<int, int> selectBestAlignment(const vector<Path<UnitigData>>& candidates, const char* ref, const size_t len_ref){

	const size_t k = Kmer::k;

	double best_edit_dist;
	int best_end_loc = -1;
	int best_cand_id = -1;

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, len_ref, config);

	best_edit_dist = static_cast<double>(align.editDistance) / max(candidate.length(), len_ref);
	best_end_loc = align.endLocations[0];
	best_cand_id = 0;

	for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();

		const size_t max_len_norm = max(candidate.length(), len_ref);

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, len_ref, config);

		if (align.editDistance >= 0){

			if ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist) {

				best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
				best_end_loc = align.endLocations[0];
				best_cand_id = i;

				for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);
			}
			else if ((static_cast<double>(align.editDistance) / max_len_norm) == best_edit_dist) {

				for (size_t j = 0; j < align.numLocations; ++j){

					if (align.endLocations[j] > best_end_loc){

						best_cand_id = i;
						best_end_loc = align.endLocations[j];
					}
				}
			}
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

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0));

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / candidate.length();
	best_cand_id = 0;
	best_end_loc = align.endLocations[0];

	for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = candidate.length();

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if (align.editDistance >= 0){

			if ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist){

				best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
				best_cand_id = i;
				best_end_loc = align.endLocations[0];

				for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);
			}
			else if ((static_cast<double>(align.editDistance) / max_len_norm) == best_edit_dist){

				for (size_t j = 0; j < align.numLocations; ++j){

					if (align.endLocations[j] > best_end_loc){

						best_cand_id = i;
						best_end_loc = align.endLocations[j];
					}
				}
			}
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

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0));
	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / candidate.length();
	best_cand_id = 0;
	best_end_loc = align.endLocations[0];

	for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i]->toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = candidate.length();

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if (align.editDistance >= 0){

			if ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist){

				best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
				best_cand_id = i;
				best_end_loc = align.endLocations[0];

				for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);
			}
			else if ((static_cast<double>(align.editDistance) / max_len_norm) == best_edit_dist){

				for (size_t j = 0; j < align.numLocations; ++j){

					if (align.endLocations[j] > best_end_loc){

						best_cand_id = i;
						best_end_loc = align.endLocations[j];
					}
				}
			}
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

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0));

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / max(candidate.length(), ref_len);
	best_cand_id = 0;
	best_end_loc_target = align.endLocations[0];
	
	for (size_t j = 1; j < align.numLocations; ++j) best_end_loc_target = max(best_end_loc_target, align.endLocations[j]);

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = max(candidate.length(), ref_len);

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_DISTANCE, NULL, 0);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if (align.editDistance >= 0){

			if ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist){

				best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
				best_cand_id = i;
				best_end_loc_target = align.endLocations[0];
				
				for (size_t j = 1; j < align.numLocations; ++j) best_end_loc_target = max(best_end_loc_target, align.endLocations[j]);
			}
			else if ((static_cast<double>(align.editDistance) / max_len_norm) == best_edit_dist){

				for (size_t j = 0; j < align.numLocations; ++j){

					if (align.endLocations[j] > best_end_loc_target){

						best_cand_id = i;
						best_end_loc_target = align.endLocations[j];
					}
				}
			}
		}

		edlibFreeAlignResult(align);

		p_cand_before = move(p_cand_curr);
	}

	if (get_end_location){

		char* cigar = nullptr;

		size_t cigar_len = 0;

		const string best_candidate = candidates[best_cand_id].toString();
		const size_t max_len_norm = max(best_candidate.length(), ref_len);

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0);
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

pair<string, string> generateConsensus(const ResultCorrection* fw_s, const ResultCorrection* bw_s, const string& ref_seq, const double max_norm_edit_distance){

	const bool haveQual = (!fw_s->getQuality().empty() && !bw_s->getQuality().empty());

	if ((bw_s->getNbCorrectedPosOldSeq() == 0) && (fw_s->getNbCorrectedPosOldSeq() != 0)) return {fw_s->getSequence(), fw_s->getQuality()};
	else if ((fw_s->getNbCorrectedPosOldSeq() == 0) && (bw_s->getNbCorrectedPosOldSeq() != 0)) return {bw_s->getSequence(), bw_s->getQuality()};
	else if (fw_s->getNbCorrectedPosOldSeq() + bw_s->getNbCorrectedPosOldSeq() == 0) return {string(), string()};

	if (bw_s->getNbCorrectedPosOldSeq() > fw_s->getNbCorrectedPosOldSeq()) swap(fw_s, bw_s);

	EdlibAlignConfig config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);

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
	size_t prev_cigar_fw_pos = 0, prev_cigar_bw_pos = 0;;
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

		if ((len_fw == 0) && (len_bw == 0)){

			len_fw = fw_s->getLengthUncorrectedRegion(i);
			len_bw = bw_s->getLengthUncorrectedRegion(i);

			if (len_fw > len_bw) len_fw = -1;
			else if (len_fw < len_bw) len_bw = -1;
		}
		
		if (len_fw >= len_bw){

			const pair<pair<size_t, size_t>, size_t> pos_read = moveIntoCIGAR(i, i + len_fw, cigar_fw_pos, prev_cigar_fw_pos, fw_s_pos, ref_pos_fw, cigar_fw, cigar_fw_len);

			if (pos_read.first.second > pos_read.first.first){

				ss << fw_s->getSequence().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);

				if (haveQual) sq << fw_s->getQuality().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);
			}

			i = pos_read.second;
		}
		else {

			const pair<pair<size_t, size_t>, size_t> pos_read = moveIntoCIGAR(i, i + len_bw, cigar_bw_pos, prev_cigar_bw_pos, bw_s_pos, ref_pos_bw, cigar_bw, cigar_bw_len);

			if (pos_read.first.second > pos_read.first.first){

				ss << bw_s->getSequence().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);

				if (haveQual) sq << bw_s->getQuality().substr(pos_read.first.first, pos_read.first.second - pos_read.first.first);
			}

			i = pos_read.second;
		}
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

/*void fixAmbiguity(string& query, string& quality, const vector<pair<size_t, char>>& v_ambiguity, const char* ref_seq, const char* ref_qual, const size_t ref_len, const size_t min_qual, const bool force_fix){

	if (!v_ambiguity.empty()) {

		const size_t quality_len = quality.length();
		const size_t query_len = query.length();

		const bool q_hasQual = (quality_len != 0);
		const bool r_hasQual = (ref_qual != nullptr);

		const size_t min_qs_char = 33 + min_qual;

		string query_tmp = query;

		for (const auto p : v_ambiguity) query_tmp[p.first] = p.second;

		const EdlibAlignConfig config_ambiguous = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
		EdlibAlignResult align = edlibAlign(query.c_str(), query_len, ref_seq, ref_len, config_ambiguous);

		char* cigar = edlibAlignmentToCigar(align.alignment, align.alignmentLength, EDLIB_CIGAR_STANDARD);

		const size_t cigar_len = strlen(cigar);

		size_t cigar_pos = 0;
		size_t prev_cigar_pos = 0;
		size_t target_pos = align.startLocations[0];
		size_t query_pos = 0;

		while (cigar_pos != cigar_len){

			if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

				if (cigar[cigar_pos] == 'M') { //match

					const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

					for (size_t query_pos_tmp = query_pos, target_pos_tmp = target_pos; query_pos_tmp < query_pos + cigar_l; ++query_pos_tmp, ++target_pos_tmp){

						if (!isDNA(query_tmp[query_pos_tmp])){

							if (isDNA(ref_seq[target_pos_tmp]) && (!r_hasQual || (ref_qual[target_pos_tmp] >= min_qs_char))){

								bool is_a_q = false, is_c_q = false, is_g_q = false, is_t_q = false;
								bool is_a_t = false, is_c_t = false, is_g_t = false, is_t_t = false;

								getAmbiguityRev(ref_seq[target_pos_tmp], is_a_t, is_c_t, is_g_t, is_t_t);
								getAmbiguityRev(query_tmp[query_pos_tmp], is_a_q, is_c_q, is_g_q, is_t_q);

								if ((is_a_q && is_a_t) || (is_c_q && is_c_t) || (is_g_q && is_g_t) || (is_t_q && is_t_t)) query_tmp[query_pos_tmp] = ref_seq[target_pos_tmp];
								//else if (q_hasQual) quality[query_pos_tmp] = getQual(0.0);
							}

							if (force_fix && !isDNA(query_tmp[query_pos_tmp])){

								query_tmp[query_pos_tmp] = query[query_pos_tmp];

								//if (q_hasQual) quality[query_pos_tmp] = getQual(0.0);
							}
						}
					}

					query_pos += cigar_l;
					target_pos += cigar_l;
				}
				else if ((cigar[cigar_pos] == 'I') || (cigar[cigar_pos] == 'S')) query_pos += atoi(&cigar[prev_cigar_pos]); //insertion or soft-clipping
				else if (cigar[cigar_pos] == 'D') target_pos += atoi(&cigar[prev_cigar_pos]);  //deletion

				prev_cigar_pos = cigar_pos + 1;
			}

			++cigar_pos;
		}

		free(cigar);

		edlibFreeAlignResult(align);

		if (force_fix){

			for (size_t i = 0; i < query_len; ++i){

				if (!isDNA(query_tmp[i])){

					query_tmp[i] = query[i];

					//if (q_hasQual) quality[i] = getQual(0.0);
				}
			}
		}

		query = move(query_tmp);
	}
}*/

void fixAmbiguity(	const CompactedDBG<UnitigData>& dbg, string& query, string& quality, const vector<pair<size_t, char>>& v_ambiguity, 
					const char* ref_seq, const char* ref_qual, const size_t ref_len, const size_t min_qual, const bool force_fix){

	if (!v_ambiguity.empty()) {

		const size_t quality_len = quality.length();
		const size_t query_len = query.length();

		const bool q_hasQual = (quality_len != 0);
		const bool r_hasQual = (ref_qual != nullptr);

		const size_t min_qs_char = 33 + min_qual;

		string query_tmp = query;

		unordered_map<size_t, char> m_ambiguity;

		for (const auto p : v_ambiguity){

			m_ambiguity.insert(p);
			query_tmp[p.first] = p.second;
		}

		const EdlibAlignConfig config_ambiguous = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
		EdlibAlignResult align = edlibAlign(query.c_str(), query_len, ref_seq, ref_len, config_ambiguous);

		char* cigar = edlibAlignmentToCigar(align.alignment, align.alignmentLength, EDLIB_CIGAR_STANDARD);

		const size_t cigar_len = strlen(cigar);

		size_t cigar_pos = 0;
		size_t prev_cigar_pos = 0;
		size_t target_pos = align.startLocations[0];
		size_t query_pos = 0;

		while (cigar_pos != cigar_len){

			if ((cigar[cigar_pos] < 0x30) || (cigar[cigar_pos] > 0x39)){ // If current char. is not a number

				if (cigar[cigar_pos] == 'M') { //match

					const size_t cigar_l = atoi(&cigar[prev_cigar_pos]);

					for (size_t query_pos_tmp = query_pos, target_pos_tmp = target_pos; query_pos_tmp < query_pos + cigar_l; ++query_pos_tmp, ++target_pos_tmp){

						if (!isDNA(query_tmp[query_pos_tmp]) && isDNA(ref_seq[target_pos_tmp]) && (!r_hasQual || (ref_qual[target_pos_tmp] >= min_qs_char))){

							bool is_a_q = false, is_c_q = false, is_g_q = false, is_t_q = false;
							bool is_a_t = false, is_c_t = false, is_g_t = false, is_t_t = false;

							getAmbiguityRev(ref_seq[target_pos_tmp], is_a_t, is_c_t, is_g_t, is_t_t);
							getAmbiguityRev(query_tmp[query_pos_tmp], is_a_q, is_c_q, is_g_q, is_t_q);

							if ((is_a_q && is_a_t) || (is_c_q && is_c_t) || (is_g_q && is_g_t) || (is_t_q && is_t_t)){

								unordered_map<size_t, char>::iterator it = m_ambiguity.find(query_pos_tmp);

								if (it != m_ambiguity.end()) it->second = ref_seq[target_pos_tmp];
							}
						}
					}

					query_pos += cigar_l;
					target_pos += cigar_l;
				}
				else if ((cigar[cigar_pos] == 'I') || (cigar[cigar_pos] == 'S')) query_pos += atoi(&cigar[prev_cigar_pos]); //insertion or soft-clipping
				else if (cigar[cigar_pos] == 'D') target_pos += atoi(&cigar[prev_cigar_pos]);  //deletion

				prev_cigar_pos = cigar_pos + 1;
			}

			++cigar_pos;
		}

		free(cigar);
		edlibFreeAlignResult(align);

		const size_t k = dbg.getK();

		struct compare_ambiguity {

			bool operator()(const pair<size_t, char>& a, const pair<size_t, char>& b) const {

				return ((a.first < b.first) || ((a.first == b.first) && (a.second < b.second)));
			}
		};

		set<pair<size_t, char>, compare_ambiguity> s_ambiguity;

		for (const auto& p : m_ambiguity){

			if (isDNA(p.second)){

				const size_t start_pos = (p.first < k) ? 0 : (p.first - k);
				const size_t len_pos = (query_len - start_pos < 2 * k - 1) ? (query_len - start_pos) : (2 * k - 1);
				const size_t amb_pos_q_sub = p.first - start_pos;

				const bool confident_snp = (p.second == query[p.first]);

				string q_sub = query.substr(start_pos, len_pos);

				q_sub[amb_pos_q_sub] = p.second;

				for (KmerIterator it_km(q_sub.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

					const pair<Kmer, int>& p_km = *it_km;
					const const_UnitigMap<UnitigData> um = dbg.findUnitig(q_sub.c_str(), p_km.second, q_sub.length());

					if (!um.isEmpty) {

						const_UnitigMap<UnitigData> um_tmp = um;

						um_tmp.dist = 0;
						um_tmp.len = um_tmp.size - k + 1;

						const string unitig_seq = um_tmp.mappedSequenceToString();
						const vector<pair<size_t, char>> v_amb = um_tmp.getData()->get_ambiguity_char(um_tmp);

						size_t amb_pos_unitig = (amb_pos_q_sub - p_km.second) + um.dist;

						if (!um.strand) amb_pos_unitig = um.size - amb_pos_unitig - 1;

						for (const auto p_amb : v_amb){

							int64_t pos = static_cast<int64_t>(p_amb.first);

							if (pos <= amb_pos_unitig) pos = static_cast<int64_t>(p.first) - static_cast<int64_t>(amb_pos_unitig - pos);
							else pos = static_cast<int64_t>(p.first) + static_cast<int64_t>(pos - amb_pos_unitig);

							if ((pos >= 0) && (pos < query_len) && (pos != p.first)) {

								const unordered_map<size_t, char>::const_iterator it = m_ambiguity.find(pos);

								if ((it != m_ambiguity.end()) && (!isDNA(it->second) || confident_snp)) s_ambiguity.insert(pair<size_t, char>(pos, unitig_seq[p_amb.first]));
							}
						}

	                    it_km += um.len - 1;
					}
				}
			}
		}

		vector<pair<size_t, char>> v(s_ambiguity.begin(), s_ambiguity.end());

		for (int64_t i = 0; i < v.size(); ++i){

			if (((i == 0) || (v[i].first != v[i-1].first)) && ((i == v.size()-1) || (v[i].first != v[i+1].first))){

				unordered_map<size_t, char>::iterator it = m_ambiguity.find(v[i].first);

				if (it != m_ambiguity.end()){

					//if (isDNA(it->second)) it->second = v[i].second;
					//else {

						bool is_a_q = false, is_c_q = false, is_g_q = false, is_t_q = false;
						bool is_a_t = false, is_c_t = false, is_g_t = false, is_t_t = false;

						getAmbiguityRev(it->second, is_a_t, is_c_t, is_g_t, is_t_t);
						getAmbiguityRev(v[i].second, is_a_q, is_c_q, is_g_q, is_t_q);

						if ((is_a_q && is_a_t) || (is_c_q && is_c_t) || (is_g_q && is_g_t) || (is_t_q && is_t_t)) it->second = v[i].second;
					//}
				}
			}
		}

		for (const auto& p : m_ambiguity) query_tmp[p.first] = (force_fix && !isDNA(p.second)) ? query[p.first] : p.second;

		query = move(query_tmp);
	}
}

pair<int, int> selectBestSubstringAlignment(const char* ref, const size_t ref_len, const vector<Path<UnitigData>>& candidates, const double cut_threshold_norm_edit){

	const size_t k = Kmer::k;

	double best_edit_dist;

	int best_cand_id = -1;
	int best_end_loc = -1;

	EdlibAlignConfig config(edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));

	Path<UnitigData>::PathOut p_cand_before = candidates[0].toStringVector(), p_cand_curr;

	const string& candidate = p_cand_before.toString();

	EdlibAlignResult align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

	best_edit_dist = static_cast<double>(align.editDistance) / candidate.length();
	best_cand_id = 0;
	best_end_loc = align.endLocations[0];

	for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);

	edlibFreeAlignResult(align);

	for (int i = 1; i < candidates.size(); ++i){

		p_cand_curr = candidates[i].toStringVector(p_cand_before);

		const string& candidate = p_cand_curr.toString();
		const size_t max_len_norm = candidate.length();

		config = edlibNewAlignConfig(best_edit_dist * max_len_norm + 1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0);
		align = edlibAlign(candidate.c_str(), candidate.length(), ref, ref_len, config);

		if (align.editDistance >= 0){

			if ((static_cast<double>(align.editDistance) / max_len_norm) < best_edit_dist){

				best_edit_dist = static_cast<double>(align.editDistance) / max_len_norm;
				best_cand_id = i;
				best_end_loc = align.endLocations[0];

				for (size_t j = 1; j < align.numLocations; ++j) best_end_loc = max(best_end_loc, align.endLocations[j]);
			}
			else if ((static_cast<double>(align.editDistance) / max_len_norm) == best_edit_dist){

				for (size_t j = 0; j < align.numLocations; ++j){

					if (align.endLocations[j] > best_end_loc){

						best_cand_id = i;
						best_end_loc = align.endLocations[j];
					}
				}
			}
		}

		edlibFreeAlignResult(align);

		p_cand_before = move(p_cand_curr);
	}

	if ((cut_threshold_norm_edit > 0.0) && (best_edit_dist > cut_threshold_norm_edit)) return {-1, -1};

	return {best_cand_id, best_end_loc};
}