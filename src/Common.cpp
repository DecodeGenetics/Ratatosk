#include "Common.hpp"

// returns maximal entropy score for random sequences which is log2|A| where A is the alphabet
// returns something close to 0 for highly repetitive sequences
double getEntropy(const char* s, const size_t len){

	double counts[4] = {0, 0, 0, 0};
	double logs[4];

	for (size_t i = 0; i != len; ++i){

        const char c = s[i] & 0xDF;

		if (isDNA(c)){

			const size_t x = (c & 4) >> 1;

			++counts[x + ((x ^ (c & 2)) >> 1)];
		}
	}

	counts[0] /= len;
	counts[1] /= len;
	counts[2] /= len;
	counts[3] /= len;

	logs[0] = (counts[0] == 0) ? 0 : log2(counts[0]);
	logs[1] = (counts[1] == 0) ? 0 : log2(counts[1]);
	logs[2] = (counts[2] == 0) ? 0 : log2(counts[2]);
	logs[3] = (counts[3] == 0) ? 0 : log2(counts[3]);

	return 0 - ((counts[0] * logs[0]) + (counts[1] * logs[1]) + (counts[2] * logs[2]) + (counts[3] * logs[3]));
}

size_t getMaxPaths(const double seq_entropy, const size_t max_len_path, const size_t k){

	const double entropy_factor = 3. - seq_entropy;
	const size_t v = rndup(max_len_path);

	return static_cast<double>(v * (__builtin_ctzll(v) + k)) * entropy_factor;
}

size_t getMaxBranch(const double seq_entropy, const size_t max_len_path, const size_t k){

	const double entropy_factor = 3. - seq_entropy;
	const double base_nb_nodes = static_cast<double>(__builtin_ctzll(rndup(max_len_path - k + 1)));

	return (base_nb_nodes + 1) * entropy_factor;
}

size_t getNumberSharedPairID(const PairID& a, const PairID& b, const size_t min_shared_ids) {

	size_t nb_shared = 0;

	const size_t max_a = a.maximum(), min_a = a.minimum();
	const size_t max_b = b.maximum(), min_b = b.minimum();

	if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

		const size_t a_card = a.cardinality();
		const size_t b_card = b.cardinality();

		const size_t log2_a = b_card * approximate_log2(a_card);
		const size_t log2_b = a_card * approximate_log2(b_card);

		const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

		if (min_a_b == (a_card + b_card)) {

			PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();
			PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

			while ((a_it_s != a_it_e) && (b_it_s != b_it_e) && (nb_shared < min_shared_ids)){

				const uint32_t val_a = *a_it_s;
				const uint32_t val_b = *b_it_s;

				if (val_a == val_b){

					++nb_shared;
					++a_it_s;
					++b_it_s;
				}
				else if (val_a < val_b) ++a_it_s;
				else ++b_it_s;
			}
		}
		else if (min_a_b == log2_a){

			PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

			while ((b_it_s != b_it_e) && (*b_it_s < min_a)) ++b_it_s;

			while ((b_it_s != b_it_e) && (*b_it_s <= max_a) && (nb_shared < min_shared_ids)){

				nb_shared += static_cast<size_t>(a.contains(*b_it_s));
				++b_it_s;
			}
		}
		else {

			PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();

			while ((a_it_s != a_it_e) && (*a_it_s < min_b)) ++a_it_s;

			while ((a_it_s != a_it_e) && (*a_it_s <= max_b) && (nb_shared < min_shared_ids)){

				nb_shared += static_cast<size_t>(b.contains(*a_it_s));
				++a_it_s;
			}
		}
	}

	return nb_shared;
}

size_t getNumberSharedPairID(const SharedPairID& a, const SharedPairID& b, const size_t min_shared_ids) {

	size_t nb_shared = 0;

	const pair<const PairID*, const PairID*> pa = a.getPairIDs();
	const pair<const PairID*, const PairID*> pb = b.getPairIDs();

	if ((pa.first != nullptr) && (pa.first == pb.first)) {

		nb_shared = pa.first->cardinality();

		if (nb_shared < min_shared_ids) nb_shared += getNumberSharedPairID(*(pa.second), *(pb.second), min_shared_ids - nb_shared);
	}
	else {

		if (pb.first != nullptr) nb_shared = getNumberSharedPairID(a, *pb.first, min_shared_ids);
		if ((pb.second != nullptr) && (nb_shared < min_shared_ids)) nb_shared += getNumberSharedPairID(a, *pb.second, min_shared_ids - nb_shared);
	}

	return nb_shared;
}

size_t getNumberSharedPairID(const SharedPairID& a, const PairID& b, const size_t min_shared_ids) {

	size_t nb_shared = 0;

	const pair<const PairID*, const PairID*> p_pid = a.getPairIDs();

	if (p_pid.first != nullptr) nb_shared = getNumberSharedPairID(*p_pid.first, b, min_shared_ids);
	if ((p_pid.second != nullptr) && (nb_shared < min_shared_ids)) nb_shared += getNumberSharedPairID(*p_pid.second, b, min_shared_ids - nb_shared);

	return nb_shared;
}

size_t getNumberSharedPairID(const PairID& a, const PairID& b) {

	size_t nb_shared = 0;

	if (!a.isEmpty() && !b.isEmpty()) {

		const size_t max_a = a.maximum(), min_a = a.minimum();
		const size_t max_b = b.maximum(), min_b = b.minimum();

		if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

			const size_t a_card = a.cardinality();
			const size_t b_card = b.cardinality();

			const size_t log2_a = b_card * approximate_log2(a_card);
			const size_t log2_b = a_card * approximate_log2(b_card);

			const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

			if (min_a_b == (a_card + b_card)) {

				PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();
				PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

				while ((a_it_s != a_it_e) && (b_it_s != b_it_e)){

					const uint32_t val_a = *a_it_s;
					const uint32_t val_b = *b_it_s;

					if (val_a == val_b){

						++nb_shared;
						++a_it_s;
						++b_it_s;
					}
					else if (val_a < val_b) ++a_it_s;
					else ++b_it_s;
				}
			}
			else if (min_a_b == log2_a){

				PairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

				while ((b_it_s != b_it_e) && (*b_it_s < min_a)) ++b_it_s;

				while ((b_it_s != b_it_e) && (*b_it_s <= max_a)){

					nb_shared += static_cast<size_t>(a.contains(*b_it_s));
					++b_it_s;
				}
			}
			else {

				PairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();

				while ((a_it_s != a_it_e) && (*a_it_s < min_b)) ++a_it_s;

				while ((a_it_s != a_it_e) && (*a_it_s <= max_b)){

					nb_shared += static_cast<size_t>(b.contains(*a_it_s));
					++a_it_s;
				}
			}
		}
	}

	return nb_shared;
}

size_t getNumberSharedPairID(const SharedPairID& a, const PairID& b) {

	size_t nb_shared = 0;

	const pair<const PairID*, const PairID*> p_pid = a.getPairIDs();

	if (p_pid.first != nullptr) nb_shared += getNumberSharedPairID(*p_pid.first, b);
	if (p_pid.second != nullptr) nb_shared += getNumberSharedPairID(*p_pid.second, b);

	return nb_shared;
}

size_t getNumberSharedPairID(const SharedPairID& a, const SharedPairID& b) {

	size_t nb_shared = 0;

	const pair<const PairID*, const PairID*> pa = a.getPairIDs();
	const pair<const PairID*, const PairID*> pb = b.getPairIDs();

	if ((pa.first != nullptr) && (pa.first == pb.first)) nb_shared = pa.first->cardinality() + getNumberSharedPairID(*(pa.second), *(pb.second));
	else {

		if (pb.first != nullptr) nb_shared += getNumberSharedPairID(a, *pb.first);
		if (pb.second != nullptr) nb_shared += getNumberSharedPairID(a, *pb.second);
	}

	return nb_shared;
}

PairID getSharedPairID(const SharedPairID& a, const SharedPairID& b, const size_t min_shared) {

	PairID pid;

	const pair<const PairID*, const PairID*> pa = a.getPairIDs();
	const pair<const PairID*, const PairID*> pb = b.getPairIDs();

	if ((pa.first != nullptr) && (pb.first != nullptr) && (pa.first == pb.first)) {

		pid = *(pa.first);

		if (pa.first->cardinality() < min_shared) pid |= getSharedPairID(*(pa.second), *(pb.second), min_shared - pid.cardinality());
	}
	else {

		const size_t a_card = a.cardinality();
		const size_t b_card = b.cardinality();

		if ((a_card >= min_shared) && (b_card >= min_shared)) {

			const size_t max_a = a.maximum(), min_a = a.minimum();
			const size_t max_b = b.maximum(), min_b = b.minimum();

			if ((min_a <= max_b) && (min_b <= max_a)) { // Check that range overlaps (both bitmaps must be non empty!)

				const size_t log2_a = b_card * approximate_log2(a_card);
				const size_t log2_b = a_card * approximate_log2(b_card);

				const size_t min_a_b = min(a_card + b_card, min(log2_a, log2_b));

				size_t nb_shared = 0;

				if (min_a_b == (a_card + b_card)) {

					SharedPairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();
					SharedPairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

					while ((a_it_s != a_it_e) && (b_it_s != b_it_e) && (nb_shared < min_shared)){

						const uint32_t val_a = *a_it_s;
						const uint32_t val_b = *b_it_s;

						if (val_a == val_b){

							pid.add(val_a);

							++nb_shared;
							++a_it_s;
							++b_it_s;
						}
						else if (val_a < val_b) ++a_it_s;
						else ++b_it_s;
					}
				}
				else if (min_a_b == log2_a){

					SharedPairID::const_iterator b_it_s = b.begin(), b_it_e = b.end();

					while ((b_it_s != b_it_e) && (nb_shared < min_shared)) {

						if (a.contains(*b_it_s)) {

							pid.add(*b_it_s);
							++nb_shared;
						}

						++b_it_s;
					}
				}
				else {

					SharedPairID::const_iterator a_it_s = a.begin(), a_it_e = a.end();

					while ((a_it_s != a_it_e)  && (nb_shared < min_shared)){

						if (b.contains(*a_it_s)) {

							pid.add(*a_it_s);
							++nb_shared;
						}

						++a_it_s;
					}
				}
			}
		}
	}

	if (pid.cardinality() < min_shared) pid.clear();

	return pid;
}

size_t approximate_log2(size_t v) {

	static const uint8_t tab64[64] = {
	    63,  0, 58,  1, 59, 47, 53,  2,
	    60, 39, 48, 27, 54, 33, 42,  3,
	    61, 51, 37, 40, 49, 18, 28, 20,
	    55, 30, 34, 11, 43, 14, 22,  4,
	    62, 57, 46, 52, 38, 26, 32, 41,
	    50, 36, 17, 19, 29, 10, 13, 21,
	    56, 45, 25, 31, 35, 16,  9, 12,
	    44, 24, 15,  8, 23,  7,  6,  5
	};

    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;

    return tab64[(static_cast<size_t>((v - (v >> 1)) * 0x07EDD5E59A4E28C2)) >> 58];
}

bool check_files(const string& fn, const bool check_files_format, const bool verbose) {

	vector<string> v_fn(1, fn);

	return check_files(v_fn, check_files_format, verbose);
}

bool check_files(vector<string>& v_fn, const bool check_files_format, const bool verbose) {

    vector<string> v_fn_tmp;

    char* buffer = new char[4096]();

    for (const auto& file : v_fn) {

        if (!check_file_exists(file)) {

            if (verbose) cerr << "Ratatosk::Ratatosk():: File " << file << " not found." << endl;
            
            return false;
        }
        else if (check_files_format) {

            const int format = FileParser::getFileFormat(file.c_str());

            if (format >= 0) v_fn_tmp.push_back(file); // File is FASTA/FASTQ/GFA
            else {

                FILE* fp = fopen(file.c_str(), "r");

                if (fp != NULL){

                    fclose(fp);

                    ifstream ifs_file_txt(file);
                    istream i_file_txt(ifs_file_txt.rdbuf());

                    size_t i = 0;

                    while (i_file_txt.getline(buffer, 4096).good()){

                        fp = fopen(buffer, "r");

                        if (fp == NULL) {

                            if (verbose) cerr << "Ratatosk::Ratatosk():: Could not open file at line " << i << " in file " << file << " for reading." << endl;

                            return false;
                        }
                        else {

                            fclose(fp);

                            v_fn_tmp.push_back(string(buffer));
                        }

                        ++i;
                    }

                    if (i_file_txt.fail() && (i == 0)) {

                        if (verbose) {

                        	cerr << "Ratatosk::Ratatosk():: File " << file << " is neither FASTA, FASTQ nor GFA." << endl;
                        	cerr << "If it is a list of files, it is either empty or has a line with >4096 characters." << endl;
                        }

                        return false;
                    }

                    ifs_file_txt.close();
                }
                else {

                    if (verbose) cerr << "Ratatosk::Ratatosk():: Could not open file " << file << " for reading." << endl;

                    return false;
                }
            }
        }
        else {

            FILE* fp = fopen(file.c_str(), "r");

            if (fp == NULL){

                if (verbose)  cerr << "Ratatosk::Ratatosk(): Could not open file " << file << " for reading." << endl;

                return false;
            }
            else {

                fclose(fp);

                v_fn_tmp.push_back(file);
            }
        }
    }

    v_fn = move(v_fn_tmp);

    delete[] buffer;

    return true;
}

PairID subsample(const SharedPairID& spid, const size_t nb_id_out) {

	PairID pid;

	if ((nb_id_out != 0) && !spid.isEmpty()) {

		vector<uint32_t> v_ids = spid.toVector();

		std::random_shuffle(v_ids.begin(), v_ids.end());

		const size_t sz = min(nb_id_out, v_ids.size());

		for (size_t i = 0; i < sz; ++i) pid.add(v_ids[i]);
	}

	return pid;
}

PairID subsample(const PairID& spid, const size_t nb_id_out) {

	PairID pid;

	if ((nb_id_out != 0) && !spid.isEmpty()) {

		vector<uint32_t> v_ids = spid.toVector();

		std::random_shuffle(v_ids.begin(), v_ids.end());

		const size_t sz = min(nb_id_out, v_ids.size());

		for (size_t i = 0; i < sz; ++i) pid.add(v_ids[i]);
	}

	return pid;
}