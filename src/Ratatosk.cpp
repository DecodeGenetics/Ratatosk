#include <iostream>
#include <sstream>

#include "CompactedDBG.hpp"

#include "Common.hpp"
#include "Correction.hpp"
#include "Graph.hpp"
#include "UnitigData.hpp"

using namespace std;

void PrintVersion() {

    cout << RATATOSK_VERSION << endl;
}

void PrintCite() {

    cout << "The paper describing this software has not been published." << endl;
}

void PrintUsage() {

    cout << endl << "Ratatosk " << RATATOSK_VERSION << endl << endl;

    cout << "Phased hybrid error correction of long reads using colored de Bruijn graphs" << endl << endl;

    cout << "Usage: Ratatosk [PARAMETERS]" << endl << endl;

    cout << "[PARAMETERS]:" << endl << endl;
    cout << "   > Mandatory with required argument:" << endl << endl <<
    "   -s, --in-short                  Input short read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input short read files to correct (one file per line)" << endl <<
    "   -l, --in-long                   Input long read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input long read files to correct (one file per line)" << endl <<
    "   -o, --out-long                  Output corrected long read file" << endl <<
    endl << "   > Optional with required argument:" << endl << endl <<
    "   -c, --cores                     Number of cores (default: 1)" << endl <<
    "   -q, --quality                   Output Quality Scores: corrected bases get QS >= t (default: t=0, no output)" << endl <<
    "   -t, --trimming                  Trim bases with quality score < t (default: t=0, no trimming)" << endl <<
    "                                   Only sub-read with length >= 63 are output if t > 0" << endl <<
    "   -d, --in-full-dataset-short     Input read file of the full short read dataset (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input read files of the full short read dataset (one file per line)" << endl <<
    "   -a, --in-accurate-long          Input high quality long read file (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input high quality long read files (one file per line)" << endl <<
    "                                   (Those reads are NOT corrected but assist the correction of reads in input)." << endl <<
    endl << "   > Optional with no argument:" << endl << endl <<
    "   -v, --verbose                   Print information messages during execution" << endl << endl;
}

int parse_ProgramOptions(int argc, char **argv, Correct_Opt& opt) {

    int option_index = 0, c;

    const char* opt_string = "s:l:o:c:q:t:d:a:v";

    static struct option long_options[] = {

        {"in-short",                required_argument,  0, 's'},
        {"in-long",                 required_argument,  0, 'l'},
        {"out-long",                required_argument,  0, 'o'},
        {"cores",                   required_argument,  0, 'c'},
        {"quality",                 required_argument,  0, 'q'},
        {"trim",                    required_argument,  0, 't'},
        {"in-full-dataset-short",   required_argument,  0, 'd'},
        {"in-accurate-long",        required_argument,  0, 'a'},
        {"verbose",                 no_argument,        0, 'v'},
        {0,                         0,                  0,  0 }
    };

    if (strcmp(argv[1], "--version") == 0) return 1; // Print version
    if (strcmp(argv[1], "--help") == 0) return 2; // print help

    while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

        switch (c) {

            case 's':
                opt.filename_seq_in.push_back(optarg);
                break;
            case 'l':
                opt.filenames_long_in.push_back(optarg);
                break;
            case 'd':
                opt.filenames_short_all.push_back(optarg);
                break;
            case 'a':
                opt.filenames_helper_long_in.push_back(optarg);
                break;
            case 'o':
                opt.filename_long_out = optarg;
                break;
            case 'c':
                opt.nb_threads = atoi(optarg);
                break;
            case 't':
                opt.trim_qual = atoi(optarg);
                break;
            case 'q':
                opt.out_qual = atoi(optarg);
                break;
            case 'v':
                opt.verbose = true;
                break;
            default: break;
        }
    }

    return 0;
}

bool check_ProgramOptions(Correct_Opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [](vector<string>& v_files) -> bool {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        bool ret = true;

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "Ratatosk::Ratatosk(): File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const string s_ext = file.substr(file.find_last_of(".") + 1);

                if ((s_ext == "txt")){

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != NULL){

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        while (i_file_txt.getline(buffer, 4096)){

                            fp = fopen(buffer, "r");

                            if (fp == NULL) {

                                cerr << "Ratatosk::Ratatosk(): Could not open file " << buffer << " for reading." << endl;
                                ret = false;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "Ratatosk::Ratatosk(): Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
                else files_tmp.push_back(file);
            }
        }

        delete[] buffer;

        v_files = move(files_tmp);

        return ret;
    };

    if (opt.nb_threads <= 0){

        cerr << "Ratatosk::Ratatosk(): Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Ratatosk::Ratatosk(): Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if ((opt.trim_qual < 0) || (opt.trim_qual > 40)){

        cerr << "Ratatosk::Ratatosk(): Quality score trimming threshold cannot be less than 0 or more than 40 (" << opt.trim_qual << " given)." << endl;
        ret = false;  
    }

    if ((opt.out_qual < 0) || (opt.out_qual > 40)){

        cerr << "Ratatosk::Ratatosk(): Minimum quality score for corrected bases cannot be less than 1 or more than 40 (" << opt.out_qual << " given)." << endl;
        ret = false;  
    }

    if (opt.filename_seq_in.empty()) {

        cerr << "Ratatosk::Ratatosk(): Missing input short read files." << endl;
        ret = false;
    }
    else ret = ret && check_files(opt.filename_seq_in);

    if (opt.filenames_long_in.empty()) {

        cerr << "Ratatosk::Ratatosk(): Missing input long read files." << endl;
        ret = false;
    }
    else ret = ret && check_files(opt.filenames_long_in);

    if (!opt.filenames_helper_long_in.empty()) ret = ret && check_files(opt.filenames_helper_long_in);
    if (!opt.filenames_short_all.empty()) ret = ret && check_files(opt.filenames_short_all);

    if (opt.filename_long_out.length() == 0){

        cerr << "Ratatosk::Ratatosk(): Missing output long read file." << endl;
        ret = false;
    }
    else {

        FILE* fp = fopen(opt.filename_long_out.c_str(), "w");

        if (fp == NULL) {

            cerr << "Ratatosk::Ratatosk(): Could not open output long read file for writing: " << opt.filename_long_out << "." << endl;
            ret = false;
        }
        else {

            fclose(fp);

            if (remove(opt.filename_long_out.c_str()) != 0) cerr << "Ratatosk::Ratatosk(): Could not remove temporary file " << opt.filename_long_out << "." << endl;
        }
    }

    return ret;
}

void writeCorrectedOutput(ostream& out, const string& name_str, const string& seq_str, const string& qual_str, const int k, const bool out_qual, const int trim = 0){

    const char header_1st_c = out_qual ? '@' : '>';

    if (trim == 0){

        out << header_1st_c << name_str << "\n" << seq_str << "\n";

        if (out_qual) out << "+" << "\n" << qual_str << "\n";
    }
    else {

        const char c_min = static_cast<char>(trim + 33);
        const int64_t l_qual = qual_str.length();

        int64_t start_pos = -1;
        int64_t len = -1;
        int64_t id_subread = 1;

        for (int64_t pos = 0; pos < l_qual; ++pos){

            if (qual_str[pos] >= c_min){ // Quality score of that base is good enough

                if (start_pos == -1) {

                    start_pos = pos;
                    len = 0;
                }

                ++len;
            }
            else {

                if (len >= k){

                    out << header_1st_c << name_str << "/" << id_subread++ << "\n" << seq_str.substr(start_pos, len) << "\n";

                    if (out_qual) out << "+" << "\n" << qual_str.substr(start_pos, len) << "\n";
                }

                start_pos = -1;
                len = -1;
            }
        }

        if (len >= k){

            out << header_1st_c << name_str << "/" << id_subread++ << "\n" << seq_str.substr(start_pos, len) << "\n";

            if (out_qual) out << "+" << "\n" << qual_str.substr(start_pos, len) << "\n";
        }
    }
}

void search(const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const bool long_read_correct, const Roaring* all_partitions) {

	const size_t k = dbg.getK();

    ofstream outfile;
    ostream out(0);

    FileParser fp(opt.filenames_long_in);

    size_t file_id;

    outfile.open(opt.filename_long_out.c_str());
    out.rdbuf(outfile.rdbuf());

    if (opt.nb_threads == 1){

    	string in_read;

	    while (fp.read(in_read, file_id)){ // For every read (a query)

            const string in_name = string(fp.getNameString());

            string in_qual = (fp.getQualityScoreString() != nullptr) ? string(fp.getQualityScoreString()) : string();

            if (!long_read_correct && !in_qual.empty()) getStdQual(in_qual); // Set all quality score from 0 to 40

	        const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(opt, dbg, in_read, long_read_correct);
	        const pair<string, string> correction = correctSequence(dbg, opt, in_read, in_qual, p.first, p.second, long_read_correct, all_partitions);

            if (long_read_correct) writeCorrectedOutput(out, in_name, correction.first, correction.second, k, static_cast<bool>(opt.out_qual), opt.trim_qual);
            else writeCorrectedOutput(out, in_name, correction.first, correction.second, k, static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual), 0);

	        if (opt.verbose){

		    	cout << fp.getNameString() << endl;
		    	cout << "Length before: " << in_read.length() << endl;
	        	cout << "Length after: " << correction.first.length() << endl;
	    	}
	    }
	}
	else {

        bool stop = false;

        size_t nb_reads_proc = 0;

        vector<thread> workers; // need to keep track of threads so we can join them

        SpinLock mutex_file_in;
        SpinLock mutex_file_out;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                	string in_read;
                    string in_qual;
                	string in_name;

                    while (true) {

                    	mutex_file_in.acquire();

                        if (stop){

                            mutex_file_in.release();
                        	return;
                        }
                        else {

                        	stop = !fp.read(in_read, file_id);

	                        if (!stop){

	                        	in_name = string(fp.getNameString());
                                in_qual = (fp.getQualityScoreString() != nullptr) ? string(fp.getQualityScoreString()) : string();

	                        	if (opt.verbose && (++nb_reads_proc % 1000 == 0)) cout << "Processed " << nb_reads_proc << " reads " << endl;

	                        	mutex_file_in.release();

                                if (!long_read_correct && !in_qual.empty()) getStdQual(in_qual); // Set all quality score from 0 to 40

						        const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(opt, dbg, in_read, long_read_correct/*, max_id_sr*/);
						        const pair<string, string> correction = correctSequence(dbg, opt, in_read, in_qual, p.first, p.second, long_read_correct, all_partitions);

	                            mutex_file_out.acquire();

                                if (long_read_correct) writeCorrectedOutput(out, in_name, correction.first, correction.second, k, static_cast<bool>(opt.out_qual), opt.trim_qual);
                                else writeCorrectedOutput(out, in_name, correction.first, correction.second, k, static_cast<bool>(opt.out_qual) || static_cast<bool>(opt.trim_qual), 0);

	                            mutex_file_out.release();

		                        in_read.clear();
		                        in_name.clear();
	                        }
							else {

								mutex_file_in.release();
								return;
							}
                    	}
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
	}

    outfile.close();
    fp.close();
}

int main(int argc, char *argv[]) {

    if (argc < 2) PrintUsage();
    else {

        Correct_Opt opt;

        const int parse = parse_ProgramOptions(argc, argv, opt);

        if (parse == 1) PrintVersion();
        else if (parse == 2) PrintUsage();
        else {

            if (!check_ProgramOptions(opt)) return 0;

            bool extra_sr = false;
            bool dbg_empty = false;

            string filename_out_extra_sr;
            string sr_graph = opt.filename_long_out + "_sr";

            // Establish if some reads are missing
            if (!opt.filenames_short_all.empty()){

                filename_out_extra_sr = retrieveMissingReads(opt);

                FILE* f = fopen(filename_out_extra_sr.c_str(), "r");

                if (f != NULL) {

                    fclose(f);

                    extra_sr = true;

                    opt.filename_seq_in.push_back(filename_out_extra_sr);
                }
            }

            Correct_Opt opt_pass1(opt);

            // Pass 1
            {
                CompactedDBG<UnitigData> dbg(opt.k);

                if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Building graph from short reads (1/2)." << endl;

                dbg.build(opt_pass1);
                //dbg.simplify(opt_pass1.deleteIsolated, opt_pass1.clipTips, opt_pass1.verbose);
                dbg.write(sr_graph, opt_pass1.nb_threads, false, opt_pass1.verbose);

                dbg_empty = (dbg.length() == 0);
                sr_graph += ".fasta";

                if (!dbg_empty){

                    Roaring* partitions = nullptr;

                    opt_pass1.k = opt_pass1.small_k;

                    opt_pass1.filename_ref_in = vector<string>(1, sr_graph);
                    opt_pass1.filename_seq_in = vector<string>();

                    dbg = CompactedDBG<UnitigData>(opt_pass1.k);

                    dbg.build(opt_pass1);

                    opt_pass1.filename_ref_in = vector<string>();
                    opt_pass1.filename_seq_in = opt.filename_seq_in;

                    opt_pass1.filename_long_out += string((opt_pass1.out_qual != 0) ? "2.fastq" : "2.fasta");

                    if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding coverage to vertices and edges (1/2)." << endl;

                    {
                        const vector<Kmer> v_km_centroids = addCoverage(dbg, opt_pass1, false, true);

                        partitions = createPartitions(dbg, v_km_centroids, dbg.nbKmers() / opt_pass1.nb_partitions, opt_pass1.nb_threads, opt_pass1.verbose);
                    }

                    if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Detecting SNPs candidates (1/2)." << endl;

                    detectSNPs(dbg, opt_pass1, partitions);

                    if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Correcting long reads (1/2)." << endl;

                    search(dbg, opt_pass1, false, partitions);

                    if (partitions != nullptr) delete[] partitions;
                }
                else if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Graph is empty, no correction can be done. Output uncorrected reads." << endl;
            }

            //const int sr_graph_format = FileParser::getFileFormat(sr_graph.c_str());

            Correct_Opt opt_pass2(opt);

            //opt_pass2.filename_ref_in.clear();

            //if (sr_graph_format != -1) opt_pass2.filename_ref_in.push_back(sr_graph);

            opt_pass2.filename_seq_in.clear();
            opt_pass2.filename_seq_in.insert(opt_pass2.filename_seq_in.end(), opt_pass2.filenames_helper_long_in.begin(), opt_pass2.filenames_helper_long_in.end());

            opt_pass2.filenames_long_in.clear();
            opt_pass2.filenames_long_in.push_back(opt_pass1.filename_long_out);

            opt_pass2.filename_long_out += string((opt_pass2.out_qual != 0) ? ".fastq" : ".fasta");

            opt_pass2.min_cov_vertices = 1;
            opt_pass2.min_cov_edges = 1;

            if (extra_sr && (remove(filename_out_extra_sr.c_str()) != 0)) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;

            // Pass 2
            if (!dbg_empty){

                if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Building graph from short reads (2/2)." << endl;

                CompactedDBG<UnitigData> dbg(opt_pass2.k);

                //if (opt_pass2.filenames_helper_long_in.empty()) dbg.read(sr_graph, opt_pass2.nb_threads, opt_pass2.verbose);
                //else dbg.build(opt_pass2);
                dbg.read(sr_graph, opt_pass2.nb_threads, opt_pass2.verbose);

                if (dbg.length() != 0) {

                    opt_pass2.filename_seq_in.push_back(opt_pass1.filename_long_out);

                    if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Adding coverage to vertices and edges (2/2)." << endl;

                    addCoverage(dbg, opt_pass2, true, true);

                    if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Detecting SNPs candidates (2/2)." << endl;

                    detectSNPs(dbg, opt_pass2, nullptr);

                    if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Correcting long reads (2/2)." << endl;

                    search(dbg, opt_pass2, true, nullptr);
                }
                else {

                    if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Graph is empty, no correction can be done. Output 1st correction pass reads. (2/2)" << endl;

                    copyFile(opt_pass2.filename_long_out, opt_pass2.filenames_long_in);
                }
            }
            else copyFile(opt_pass2.filename_long_out, opt.filenames_long_in); // Graph is empty, LR cannot be corrected

            // Clean up data on disk
            if (remove(sr_graph.c_str()) != 0) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;
            if (!dbg_empty && remove(opt_pass1.filename_long_out.c_str()) != 0) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;
        }
    }
}