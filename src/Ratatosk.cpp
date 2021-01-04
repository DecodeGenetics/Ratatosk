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

    Correct_Opt opt;

    cout << endl << "Ratatosk " << RATATOSK_VERSION << endl << endl;

    cout << "Hybrid error correction of long reads using colored de Bruijn graphs" << endl << endl;

    cout << "Usage: Ratatosk [PARAMETERS]" << endl << endl;

    cout << "[PARAMETERS]:" << endl << endl;
    cout << "   > Mandatory with required argument:" << endl << endl <<
    "   -s, --in-short                  Input short read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input short read files to correct (one file per line)" << endl <<
    "   -l, --in-long                   Input long read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input long read files to correct (one file per line)" << endl <<
    "   -o, --out-long                  Output file prefix" << endl <<
    endl << "   > Optional with required argument:" << endl << endl <<
    "   -c, --cores                     Number of cores (default: " << opt.nb_threads << ")" << endl <<
    "   -t, --trim-split                Trim and split bases with quality score < t (default: no trim/split)" << endl <<
    "                                   Only sub-read with length >= " << opt.k << " are output if used" << endl <<
    "   -u, --in-unmapped-short         Input read file of the unmapped short reads (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input read files of the unmapped short reads (one file per line)" << endl <<
    "   -a, --in-accurate-long          Input high quality long read file (FASTA/FASTQ possibly gzipped)" << endl <<
    "                                   List of input high quality long read files (one file per line)" << endl <<
    "                                   (Those reads are NOT corrected but assist the correction of reads in input)" << endl <<
    endl << "   > Optional with no argument:" << endl << endl <<
    "   -v, --verbose                   Print information" << endl << endl;

    cout << "[ADVANCED PARAMETERS]:" << endl << endl;
    cout << "   > Optional with required argument:" << endl << endl <<
    "   -i, --insert-sz                 Insert size of the input paired-end short reads (default: " << opt.insert_sz << ")" << endl <<
    "   -k, --k1                        Length of short k-mers for 1st pass (default: " << opt.small_k << ")" << endl <<
    "   -K, --k2                        Length of long k-mers for 2nd pass (default: " << opt.k << ")" << endl <<
    "   -w, --max-len-weak1             Do not correct non-solid regions >= w bases during 1st pass (default: " << opt.max_len_weak_region1 << ")" << endl <<
    "   -W, --max-len-weak2             Do not correct non-solid regions >= w bases during 2nd pass (default: " << opt.max_len_weak_region2 << ")" << endl << endl;
    cout << "   > Optional with no argument:" << endl << endl <<
    "   -1, --1st-pass-only             Perform *only* the 1st correction pass (default: " << (opt.pass1_only ? "true" : "false") << ")" << endl <<
    "   -2, --2nd-pass-only             Perform *only* the 2nd correction pass (default: " << (opt.pass2_only ? "true" : "false") << ")" << endl << endl;

    cout << "[EXPERIMENTAL PARAMETERS]:" << endl << endl;
    cout << "   > Optional with required argument:" << endl << endl <<
    "   -p, --in-short-phase            Input short read phasing file" << endl <<
    "                                   List of input short read phasing files (one file per line)" << endl <<
    "   -P, --in-long-phase             Input long read phasing file" << endl <<
    "                                   List of input long read phasing files (one file per line)" << endl << endl;
}

int parse_ProgramOptions(int argc, char **argv, Correct_Opt& opt) {

    int option_index = 0, c;

    const char* opt_string = "s:l:o:c:t:u:a:p:P:i:k:K:w:W:12v";

    static struct option long_options[] = {

        {"in-short",                required_argument,  0, 's'},
        {"in-long",                 required_argument,  0, 'l'},
        {"out-long",                required_argument,  0, 'o'},
        {"cores",                   required_argument,  0, 'c'},
        {"trim-split",              required_argument,  0, 't'},
        {"in-unmapped-short",       required_argument,  0, 'u'},
        {"in-accurate-long",        required_argument,  0, 'a'},
        {"in-short-phase",          required_argument,  0, 'p'},
        {"in-long-phase",           required_argument,  0, 'P'},
        {"insert_sz",               required_argument,  0, 'i'},
        {"k1",                      required_argument,  0, 'k'},
        {"k2",                      required_argument,  0, 'K'},
        {"max-len-weak1",           required_argument,  0, 'w'},
        {"max-len-weak2",           required_argument,  0, 'W'},
        {"1st-pass-only",           no_argument,        0, '1'},
        {"2nd-pass-only",           no_argument,        0, '2'},
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
            case 'u':
                opt.filenames_short_all.push_back(optarg);
                break;
            case 'a':
                opt.filenames_helper_long_in.push_back(optarg);
                break;
            case 'p':
                opt.filenames_short_phase.push_back(optarg);
                break;
            case 'P':
                opt.filenames_long_phase.push_back(optarg);
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
            case 'i':
                opt.insert_sz = atoi(optarg);
                break;
            case 'k':
                opt.small_k = atoi(optarg);
                break;
            case 'K':
                opt.k = atoi(optarg);
                break;
            case 'w':
                opt.max_len_weak_region1 = atoi(optarg);
                break;
            case 'W':
                opt.max_len_weak_region2 = atoi(optarg);
                break;
            case '1':
                opt.pass1_only = true;
                break;
            case '2':
                opt.pass2_only = true;
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

                const int format = FileParser::getFileFormat(file.c_str());

                if (format == 0){

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp == NULL){

                        cerr << "Ratatosk::Ratatosk(): Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                    else {

                        fclose(fp);

                        files_tmp.push_back(file);
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

    if (opt.k < opt.small_k){

        cerr << "Ratatosk::Ratatosk(): Length of long k-mers for 2nd pass correction cannot be less than or equal to length of short k-mers for 1st pass correction (" << opt.k << " and " << opt.small_k << " given)." << endl;
        ret = false;  
    }

    if (opt.insert_sz == 0){

        cerr << "Ratatosk::Ratatosk(): Insert size of short reads cannot be less than or equal to 0" << endl;
        ret = false;  
    }

    if ((opt.max_len_weak_region1 == 0) || (opt.max_len_weak_region2 == 0)){

        cerr << "Ratatosk::Ratatosk(): Maximum length of a weak region to correct cannot be less than or equal to 0" << endl;
        ret = false;  
    }

    if (opt.pass1_only && opt.pass2_only){

        cerr << "Ratatosk::Ratatosk(): -1 and -2 are mutually exclusive (perform *only* one of the two correction passes). To perform both, remove -1 and -2 from your command line." << endl;
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
    if (!opt.filenames_short_phase.empty()) ret = ret && check_files(opt.filenames_short_phase);
    if (!opt.filenames_long_phase.empty()) ret = ret && check_files(opt.filenames_long_phase);

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

void writeCorrectedOutput(ostream& out, const string& name_str, const string& seq_str, const string& qual_str, const int k, const int trim = 0){

    const char header_1st_c = '@';

    if (trim == 0){

        out << header_1st_c << name_str << "\n" << seq_str << "\n";
        out << "+" << "\n" << qual_str << "\n";
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
                    out << "+" << "\n" << qual_str.substr(start_pos, len) << "\n";
                }

                start_pos = -1;
                len = -1;
            }
        }

        if (len >= k){

            out << header_1st_c << name_str << "/" << id_subread++ << "\n" << seq_str.substr(start_pos, len) << "\n";
            out << "+" << "\n" << qual_str.substr(start_pos, len) << "\n";
        }
    }
}

void search(const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const bool long_read_correct, const Roaring* partitions, const pair<HapReads, HapReads>& hap_reads) {

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

            uint64_t hap_id = 0xffffffffffffffffULL;

            std::transform(in_read.begin(), in_read.end(), in_read.begin(), ::toupper); // Read sequence in upper case characters

            string in_qual = (fp.getQualityScoreString() != nullptr) ? string(fp.getQualityScoreString()) : string();

            if (!long_read_correct && !in_qual.empty()) getStdQual(in_qual); // Set all quality score from 0 to 40

            if (!hap_reads.second.read2hap.empty()){

                const uint64_t in_name_h = XXH64(in_name.c_str(), in_name.length(), opt.h_seed);
                const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_read2hap = hap_reads.second.read2hap.find(in_name_h);

                if (it_read2hap != hap_reads.second.read2hap.end()) hap_id = it_read2hap->second;
            }

            const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(opt, dbg, in_read, in_qual, long_read_correct, hap_id);
            const pair<string, string> correction = correctSequence(dbg, opt, in_read, in_qual, p.first, p.second, long_read_correct, partitions, hap_id, hap_reads);

            writeCorrectedOutput(out, in_name, correction.first, correction.second, k, (long_read_correct ? opt.trim_qual : 0));

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

                    vector<string> v_in_read, v_in_qual, v_in_name;

                    string in_read;

                    size_t in_read_len = 0;

                    while (true) {

                        mutex_file_in.acquire();

                        if (stop){

                            mutex_file_in.release();
                            return;
                        }
                        else {

                            while (in_read_len < opt.buffer_sz) {

                                stop = !fp.read(in_read, file_id);

                                if (!stop){

                                    in_read_len += in_read.length();

                                    v_in_read.push_back(move(in_read));
                                    v_in_name.push_back(string(fp.getNameString()));
                                    v_in_qual.push_back((fp.getQualityScoreString() != nullptr) ? string(fp.getQualityScoreString()) : string());

                                    if (opt.verbose && (++nb_reads_proc % 1000 == 0)) cout << "Ratatosk::correct(): Processed " << nb_reads_proc << " reads " << endl;
                                }
                                else break;
                            }

                            mutex_file_in.release();

                            if (!v_in_read.empty()){

                                for (size_t i = 0; i < v_in_read.size(); ++i){

                                    uint64_t hap_id = 0xffffffffffffffffULL;

                                    std::transform(v_in_read[i].begin(), v_in_read[i].end(), v_in_read[i].begin(), ::toupper); // Read sequence in upper case characters

                                    if (!long_read_correct && !v_in_qual[i].empty()) getStdQual(v_in_qual[i]); // Set all quality score from 0 to 40

                                    if (!hap_reads.second.read2hap.empty()){ // Fetch the haploblock ID of that read, if it exists or the read has one

                                        const uint64_t in_name_h = XXH64(v_in_name[i].c_str(), v_in_name[i].length(), opt.h_seed);
                                        const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_read2hap = hap_reads.second.read2hap.find(in_name_h);

                                        if (it_read2hap != hap_reads.second.read2hap.end()) hap_id = it_read2hap->second;
                                    }

                                    const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(opt, dbg, v_in_read[i], v_in_qual[i], long_read_correct, hap_id);

                                    pair<string, string> correction = correctSequence(dbg, opt, v_in_read[i], v_in_qual[i], p.first, p.second, long_read_correct, partitions, hap_id, hap_reads);

                                    v_in_read[i] = move(correction.first);
                                    v_in_qual[i] = move(correction.second);
                                }

                                mutex_file_out.acquire();

                                for (size_t i = 0; i < v_in_read.size(); ++i) writeCorrectedOutput(out, v_in_name[i], v_in_read[i], v_in_qual[i], k, (long_read_correct ? opt.trim_qual : 0));

                                mutex_file_out.release();
                            }

                            in_read_len = 0;

                            v_in_read.clear();
                            v_in_qual.clear();
                            v_in_name.clear();

                            in_read.clear();
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

            Correct_Opt opt_pass1(opt);
            Correct_Opt opt_pass2(opt);

            bool dbg_empty = false;

            string filename_out_extra_sr, sr_graph_long;

            // Establish if some reads are missing
            if (!opt.filenames_short_all.empty()){

                filename_out_extra_sr = retrieveMissingReads(opt);

                if (filename_out_extra_sr.length() != 0){

                    FILE* f = fopen(filename_out_extra_sr.c_str(), "r");

                    if (f != NULL) {

                        fclose(f);

                        opt_pass1.filename_seq_in.push_back(filename_out_extra_sr);
                    }
                    else filename_out_extra_sr = "";
                }
            }

            if (opt.pass1_only || !opt.pass2_only) { // Correction pass 1

                CompactedDBG<UnitigData> dbg(opt_pass1.k);

                if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Building graph from short reads (1/2)." << endl;

                sr_graph_long = opt.filename_long_out + "_k" + std::to_string(opt.k);

                dbg.build(opt_pass1);
                dbg.write(sr_graph_long, opt_pass1.nb_threads, false, opt_pass1.verbose);

                dbg_empty = (dbg.length() == 0);
                sr_graph_long += ".fasta";

                if (!dbg_empty){

                    pair<HapReads, HapReads> hapPass1; // <short reads, long reads>

                    Roaring* partitions = nullptr;

                    opt_pass1.k = opt_pass1.small_k;

                    opt_pass1.filename_ref_in = vector<string>(1, sr_graph_long);
                    opt_pass1.filename_seq_in = vector<string>();

                    dbg = CompactedDBG<UnitigData>(opt_pass1.k);

                    dbg.build(opt_pass1);

                    opt_pass1.filename_ref_in = vector<string>();
                    opt_pass1.filename_seq_in = opt.filename_seq_in;

                    if (filename_out_extra_sr.length() != 0) opt_pass1.filename_seq_in.push_back(filename_out_extra_sr);

                    opt_pass1.filename_long_out += "2.fastq";

                    if (!opt.filenames_long_phase.empty()){

                        if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding phasing to graph (1/2)." << endl;

                        addPhasing(dbg, opt, hapPass1.first, false, false);

                        hapPass1.second.hapBlock2id = hapPass1.first.hapBlock2id;
                        hapPass1.second.hapType2id = hapPass1.first.hapType2id;

                        addPhasing(dbg, opt, hapPass1.second, true, false);
                    }

                    if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding colors and coverage to graph (1/2)." << endl;

                    const vector<Kmer> v_km_centroids = addCoverage(dbg, opt_pass1, hapPass1.first, false, true);

                    if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding SNPs candidates to graph (1/2)." << endl;

                    detectSNPs(dbg, opt_pass1, partitions);

                    if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Correcting long reads (1/2)." << endl;

                    search(dbg, opt_pass1, false, partitions, hapPass1);

                    if (partitions != nullptr) delete[] partitions;
                }
                else if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Graph is empty, no correction can be done. Output uncorrected reads." << endl;
            }

            if (opt.pass2_only || !opt.pass1_only) {

                if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Building graph from short reads (2/2)." << endl;

                pair<HapReads, HapReads> hapPass2;

                CompactedDBG<UnitigData> dbg(opt_pass2.k);

                if (opt.pass2_only){

                    dbg.build(opt_pass2);

                    opt_pass2.filename_seq_in.clear();
                    opt_pass2.filename_seq_in.insert(opt_pass2.filename_seq_in.end(), opt_pass2.filenames_helper_long_in.begin(), opt_pass2.filenames_helper_long_in.end());
                    opt_pass2.filename_seq_in.insert(opt_pass2.filename_seq_in.end(), opt_pass2.filenames_long_in.begin(), opt_pass2.filenames_long_in.end());
                }
                else {

                    opt_pass2.filename_seq_in.clear();
                    opt_pass2.filename_seq_in.insert(opt_pass2.filename_seq_in.end(), opt_pass2.filenames_helper_long_in.begin(), opt_pass2.filenames_helper_long_in.end());
                    opt_pass2.filename_seq_in.push_back(opt_pass1.filename_long_out);

                    opt_pass2.filenames_long_in.clear();
                    opt_pass2.filenames_long_in.push_back(opt_pass1.filename_long_out);

                    dbg.read(sr_graph_long, opt_pass2.nb_threads, opt_pass2.verbose);

                    // Remove long k-mer graph and unmapped reads temporary files from disk
                    if ((sr_graph_long.length() != 0) && (remove(sr_graph_long.c_str()) != 0)) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;
                    if ((filename_out_extra_sr.length() != 0) && (remove(filename_out_extra_sr.c_str()) != 0)) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;

                    sr_graph_long = ""; // Graph was read and deleted
                    filename_out_extra_sr = ""; // Unmapped reads were read and deleted
                }

                opt_pass2.filename_long_out += ".fastq";

                // Pass 2
                if (!dbg_empty || opt.pass2_only){

                    if (dbg.length() != 0) {

                        if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Adding colors and coverage to graph (2/2)." << endl;

                        addCoverage(dbg, opt_pass2, hapPass2.second, true, true);

                        if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Adding SNPs candidates to graph (2/2)." << endl;

                        detectSNPs(dbg, opt_pass2, nullptr);

                        if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Correcting long reads (2/2)." << endl;

                        search(dbg, opt_pass2, true, nullptr, hapPass2);
                    }
                    else {

                        if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Graph is empty, no correction can be done. Output 1st correction pass reads. (2/2)" << endl;

                        copyFile(opt_pass2.filename_long_out, opt_pass2.filenames_long_in); // Graph is empty, LR cannot be corrected
                    }
                }
                else copyFile(opt_pass2.filename_long_out, opt.filenames_long_in); // Graph is empty, LR cannot be corrected
            }

            if ((sr_graph_long.length() != 0) && (remove(sr_graph_long.c_str()) != 0)) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl; // Clean up disk data
            if (!opt.pass1_only && !opt.pass2_only && !dbg_empty && remove(opt_pass1.filename_long_out.c_str()) != 0) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;
        }
    }
}