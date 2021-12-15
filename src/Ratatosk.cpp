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

void PrintCitation() {

    cout << endl << "Holley, G. et al. Ratatosk: hybrid error correction of long reads enables accurate variant calling and assembly. Genome Biol. 22, 28 (2021)." << endl << endl;
}

void PrintUsage(const Correct_Opt& opt) {

    // Header, valid for all subcommands
    {
        cout << endl << "Ratatosk " << RATATOSK_VERSION << endl << endl;

        cout << "Hybrid error correction of long reads using colored de Bruijn graphs" << endl << endl;

        cout << "Usage: Ratatosk [COMMAND] [PARAMETERS]" << endl <<
        "Usage: Ratatosk --help" << endl <<
        "Usage: Ratatosk --version" << endl <<
        "Usage: Ratatosk --cite" << endl << endl;
    }

    if (opt.index == opt.correct){

        cout << "[COMMAND]:" << endl << endl;
        cout << "   correct                         Correct long reads with short reads" << endl;
        cout << "   index                           Prepare a Ratatosk index (advanced)" << endl << endl;

        cout << "Use \"Ratatosk [COMMAND] --help\" to get a specific command help" << endl << endl;
    }
    else if (opt.correct) {

        cout << "[COMMAND]: correct" << endl << endl;
        cout << "[PARAMETERS]:" << endl << endl;
        cout << "   > Mandatory with required argument:" << endl << endl <<
        "   -s, --in-short                  Input short read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input short read files to correct (one file per line)" << endl <<
        "   -l, --in-long                   Input long read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input long read files to correct (one file per line)" << endl <<
        "   -o, --out-long                  Output corrected long read file" << endl <<
        endl << "   > Optional with required argument:" << endl << endl <<
        "   -c, --cores                     Number of cores (default: " << opt.nb_threads << ")" << endl <<
        "   -S, --subsampling               Rate of short reads subsampling (default: " << ((opt.sampling_rate != 1.0) ? std::to_string(opt.sampling_rate) : string("Auto")) << ")" << endl <<
        "   -t, --trim-split                Trim and split bases with quality score < t (default: no trim/split)" << endl <<
        "                                   Only sub-read with length >= " << opt.k << " are output if used" << endl <<
        "   -u, --in-unmapped-short         Input read file of the unmapped short reads (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input read files of the unmapped short reads (one file per line)" << endl <<
        "   -a, --in-accurate-long          Input high quality long read file (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input high quality long read files (one file per line)" << endl <<
        "                                   (Those reads are NOT corrected but assist the correction of reads in input)" << endl <<
        "   -g, --in-graph                  Load graph file prepared with the index command" << endl <<
        "   -d, --in-unitig-data            Load unitig data file prepared with the index command" << endl <<
        endl << "   > Optional with no argument:" << endl << endl <<
        "   -v, --verbose                   Print information" << endl << endl;

        cout << "[ADVANCED PARAMETERS]:" << endl << endl;
        cout << "   > Optional with required argument:" << endl << endl <<
        "   -m, --min-conf-snp-corr         Minimum confidence threshold to correct a SNP (default: " << opt.min_confidence_snp_corr << ")" << endl <<
        "   -M, --min-conf-color2           Minimum confidence threshold to color vertices for 2nd pass (default: " << opt.min_confidence_2nd_pass << ")" << endl <<
        "   -C, --min-len-color2            Minimum length of a long read to color vertices for 2nd pass (default: " << opt.min_len_2nd_pass << ")" << endl <<
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
        //"   -r, --correction-rounds         Number of short read correction rounds (default: " << opt.nb_correction_rounds << ")" << endl <<
        "   -L, --in-long_raw               Input long read file from 1st pass (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input long read files to correct (one file per line)" << endl <<
        "   -p, --in-short-phase            Input short read phasing file (diploid only)" << endl <<
        "                                   List of input short read phasing files (one file per line)" << endl <<
        "   -P, --in-long-phase             Input long read phasing file (diploid only)" << endl <<
        "                                   List of input long read phasing files (one file per line)" << endl << endl;
        //cout << "   > Optional with no argument:" << endl << endl <<
        //"   -f, --force-correct-snp         Force unresolved SNP candidates correction (default: " << (opt.force_unres_snp_corr ? "true" : "false") << ")" << endl << endl;
    }
    else if (opt.index) {

        cout << "[COMMAND]: index" << endl << endl;
        cout << "[PARAMETERS]:" << endl << endl;
        cout << "   > Mandatory with required argument:" << endl << endl <<
        "   -s, --in-short                  Input short read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input short read files to correct (one file per line)" << endl <<
        "   -l, --in-long                   Input long read file to correct (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input long read files to correct (one file per line)" << endl <<
        "   -o, --out-long                  Prefix of the output index file" << endl <<
        endl << "   > Mandatory with no argument:" << endl << endl <<
        "   -1, --1st-pass-only             Prepare index for the 1st correction pass (default: " << (opt.pass1_only ? "true" : "false") << ")" << endl <<
        "   -2, --2nd-pass-only             Prepare index for the the 2nd correction pass (default: " << (opt.pass2_only ? "true" : "false") << ")" << endl <<
        endl << "   > Optional with required argument:" << endl << endl <<
        "   -c, --cores                     Number of cores (default: " << opt.nb_threads << ")" << endl <<
        "   -S, --subsampling               Rate of short reads subsampling (default: " << ((opt.sampling_rate != 1.0) ? std::to_string(opt.sampling_rate) : string("No subsampling")) << ")" << endl <<
        "   -u, --in-unmapped-short         Input read file of the unmapped short reads (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input read files of the unmapped short reads (one file per line)" << endl <<
        "   -a, --in-accurate-long          Input high quality long read file (FASTA/FASTQ possibly gzipped)" << endl <<
        "                                   List of input high quality long read files (one file per line)" << endl <<
        "                                   (Those reads are NOT corrected but assist the correction of reads in input)" << endl <<
        "   -g, --in-graph                  Load graph file prepared with the index command" << endl <<
        endl << "   > Optional with no argument:" << endl << endl <<
        "   -v, --verbose                   Print information" << endl << endl;

        cout << "[ADVANCED PARAMETERS]:" << endl << endl;
        cout << "   > Optional with required argument:" << endl << endl <<
        "   -M, --min-conf-color2           Minimum confidence threshold to color vertices for 2nd pass (default: " << opt.min_confidence_2nd_pass << ")" << endl <<
        "   -C, --min-len-color2            Minimum length of a long read to color vertices for 2nd pass (default: " << opt.min_len_2nd_pass << ")" << endl <<
        "   -i, --insert-sz                 Insert size of the input paired-end short reads (default: " << opt.insert_sz << ")" << endl <<
        "   -k, --k1                        Length of short k-mers for 1st pass (default: " << opt.small_k << ")" << endl <<
        "   -K, --k2                        Length of long k-mers for 2nd pass (default: " << opt.k << ")" << endl << endl;
    }
}

int parse_ProgramOptions(int argc, char **argv, Correct_Opt& opt) {

    int option_index = 0, c;

    const char* opt_string = "s:l:o:c:S:t:u:a:g:d:m:M:C:i:k:K:w:W:r:L:p:P:12fv";

    static struct option long_options[] = {

        {"in-short",                required_argument,  0, 's'},
        {"in-long",                 required_argument,  0, 'l'},
        {"out-long",                required_argument,  0, 'o'},
        {"cores",                   required_argument,  0, 'c'},
        {"sampling",                required_argument,  0, 'S'},
        {"trim-split",              required_argument,  0, 't'},
        {"in-unmapped-short",       required_argument,  0, 'u'},
        {"in-accurate-long",        required_argument,  0, 'a'},
        {"in-graph",                required_argument,  0, 'g'},
        {"in-unitig-data",          required_argument,  0, 'd'},
        {"min-conf-snp-corr",       required_argument,  0, 'm'},
        {"min-conf-color2",         required_argument,  0, 'M'},
        {"min-len-color2",          required_argument,  0, 'C'},
        {"insert-sz",               required_argument,  0, 'i'},
        {"k1",                      required_argument,  0, 'k'},
        {"k2",                      required_argument,  0, 'K'},
        {"max-len-weak1",           required_argument,  0, 'w'},
        {"max-len-weak2",           required_argument,  0, 'W'},
        {"correction-rounds",       required_argument,  0, 'r'},
        {"in-long-raw",             required_argument,  0, 'L'},
        {"in-short-phase",          required_argument,  0, 'p'},
        {"in-long-phase",           required_argument,  0, 'P'},
        {"1st-pass-only",           no_argument,        0, '1'},
        {"2nd-pass-only",           no_argument,        0, '2'},
        {"force-correct-snp",       no_argument,        0, 'f'},
        {"verbose",                 no_argument,        0, 'v'},
        {0,                         0,                  0,  0 }
    };

    if (argc <= 1) return 2; // No command or arguments -> print help

    if (strcmp(argv[1], "--version") == 0) return 1; // Print version
    if (strcmp(argv[1], "--help") == 0) return 2; // print help
    if (strcmp(argv[1], "--cite") == 0) return 3; // print citation

    if (strcmp(argv[1], "correct") == 0) opt.correct = true; // Correction command
    else if (strcmp(argv[1], "index") == 0) opt.index = true; // Indexing command
    else return 2; // No command known -> print help

    if ((argc == 2) || ((argc > 2) && (strcmp(argv[2], "--help") == 0))) return 2; // print command help

    while ((c = getopt_long(argc - 1, argv + 1, opt_string, long_options, &option_index)) != -1) {

        switch (c) {

            case 's':
                opt.filename_seq_in.push_back(string(optarg));
                break;
            case 'l':
                opt.filenames_long_in.push_back(string(optarg));
                break;
            case 'u':
                opt.filenames_short_all.push_back(string(optarg));
                break;
            case 'a':
                opt.filenames_helper_long_in.push_back(string(optarg));
                break;
            case 'p':
                opt.filenames_short_phase.push_back(string(optarg));
                break;
            case 'P':
                opt.filenames_long_phase.push_back(string(optarg));
                break;
            case 'o':
                opt.filename_long_out = string(optarg);
                break;
            case 'g':
                opt.fn_graph_in = string(optarg);
                break;
            case 'd':
                opt.fn_index_in = string(optarg);
                break;
            case 'c':
                opt.nb_threads = atoi(optarg);
                break;
            case 'S':
                opt.sampling_rate = atof(optarg);
                break;
            case 't':
                opt.trim_qual = atoi(optarg);
                break;
            case 'm':
                opt.min_confidence_snp_corr = atof(optarg);
                break;
            case 'M':
                opt.min_confidence_2nd_pass = atof(optarg);
                break;
            case 'C':
                opt.min_len_2nd_pass = atoi(optarg);
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
            case 'r':
                opt.nb_correction_rounds = atoi(optarg);
                break;
            case 'L':
                opt.filenames_long_raw.push_back(string(optarg));
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
            case 'f':
                opt.force_unres_snp_corr = true;
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

    const size_t max_threads = std::thread::hardware_concurrency();

    const bool hasGraph = (opt.fn_graph_in.length() != 0);
    const bool hasIndex = (opt.fn_index_in.length() != 0);

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

    if (opt.nb_correction_rounds < 1){

        cerr << "Ratatosk::Ratatosk(): At least one short read correction round is required." << endl;
        ret = false;  
    }

    if (opt.sampling_rate <= 0.0){

        cerr << "Ratatosk::Ratatosk(): Sampling rate must be greater than 0.0." << endl;
        ret = false;  
    }

    if (opt.sampling_rate > 1.0){

        cerr << "Ratatosk::Ratatosk(): Sampling rate must be lower or equal to 1.0." << endl;
        ret = false;  
    }

    if (opt.min_confidence_snp_corr < 0.0){

        cerr << "Ratatosk::Ratatosk(): Minimum confidence threshold to correct a SNP must be greater or equal to 0.0." << endl;
        ret = false;  
    }

    if (opt.min_confidence_snp_corr > 1.0){

        cerr << "Ratatosk::Ratatosk(): Minimum confidence threshold to correct a SNP must be lower or equal to 1.0." << endl;
        ret = false;  
    }

    if (opt.min_confidence_2nd_pass < 0.0){

        cerr << "Ratatosk::Ratatosk(): Minimum confidence threshold of a 1st pass corrected base to be used for coloring in 2nd pass must be greater or equal to 0.0." << endl;
        ret = false;  
    }

    if (opt.min_confidence_2nd_pass > 1.0){

        cerr << "Ratatosk::Ratatosk(): Minimum confidence threshold of a 1st pass corrected base to be used for coloring in 2nd pass must be lower or equal to 1.0." << endl;
        ret = false;  
    }

    if (opt.min_len_2nd_pass < 0){

        cerr << "Ratatosk::Ratatosk(): Minimum length of a 1st-pass-corrected read to be used for coloring in 2nd pass cannot be less than 0." << endl;
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
    else if (opt.index && !opt.pass1_only && !opt.pass2_only) {

        cerr << "Ratatosk::Ratatosk(): A correction pass (-1 or -2) must be selected to prepare an index for it." << endl;
        ret = false;
    }

    if (opt.correct) {

        if (hasGraph != hasIndex) {

            cerr << "Ratatosk::Ratatosk(): One of the input index files is missing (either the graph or the data)." << endl;
            ret = false;
        }
        else if (hasGraph && hasIndex) {

            if (!opt.pass1_only && !opt.pass2_only) {

                cerr << "Ratatosk::Ratatosk(): A correction pass (-1 or -2) must be selected to load its corresponding index." << endl;
                ret = false;
            }
            else {

                vector<string> v_fn_graph_in(1, opt.fn_graph_in);
                vector<string> v_fn_index_in(1, opt.fn_index_in);

                ret = ret && check_files(v_fn_graph_in, true, true);
                ret = ret && check_files(v_fn_index_in, false, true);
            }
        }
        else if (opt.filename_seq_in.empty()) {

            cerr << "Ratatosk::Ratatosk(): Missing input short read files." << endl;
            ret = false;
        }
        else ret = ret && check_files(opt.filename_seq_in, true, true);
    }
    else if (opt.index) {

        if (!opt.pass1_only && !opt.pass2_only) {

            cerr << "Ratatosk::Ratatosk(): A correction pass (-1 or -2) must be selected to compute its corresponding index." << endl;
            ret = false;
        }
        else if (hasGraph && opt.pass1_only) {

            cerr << "Ratatosk::Ratatosk(): The index for the 1st correction pass must be computed from the short reads, no input graph can be provided." << endl;
            ret = false;
        }
        else if (hasGraph && opt.pass2_only) {

            vector<string> v_fn_graph_in(1, opt.fn_graph_in);

            ret = ret && check_files(v_fn_graph_in, true, true);
        }
        else if (opt.filename_seq_in.empty()) {

            cerr << "Ratatosk::Ratatosk(): Missing input short read files." << endl;
            ret = false;
        }
        else ret = ret && check_files(opt.filename_seq_in, true, true);
    }

    if (opt.filenames_long_in.empty()) {

        cerr << "Ratatosk::Ratatosk(): Missing input long read files." << endl;
        ret = false;
    }
    else ret = ret && check_files(opt.filenames_long_in, true, true);

    if (!opt.filenames_helper_long_in.empty()) ret = ret && check_files(opt.filenames_helper_long_in, true, true);
    if (!opt.filenames_short_all.empty()) ret = ret && check_files(opt.filenames_short_all, true, true);
    if (!opt.filenames_long_raw.empty()) ret = ret && check_files(opt.filenames_long_raw, true, true);
    if (!opt.filenames_short_phase.empty()) ret = ret && check_files(opt.filenames_short_phase, true, true);
    if (!opt.filenames_long_phase.empty()) ret = ret && check_files(opt.filenames_long_phase, true, true);

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

size_t writeCorrectedOutput(ostream& out, const string& name_str, const string& seq_str, const string& qual_str, const int k, const int trim = 0){

    const char header_1st_c = '@';

    const streampos spos_b = out.tellp();

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

    return static_cast<size_t>(out.tellp() - spos_b);
}

void search(const CompactedDBG<UnitigData>& dbg, const Correct_Opt& opt, const bool long_read_correct, const Roaring* partitions, const pair<HapReads, HapReads>& hap_reads) {

    const size_t k = dbg.getK();
    const string fn_out = opt.filename_long_out + ".fastq";

    const size_t max_km_cov = max(getMaxKmerCoverage(dbg, opt.top_km_cov_ratio), opt.max_km_cov);

    ofstream outfile;
    ostream out(0);

    vector<pair<size_t, streampos>> v_block_w; // Blocks order and position if reordering needed

    FileParser fp(opt.filenames_long_in);

    size_t file_id;

    outfile.open(fn_out.c_str());
    out.rdbuf(outfile.rdbuf());

    if (long_read_correct && opt.verbose && opt.force_unres_snp_corr) cerr << "Ratatosk::search(): Force unresolved SNP correction is activated." << endl;

    if (opt.nb_threads == 1){

        string in_read;

        while (fp.read(in_read, file_id)){ // For every read (a query)

            const string in_name = string(fp.getNameString());
            const size_t in_len = in_read.length();

            uint64_t hap_id = 0xffffffffffffffffULL;

            unordered_map<Kmer, vector<const_UnitigMap<UnitigData>>, KmerHash> m_km_um;

            std::transform(in_read.begin(), in_read.end(), in_read.begin(), ::toupper); // Read sequence in upper case characters

            string in_qual = (fp.getQualityScoreString() != nullptr) ? string(fp.getQualityScoreString()) : string();

            if (!long_read_correct && !in_qual.empty()) getStdQual(in_qual); // Set all quality score from 0 to 40

            if (!hap_reads.second.read2hap.empty()){

                const uint64_t in_name_h = XXH64(in_name.c_str(), in_name.length(), opt.h_seed);
                const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_read2hap = hap_reads.second.read2hap.find(in_name_h);

                if (it_read2hap != hap_reads.second.read2hap.end()) hap_id = it_read2hap->second;
            }

            if (long_read_correct){

                if (opt.force_unres_snp_corr) in_read = fixSNPs(opt, dbg, in_read);

                //const auto start = std::chrono::high_resolution_clock::now();

                const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(opt, dbg, in_read, in_qual, long_read_correct, hap_id, m_km_um, false);

                //const auto middle = std::chrono::high_resolution_clock::now();

                pair<string, string> correction = correctSequence(dbg, opt, in_read, in_qual, p.first, p.second, long_read_correct, partitions, hap_id, hap_reads, max_km_cov);

                //const auto end = std::chrono::high_resolution_clock::now();

                //cout << std::chrono::duration_cast<std::chrono::milliseconds>(middle - start).count() << " " << std::chrono::duration_cast<std::chrono::milliseconds>(end - middle).count()
                //<< " " << (static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - middle).count()) / static_cast<double>(in_read.length())) << endl;

                in_read = move(correction.first);
                in_qual = move(correction.second);
            }
            else {

                const double step_min_score = 1.00 / static_cast<double>(opt.nb_correction_rounds);
                const double step_weak_region_len_factor = (opt.nb_correction_rounds == 1) ? 0.0 : ((opt.weak_region_len_factor - 0.10) / static_cast<double>(opt.nb_correction_rounds - 1));
                const size_t step_max_len_weak_region1 = opt.max_len_weak_region1 / opt.nb_correction_rounds;

                for (size_t j = 0; j < opt.nb_correction_rounds; ++j) {

                    Correct_Opt l_opt = opt;

                    l_opt.min_score = 1.00 - (j+1) * step_min_score;
                    l_opt.weak_region_len_factor = opt.weak_region_len_factor - (opt.nb_correction_rounds - j - 1) * step_weak_region_len_factor;
                    l_opt.max_len_weak_region1 = (j+1) * step_max_len_weak_region1;

                    const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(l_opt, dbg, in_read, in_qual, long_read_correct, hap_id, m_km_um, (j+1) != opt.nb_correction_rounds);

                    pair<string, string> correction = correctSequence(dbg, l_opt, in_read, in_qual, p.first, p.second, long_read_correct, partitions, hap_id, hap_reads, max_km_cov);

                    in_read = move(correction.first);
                    in_qual = move(correction.second);
                }
            }

            writeCorrectedOutput(out, in_name, in_read, in_qual, k, (long_read_correct ? opt.trim_qual : 0));

            if (opt.verbose) cout << fp.getNameString() << ", " << in_len << " bp raw, " << in_read.length() << " bp corrected" << endl;
        }
    }
    /*else {

        bool stop = false;

        size_t nb_reads_proc = 0;
        size_t ticket_dispenser = 0;

        vector<thread> workers; // need to keep track of threads so we can join them

        SpinLock mutex_file_in;
        SpinLock mutex_file_out;

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    vector<string> v_in_read, v_in_qual, v_in_name;

                    string in_read;

                    size_t in_read_len = 0;
                    size_t ticket_id = 0;

                    while (true) {

                        mutex_file_in.acquire();

                        if (stop){

                            mutex_file_in.release();
                            return;
                        }
                        else {

                            ticket_id = ticket_dispenser++;

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

                                    unordered_map<Kmer, vector<const_UnitigMap<UnitigData>>, KmerHash> m_km_um;

                                    std::transform(v_in_read[i].begin(), v_in_read[i].end(), v_in_read[i].begin(), ::toupper); // Read sequence in upper case characters

                                    if (!long_read_correct && !v_in_qual[i].empty()) getStdQual(v_in_qual[i]); // Set all quality score from 0 to 40

                                    if (!hap_reads.second.read2hap.empty()){ // Fetch the haploblock ID of that read, if it exists or the read has one

                                        const uint64_t in_name_h = XXH64(v_in_name[i].c_str(), v_in_name[i].length(), opt.h_seed);
                                        const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_read2hap = hap_reads.second.read2hap.find(in_name_h);

                                        if (it_read2hap != hap_reads.second.read2hap.end()) hap_id = it_read2hap->second;
                                    }

                                    if (long_read_correct){

                                        if (opt.force_unres_snp_corr) v_in_read[i] = fixSNPs(opt, dbg, v_in_read[i]);

                                        const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(opt, dbg, v_in_read[i], v_in_qual[i], long_read_correct, hap_id, m_km_um, false);

                                        pair<string, string> correction = correctSequence(dbg, opt, v_in_read[i], v_in_qual[i], p.first, p.second, long_read_correct, partitions, hap_id, hap_reads, max_km_cov);

                                        v_in_read[i] = move(correction.first);
                                        v_in_qual[i] = move(correction.second);

                                        //v_in_read[i] = phasing(dbg, opt, v_in_read[i], p.first);
                                    }
                                    else {

                                        const double step_min_score = 1.00 / static_cast<double>(opt.nb_correction_rounds);
                                        const double step_weak_region_len_factor = (opt.nb_correction_rounds == 1) ? 0.0 : ((opt.weak_region_len_factor - 0.10) / static_cast<double>(opt.nb_correction_rounds - 1));
                                        const size_t step_max_len_weak_region1 = opt.max_len_weak_region1 / opt.nb_correction_rounds;

                                        for (size_t j = 0; j < opt.nb_correction_rounds; ++j) {

                                            Correct_Opt l_opt = opt;

                                            l_opt.min_score = 1.00 - (j+1) * step_min_score;
                                            l_opt.weak_region_len_factor = opt.weak_region_len_factor - (opt.nb_correction_rounds - j - 1) * step_weak_region_len_factor;
                                            l_opt.max_len_weak_region1 = (j+1) * step_max_len_weak_region1;

                                            const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(l_opt, dbg, v_in_read[i], v_in_qual[i], long_read_correct, hap_id, m_km_um, (j+1) != opt.nb_correction_rounds);

                                            pair<string, string> correction = correctSequence(dbg, l_opt, v_in_read[i], v_in_qual[i], p.first, p.second, long_read_correct, partitions, hap_id, hap_reads, max_km_cov);

                                            v_in_read[i] = move(correction.first);
                                            v_in_qual[i] = move(correction.second);
                                        }
                                    }
                                }

                                mutex_file_out.acquire();

                                size_t block_sz = 0;

                                const streampos sp = out.tellp();

                                for (size_t i = 0; i < v_in_read.size(); ++i) block_sz += writeCorrectedOutput(out, v_in_name[i], v_in_read[i], v_in_qual[i], k, (long_read_correct ? opt.trim_qual : 0));

                                v_block_w.push_back({(ticket_id << 32) | block_sz, sp});

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
    }*/
    else {

        bool stop = false;

        size_t nb_reads_proc = 0;
        size_t ticket_dispenser = 0;
        size_t file_id_raw = 0;

        vector<thread> workers; // need to keep track of threads so we can join them

        SpinLock mutex_file_in;
        SpinLock mutex_file_out;

        FileParser fp_raw(long_read_correct ? opt.filenames_long_raw : opt.filenames_long_in);

        for (size_t t = 0; t < opt.nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    vector<string> v_in_read, v_in_name, v_in_qual;

                    string in_read;

                    vector<string> v_in_read_raw, v_in_name_raw;

                    string in_read_raw;

                    size_t in_read_len = 0;
                    size_t ticket_id = 0;

                    while (true) {

                        mutex_file_in.acquire();

                        if (stop){

                            mutex_file_in.release();
                            return;
                        }
                        else {

                            ticket_id = ticket_dispenser++;

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

                            if (long_read_correct) {

                                while (v_in_read_raw.size() < v_in_read.size()) {

                                    if (fp_raw.read(in_read_raw, file_id_raw)) {

                                        v_in_read_raw.push_back(move(in_read_raw));
                                        v_in_name_raw.push_back(string(fp_raw.getNameString()));
                                    }
                                    else break;
                                }

                                if (v_in_read_raw.size() != v_in_read.size()) {

                                    cerr << "Ratatosk::correct(): Corrected read file is not in the same order as input long read file. Abort." << endl;
                                    exit(1);
                                }
                                else {

                                    for (size_t i = 0; i < v_in_read.size(); ++i) {

                                        if (v_in_name[i].compare(1, v_in_name[i].length() - 1, v_in_name_raw[i], 1, v_in_name_raw[i].length() - 1) != 0) {

                                            cerr << "Ratatosk::correct(): Corrected read file is not in the same order as input long read file. Abort." << endl;
                                            exit(1);
                                        }
                                    }
                                }
                            }

                            mutex_file_in.release();

                            if (!v_in_read.empty()){

                                for (size_t i = 0; i < v_in_read.size(); ++i){

                                    uint64_t hap_id = 0xffffffffffffffffULL;

                                    unordered_map<Kmer, vector<const_UnitigMap<UnitigData>>, KmerHash> m_km_um;

                                    std::transform(v_in_read[i].begin(), v_in_read[i].end(), v_in_read[i].begin(), ::toupper); // Read sequence in upper case characters

                                    if (!long_read_correct && !v_in_qual[i].empty()) getStdQual(v_in_qual[i]); // Set all quality score from 0 to 40

                                    if (!hap_reads.second.read2hap.empty()){ // Fetch the haploblock ID of that read, if it exists or the read has one

                                        const uint64_t in_name_h = XXH64(v_in_name[i].c_str(), v_in_name[i].length(), opt.h_seed);
                                        const unordered_map<uint64_t, uint64_t, CustomHashUint64_t>::const_iterator it_read2hap = hap_reads.second.read2hap.find(in_name_h);

                                        if (it_read2hap != hap_reads.second.read2hap.end()) hap_id = it_read2hap->second;
                                    }

                                    if (long_read_correct){

                                        if (opt.force_unres_snp_corr) v_in_read[i] = fixSNPs(opt, dbg, v_in_read[i]);

                                        // ========== TEST ==========
                                        {
                                            pair<string, string> new_reads = phasing(dbg, opt, v_in_read_raw[i], v_in_read[i], v_in_qual[i]);

                                            v_in_read[i] = move(new_reads.first);
                                            v_in_qual[i] = move(new_reads.second);
                                        }
                                        // ========== TEST ==========

                                        const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(opt, dbg, v_in_read[i], v_in_qual[i], long_read_correct, hap_id, m_km_um, false);

                                        pair<string, string> correction = correctSequence(dbg, opt, v_in_read[i], v_in_qual[i], p.first, p.second, long_read_correct, partitions, hap_id, hap_reads, max_km_cov);

                                        v_in_read[i] = move(correction.first);
                                        v_in_qual[i] = move(correction.second);
                                    }
                                    else {

                                        const double step_min_score = 1.00 / static_cast<double>(opt.nb_correction_rounds);
                                        const double step_weak_region_len_factor = (opt.nb_correction_rounds == 1) ? 0.0 : ((opt.weak_region_len_factor - 0.10) / static_cast<double>(opt.nb_correction_rounds - 1));
                                        const size_t step_max_len_weak_region1 = opt.max_len_weak_region1 / opt.nb_correction_rounds;

                                        for (size_t j = 0; j < opt.nb_correction_rounds; ++j) {

                                            Correct_Opt l_opt = opt;

                                            l_opt.min_score = 1.00 - (j+1) * step_min_score;
                                            l_opt.weak_region_len_factor = opt.weak_region_len_factor - (opt.nb_correction_rounds - j - 1) * step_weak_region_len_factor;
                                            l_opt.max_len_weak_region1 = (j+1) * step_max_len_weak_region1;

                                            const pair<vector<pair<size_t, const_UnitigMap<UnitigData>>>, vector<pair<size_t, const_UnitigMap<UnitigData>>>> p = getSeeds(l_opt, dbg, v_in_read[i], v_in_qual[i], long_read_correct, hap_id, m_km_um, (j+1) != opt.nb_correction_rounds);

                                            pair<string, string> correction = correctSequence(dbg, l_opt, v_in_read[i], v_in_qual[i], p.first, p.second, long_read_correct, partitions, hap_id, hap_reads, max_km_cov);

                                            v_in_read[i] = move(correction.first);
                                            v_in_qual[i] = move(correction.second);
                                        }
                                    }
                                }

                                mutex_file_out.acquire();

                                size_t block_sz = 0;

                                const streampos sp = out.tellp();

                                for (size_t i = 0; i < v_in_read.size(); ++i) block_sz += writeCorrectedOutput(out, v_in_name[i], v_in_read[i], v_in_qual[i], k, (long_read_correct ? opt.trim_qual : 0));

                                v_block_w.push_back({(ticket_id << 32) | block_sz, sp});

                                mutex_file_out.release();
                            }

                            in_read_len = 0;

                            v_in_read.clear();
                            v_in_qual.clear();
                            v_in_name.clear();

                            v_in_read_raw.clear();
                            v_in_name_raw.clear();

                            in_read.clear();
                            in_read_raw.clear();
                        }
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        fp_raw.close();
    }

    outfile.close();
    fp.close();

    // File reordering -> makes sure output reads are in the same order as input reads. It's more IO but meh, what can I do...
    if (!long_read_correct && (opt.nb_threads != 1)) { // If one thread used, output file is same order as input file so no reordering needed

        if (opt.verbose) cout << "Ratatosk::search(): Reordering reads" << endl;

        const string fn_tmp_out = fn_out + ".tmp";

        ofstream outfile;
        ostream out(0);

        ifstream infile;
        istream in(0);

        size_t max_block_sz = 0;

        char* buffer = nullptr;

        auto cmpBlocks = [](const pair<size_t, streampos>& a, const pair<size_t, streampos>& b) {

            return ((a.first >> 32) < (b.first >> 32)); // 32 bits shifting is not necessary here but kept for readability
        };

        sort(v_block_w.begin(), v_block_w.end(), cmpBlocks);

        outfile.open(fn_tmp_out.c_str());
        out.rdbuf(outfile.rdbuf());

        if (check_files(fn_out, false, false)) {

            infile.open(fn_out.c_str(), std::ifstream::in);
            in.rdbuf(infile.rdbuf());
        }
        else {

            cerr << "Ratatosk::search(): Reordering could not happen because input file " << fn_out << " could not be open. Abort." << endl;
            exit(1);
        }

        for (const auto p : v_block_w) max_block_sz = max(max_block_sz, static_cast<size_t>(p.first & 0x00000000ffffffffULL));

        buffer = new char[max_block_sz];

        for (const auto p : v_block_w) {

            const size_t block_sz = p.first & 0x00000000ffffffffULL;

            in.seekg(p.second);
            in.read(buffer, block_sz);
            out.write(buffer, block_sz);
        }

        delete buffer;

        outfile.close();
        infile.close();

        if (remove(fn_out.c_str()) != 0) {

            cerr << "Ratatosk::search(): Reordering could not happen because file " << fn_out << " could not be removed. Abort." << endl;
            exit(1);
        }
        else if (rename(fn_tmp_out.c_str(), fn_out.c_str()) != 0) {

            cerr << "Ratatosk::search(): Reordering could not happen because file " << fn_tmp_out << " could not be renamed. Abort." << endl;
            exit(1);
        }
    }
}

int main(int argc, char *argv[]) {

    srand(time(NULL));

    Correct_Opt opt;

    if (argc < 2) PrintUsage(opt);
    else {

        const int parse = parse_ProgramOptions(argc, argv, opt);

        if (parse == 1) PrintVersion();
        else if (parse == 2) PrintUsage(opt);
        else if (parse == 3) PrintCitation();
        else {

            if (!check_ProgramOptions(opt)) return 0;

            Correct_Opt opt_pass1(opt);
            Correct_Opt opt_pass2(opt);

            bool dbg_empty = false;

            const bool hasGraph = (opt.fn_graph_in.length() != 0);
            const bool hasData = (opt.fn_index_in.length() != 0);
            const bool hasIndex = (hasGraph && hasData);

            string filename_out_extra_sr;

            string graph_long = opt.filename_long_out + ".index.k" + std::to_string(opt.k);
            string graph_short = opt.filename_long_out + ".index.k" + std::to_string(opt.small_k);

            // Establish if some reads are missing
            if (!opt.filenames_short_all.empty() && !hasIndex) {

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

                if (!hasIndex) {

                    if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Building graph from short reads (1/2)." << endl;

                    dbg.build(opt_pass1);
                    dbg.write(graph_long, opt_pass1.nb_threads, false, opt_pass1.verbose);

                    graph_long += ".fasta";
                    dbg_empty = (dbg.length() == 0);
                }

                if (!dbg_empty){

                    pair<HapReads, HapReads> hapPass1; // <short reads, long reads>

                    Roaring* partitions = nullptr;

                    opt_pass1.k = opt_pass1.small_k;
                    opt_pass1.filename_long_out += ".2";

                    dbg = CompactedDBG<UnitigData>(opt_pass1.k);

                    if (hasIndex) {

                        if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Loading graph (1/2)." << endl;

                        dbg.read(opt.fn_graph_in, opt_pass1.nb_threads, opt_pass1.verbose);

                        if (!readGraphData(opt.fn_index_in, dbg, false, opt_pass1.verbose)) {

                            cerr << "Ratatosk::Ratatosk(): Graph data could not be loaded. Abort." << endl;
                            dbg_empty = true;
                        }
                        else dbg_empty = (dbg.length() == 0);
                    }
                    else {

                        opt_pass1.filename_ref_in = vector<string>(1, graph_long);
                        opt_pass1.filename_seq_in = vector<string>();

                        dbg.build(opt_pass1);

                        opt_pass1.filename_ref_in = vector<string>();
                        opt_pass1.filename_seq_in = opt.filename_seq_in;

                        if (filename_out_extra_sr.length() != 0) opt_pass1.filename_seq_in.push_back(filename_out_extra_sr);

                        if (!opt.filenames_long_phase.empty()){

                            if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding phasing to graph (1/2)." << endl;

                            addPhasing(dbg, opt, hapPass1.first, false, false);

                            hapPass1.second.hapBlock2id = hapPass1.first.hapBlock2id;
                            hapPass1.second.hapType2id = hapPass1.first.hapType2id;

                            addPhasing(dbg, opt, hapPass1.second, true, false);
                        }

                        if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding colors and coverage to graph (1/2)." << endl;

                        addCoverage(dbg, opt_pass1, hapPass1.first, false);

                        if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding SNPs candidates to graph (1/2)." << endl;

                        detectSNPs(dbg, opt_pass1);

                        if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Adding micro/mini-satellites motif candidates to graph (1/2)." << endl;

                        detectShortCycles(dbg, opt_pass1);
                    }

                    if (opt.index) {

                        if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Writing index to disk (1/2)." << endl;

                        dbg.write(graph_short, opt_pass1.nb_threads, false, opt_pass1.verbose); // Write graph

                        writeGraphData(graph_short + ".rtsk", dbg, opt_pass1.verbose); // Write UnitigData
                    }
                    else if (!dbg_empty){ // correct

                        if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Correcting long reads (1/2)." << endl;

                        //opt_pass1.nb_threads = 1;

                        search(dbg, opt_pass1, false, partitions, hapPass1);
                    }
                    else if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Graph is empty, no correction can be done. Output uncorrected reads." << endl;

                    if (partitions != nullptr) delete[] partitions;
                }
                else if (opt_pass1.verbose) cout << "Ratatosk::Ratatosk(): Graph is empty, no correction can be done. Output uncorrected reads." << endl;
            }

            if (opt.pass2_only || !opt.pass1_only) {

                pair<HapReads, HapReads> hapPass2;

                CompactedDBG<UnitigData> dbg(opt_pass2.k);

                if (opt.pass2_only){

                    if (hasGraph) {

                        if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Loading graph (2/2)." << endl;

                        dbg.read(opt.fn_graph_in, opt_pass2.nb_threads, opt_pass2.verbose);
                    }
                    else {

                        if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Building graph (2/2)." << endl;

                        dbg.build(opt_pass2);
                    }

                    opt_pass2.filename_seq_in.clear();
                    opt_pass2.filename_seq_in.insert(opt_pass2.filename_seq_in.end(), opt_pass2.filenames_helper_long_in.begin(), opt_pass2.filenames_helper_long_in.end());
                    opt_pass2.filename_seq_in.insert(opt_pass2.filename_seq_in.end(), opt_pass2.filenames_long_in.begin(), opt_pass2.filenames_long_in.end());
                }
                else {

                    if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Building graph (2/2)." << endl;

                    opt_pass2.filename_seq_in.clear();
                    opt_pass2.filename_seq_in.insert(opt_pass2.filename_seq_in.end(), opt_pass2.filenames_helper_long_in.begin(), opt_pass2.filenames_helper_long_in.end());
                    opt_pass2.filename_seq_in.push_back(opt_pass1.filename_long_out + ".fastq");

                    opt_pass2.filenames_long_in.clear();
                    opt_pass2.filenames_long_in.push_back(opt_pass1.filename_long_out + ".fastq");

                    dbg.read(graph_long, opt_pass2.nb_threads, opt_pass2.verbose);

                    { // Remove long k-mer graph and unmapped reads temporary files from disk
                        if (!opt.pass1_only && !opt.pass2_only && (remove(graph_long.c_str()) != 0)) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;
                        if ((filename_out_extra_sr.length() != 0) && (remove(filename_out_extra_sr.c_str()) != 0)) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;

                        graph_long = ""; // Graph was read and deleted
                        filename_out_extra_sr = ""; // Unmapped reads were read and deleted
                    }
                }

                dbg_empty = (dbg.length() == 0);

                // Pass 2
                if (!dbg_empty || opt.pass2_only){

                    if (dbg.length() != 0) {

                        if (hasData) {

                            if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Loading unitig data (2/2)." << endl;

                            if (!readGraphData(opt.fn_index_in, dbg, false, opt_pass2.verbose)) cerr << "Ratatosk::Ratatosk(): Graph data could not be loaded. Abort." << endl;
                        }
                        else {

                            if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Adding colors and coverage to graph (2/2)." << endl;

                            addCoverage(dbg, opt_pass2, hapPass2.second, true);

                            if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Adding SNPs candidates to graph (2/2)." << endl;

                            detectSNPs(dbg, opt_pass2);

                            if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Adding micro/mini-satellites motif candidates to graph (2/2)." << endl;

                            detectShortCycles(dbg, opt_pass2);
                        }

                        if (opt.index) {

                            if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Writing index to disk (1/2)." << endl;

                            dbg.write(graph_long, opt_pass2.nb_threads, false, opt_pass2.verbose); // Write graph

                            writeGraphData(graph_long + ".rtsk", dbg, opt_pass2.verbose); // Write UnitigData
                        }
                        else { // correct

                            if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Correcting long reads (2/2)." << endl;

                            /*{
                                const string long_in_tmp = string(opt_pass2.filename_long_out + ".tmp");

                                phasing_test(dbg, opt_pass2, vector<string>(1, string("/nfs/odinn/tmp/guillaumeh/data/ont/test6/HG002_ONT.fq")), opt_pass2.filenames_long_in, long_in_tmp);

                                opt_pass2.filenames_long_in.clear();
                                opt_pass2.filenames_long_in.push_back(long_in_tmp);
                            }*/

                            search(dbg, opt_pass2, true, nullptr, hapPass2);
                        }
                    }
                    else {

                        if (opt_pass2.verbose) cout << "Ratatosk::Ratatosk(): Graph is empty, no correction can be done. Output 1st correction pass reads. (2/2)" << endl;

                        copyFile(opt_pass2.filename_long_out + ".fastq", opt_pass2.filenames_long_in); // Graph is empty, LR cannot be corrected
                    }
                }
                else copyFile(opt_pass2.filename_long_out + ".fastq", opt.filenames_long_in); // Graph is empty, LR cannot be corrected
            }

            if (!opt.pass1_only && !opt.pass2_only) { // Clean up tmp files if correction was made all at once from the raw files

                if (remove(graph_long.c_str()) != 0) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;
                if (!dbg_empty && remove((opt_pass1.filename_long_out + ".fastq").c_str()) != 0) cerr << "Ratatosk::Ratatosk(): Couldn't remove temporary file" << endl;
            }
        }
    }
}