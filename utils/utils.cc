#include <iostream>
#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <optional>
#include <regex>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "aa.h"
#include "cdn.h"
#include "dna.h"
#include "defines.h"

#ifdef DSA_TARGET_WIN64
#include "local-getopt.h"
#elif defined(DSA_TARGET_LINUX)
#include <getopt.h>
#endif

namespace fs = std::filesystem;

using bio::Aas;
using bio::Cdns;

void print_usage(std::ostream &);
void run_extract_aas (int, char *argv[]);
void run_extract_cdns(int, char *argv[]);
void run_venn_diagram(int, char *argv[]);
void run_print_help(int, char *argv[]);

std::optional<bio::Cdns> format_ascii(std::string &&, const std::optional<std::regex> &);
std::optional<bio::Cdns> format_horizontal(std::string &&, const std::optional<std::regex> &);

template<typename InputIter>
std::string
join(std::string delim, InputIter begin, InputIter end) {
    std::stringstream ss;
    if (begin == end) return ss.str();
    for (;;) {
        ss << *begin;
        ++begin;
        if (begin == end) break;
        ss << delim;
    }
    return ss.str();
}

void
run_print_help(int argc, char *argv[]) {
std::cout <<
"dsa-utils COMMAND [OPTIONS...]\n"
"  Recognized COMMANDs: extract_aas, venn\n"
"    extract_aas : Open one or more dsa output files compile create a list of\n"
"                  unique amino sequences. Optionally filter and/or capture\n"
"                  using a regular expression. Optionally add a label column\n"
"                  to the output.\n"
"      OPTIONS:\n" 
"        --label (-l) LABEL\n"
"          Output will be two columns separated by a tab character.\n"
"            Column 1 will contain LABEL.\n"
"            Column 2 will contain the amino acid sequences.\n"
"        --regex (-r) REGEX\n"
"          Search sequences for REGEX and discard if REGEX is not found.\n"
"          REGEX may optinally contain exactly 1 capture group. In this case\n"
"          the contents of the capture group will be returned in the sequence\n"
"          column.\n"
"        --output (-o) FILENAME\n"
"          Write output to FILENAME. If -o is not used, output will be printed\n"
"          to stdout.\n"
"      EXAMPLE: Combine HCDR3s from two dsa files, label as 'control', and\n"
"               print list of unique sequences.\n"
"        dsa-util -l control -r \"[YF][YF]C(.*)WG.G\" dsa1.csv dsa2.csv\n" 
"\n"
"    venn : Accepts sets of labeled sequences (e.g., as created by extract_aas)\n"
"           from stdin and calculates all set intersections. Results are printed\n"
"           to stdout in a multicolumn format where headers are the labels from\n"
"           the input files and the label columns contain 0 or 1 depending on\n"
"           whether the sequence was found in that labeled population.\n"
"      OPTIONS:\n"
"        --include_summary\n"
"          Output will include a summary table with the count and percent of\n"
"          seqeunces shared among all combinations of labeled populations.\n"
"        --omit_sequences\n"
"          Suppress printing out the sequences themselves (if, for exampe, only\n"
"          the summary information is required). The summary table has a 0/1\n"
"          column for each label, an 'of N' colum where N is the total number\n"
"          of unique sequences, and a 'percent' column showing 100*N/total.\n"
"      EXAMPLE: extract HCDR3s from control and experimental datasets and find\n"
"             which are shared between them and which are not.\n"
"        { dsa-util -l control -r \"[YF][YF]C(.*)WG.G\" dsa1.csv ; \\\n"
"          dsa-util -l exptl   -r \"[YF][YF]C(.*)WG.G\" dsa2.csv ; } | \\\n"
"          dsa-util --include_summary venn" << std::endl;
}

const static std::map<std::string_view, void(*)(int, char *[])> COMMAND_RUNNERS = {
    {"extract_aas",  run_extract_aas },
    //{"extract_cdns", run_extract_cdns},
    {"venn",         run_venn_diagram},
    {"--help",       run_print_help}
};

const static std::map<std::string_view, 
    std::optional<bio::Cdns>(*)(std::string &&, const std::optional<std::regex> &)> CDN_FORMATTERS = {
    {"ascii",      nullptr},
    {"horizontal", nullptr}
};

void
run_venn_diagram(int argc, char *argv[]) {
    int include_summary = 0, omit_sequences = 0;
    
    const char *opt_chars = "";

    static struct option long_options[] = {
        {"include_summary", no_argument, &include_summary, 1},
        {"omit_sequences",  no_argument, &omit_sequences,  1},
        {               0,            0,               0,  0}
    };

    for (;;) {
        int option_index = 0;
        int c = getopt_long(argc, argv, opt_chars, long_options, &option_index);
        if (c == -1) break;
        switch (c) {
            case 0: // a long option with no single char version 
                break;
            case '?':
                std::cerr << "unrecognized option: -" << optopt << std::endl;
                break;
            case ':':
                std::cerr << "missing required argument for -" << optopt << std::endl;
                exit(EXIT_FAILURE);
                break;
            default:
                std::cerr << "unrecognized option -" << optopt << std::endl;
                exit(EXIT_FAILURE); //should never happen
        }
    }

    std::unordered_map<std::string, size_t> labels;
    std::unordered_map<std::string, std::unordered_set<size_t>> venn;
    for (std::string line; std::getline(std::cin, line); ) {
        while (!line.empty() && std::isspace(line.back())) line.pop_back();
        size_t tab = line.find('\t');
        std::string label    = line.substr(0, tab);
        std::string sequence = line.substr(tab+1);

        auto [ll, _] = labels.insert({label, labels.size()});
        size_t label_index = ll->second;

        auto vv = venn.find(sequence);
        if (vv == venn.end()) {
            venn.insert(std::make_pair(std::move(sequence), std::unordered_set<size_t>{label_index}));
        } else {
            vv->second.insert(label_index);
        }
    }

    std::vector<std::pair<std::string, std::vector<bool>>> sorted;
    while (!venn.empty()) {
        auto nh = venn.extract(venn.begin());
        std::vector<bool> inclusion(labels.size(), false);
        for (size_t label : nh.mapped()) inclusion[label] = 1;
        sorted.push_back(std::make_pair(std::move(nh.key()), std::move(inclusion)));
    }

    struct inclusion_cmp {
        bool operator()(const std::pair<std::string, std::vector<bool>> &a,
                        const std::pair<std::string, std::vector<bool>> &b) const {
            auto sum_a = std::count(a.second.begin(), a.second.end(), true);
            auto sum_b = std::count(b.second.begin(), b.second.end(), true);
            if (sum_b < sum_a) return true;
            if (sum_a < sum_b) return false;
            return b.second < a.second;
        }
    };

    std::sort(
        sorted.begin(),
        sorted.end(),
        inclusion_cmp()
    );

    struct DerefHash {
        size_t operator()(const std::vector<bool> *p) const { return std::hash<std::vector<bool>>{}(*p); }
    };

    struct DerefEquals {
        bool operator()(const std::vector<bool> *p, const std::vector<bool> *q) const { return *p == *q; }
    };

    std::vector<std::string> indexed_labels(labels.size());
    for ( const auto & [k, v] : labels ) indexed_labels [ v ] = k;

    if (include_summary) {
        std::unordered_map<const std::vector<bool> *, size_t, DerefHash, DerefEquals> counts;
        size_t total = 0;
        for (const auto &[_, inclusion] : sorted) {
            counts[&inclusion] += 1;
            total += 1;
        }

        std::cout << join("\t", indexed_labels.cbegin(), indexed_labels.cend()) << "\tof " << total << "\tpercent" << std::endl;

        std::vector<bool> perms(labels.size(), false);
        for (size_t n=0; n<labels.size(); ++n) {
            std::fill(perms.begin(),   perms.begin()+n, false);
            std::fill(perms.begin()+n, perms.end(),     true );
            do {
                for (size_t i=0; i<perms.size(); ++i) std::cout << perms[i] << '\t';
                std::cout <<  counts[&perms] << '\t'
                          << (counts[&perms] / static_cast<double>(total) * 100) << '%' 
                          << std::endl;
            } while (std::next_permutation(perms.begin(), perms.end()));
        }
        std::cout << std::endl;
    }

    if (!omit_sequences) {
        std::cout << join("\t", indexed_labels.cbegin(), indexed_labels.cend()) << "\tsequence" << std::endl;
        for (const auto &[sequence, inclusion] : sorted) {
            std::cout << join("\t", inclusion.begin(), inclusion.end()) << "\t" << sequence << std::endl;
        }
    }
}

void
run_extract_aas(int argc, char *argv[]) {
    const char * opt_chars = "l:o:r:";
    int omit_label = 0;
    std::string regex_string, output_filename, label;
    std::vector<std::string> input_filenames;

    static struct option long_options[] = {
        //options
        {"label",       required_argument, 0, 'l'},
        {"regex",       required_argument, 0, 'r'},
        {"output",      required_argument, 0, 'o'}
    };

    for (;;) {
        int option_index = 0;
        int c = getopt_long(argc, argv, opt_chars, long_options, &option_index);
        if (c == -1) break;
        switch (c) {
            //case 0 indicates a long option with no corresponding single-letter flag
            //case 0:
            //    break;
            case 'l':
                label = optarg;
                break;
            case 'r':
                regex_string = optarg;
                break;
            case 'o':
                output_filename = optarg;
                break;
            case '?':
                std::cerr << "unrecognized option: -" << optopt << std::endl;
                break;
            case ':':
                std::cerr << "missing required argument for -" << optopt << std::endl;
                exit (EXIT_FAILURE);
                break;
            default:
                exit (EXIT_FAILURE); //should never happen
        }
    }

    for (; optind != argc; ++optind) input_filenames.push_back(argv[optind]);

    if (input_filenames.empty()) {
        std::cerr << "No dsa input files specified" << std::endl << std::endl;
        print_usage(std::cerr);
        exit (EXIT_FAILURE);
    }

    std::optional<std::regex> rgx;
    if (!regex_string.empty()) {
        try {
            rgx = std::regex(regex_string);
            if (rgx->mark_count() > 1) {
                std::cerr << "Regexes are limited to only one capture group" << std::endl;
                exit (EXIT_FAILURE);
            }
        } catch (std::regex_error &ex) {
            (void)ex;
            std::cerr << "-r " << regex_string << " could not be interpreted as a regular expression" << std::endl << std::endl;
            print_usage(std::cerr);
            exit (EXIT_FAILURE);
        }
    }

    std::unordered_set<Aas> unique_aas;
    for (size_t i=0; i<input_filenames.size(); ++i) {
        const size_t label = i;
        const std::string &filename = input_filenames[i];

        std::ifstream ifs(filename);
        if (!ifs) {
            std::cerr << "Could not open '" << filename << "' for reading" << std::endl << std::endl;
            exit (EXIT_FAILURE);
        }

        size_t line_no = 0;
        for (std::string line; std::getline(ifs, line); ) {
            ++line_no;
            if (line.find("#Unique Amino Acids") != std::string::npos) {
                std::getline(ifs, line); //skip headers
                break;
            }
        }

        for (std::string line; std::getline(ifs, line); ) {
            ++line_no;
            if (line.find('#') != std::string::npos) break;
            while (!line.empty() && std::isspace(line.back())) line.pop_back();
            size_t tab = line.find('\t', 0);
            tab = line.find('\t', tab+1);
            if (tab == std::string::npos) {
                std::cerr << "Bad formatting at '" << filename << "' line " << line_no << ":" << std::endl;
                std::cerr << "  Expected protein sequence in column 3 but got '" << line << "' instead" << std::endl;
                exit (EXIT_FAILURE);
            }
            Aas aas(line.c_str()+tab+1);
            if (line.size()-(tab+1) != aas.size()) {
                std::cerr << "Bad formatting at '" << filename << "' line " << line_no << ":" << std::endl;
                std::cerr << "  Protein sequence contained invalid characters" << std::endl;
                exit(EXIT_FAILURE);
            }
            if (rgx) {
                std::match_results<std::string::const_iterator> match;
                const std::string::const_iterator begin = line.cbegin()+tab+1;
                const std::string::const_iterator end   = line.cend();
                if (!std::regex_search(begin, end, match, *rgx)) continue;
                if (rgx->mark_count() != 0) {
                    auto sm = match[1];
                    auto left = sm.first - begin, right = end - sm.second;
                    aas.exo(left, right);
                }
            }
            unique_aas.insert(std::move(aas));
        }
        ifs.close();
    }

    std::ofstream ofs;
    if (!output_filename.empty()) {
        ofs.open(output_filename);
        if (!ofs) {
            std::cerr << "Could not open '" << output_filename << "' for writing" << std::endl;
            exit (EXIT_FAILURE);
        }
    }

    std::ostream &os = output_filename.empty() ? std::cout : ofs;
    for (const Aas &aas : unique_aas) {
        if (!label.empty()) os << label << '\t';
        os << aas << std::endl;
    }
    if (ofs.is_open()) ofs.close();
}

void
run_extract_cdns(int argc, char *argv[]) {
    static struct option long_options [ ] = {
        {"format", required_argument, 0, 'f'},
        {"output", required_argument, 0, 'o'},
        {"regex",  required_argument, 0, 'r'}
    };

    const char * opt_chars = "f:o:r:";
    std::string format_string = "ascii", regex_string, output_filename;
    std::vector<std::string> input_filenames;

    for (;;) {
        int option_index = 0;
        int c = getopt_long(argc, argv, opt_chars, long_options, &option_index);
        if ( c == -1 ) break;
        switch ( c ) {
            //case 0 indicates a long option with no corresponding single-letter flag
            //case 0:
            //    break;
        case 'f':
            format_string = optarg;
            break;
        case 'o':
            output_filename = optarg;
            break;
        case 'r':
            regex_string = optarg;
            break;
        case '?':
            std::cerr << "unrecognized option: -" << optopt << std::endl;
            break;
        case ':':
            std::cerr << "missing required argument for -" << optopt << std::endl;
            exit(EXIT_FAILURE);
            break;
        default:
            exit(EXIT_FAILURE); //should never happen
        }
    }


    if (CDN_FORMATTERS.find(format_string) == CDN_FORMATTERS.end()) {
        std::cerr << "Unrecognized format type: '" << format_string << "'" << std::endl;
        exit (EXIT_FAILURE);
    }
    auto formatter = CDN_FORMATTERS.at(format_string);

    for ( ; optind != argc; ++optind ) input_filenames.push_back(argv [ optind ]);

    if (input_filenames.empty()) {
        std::cerr << "No dsa input files specified" << std::endl << std::endl;
        print_usage(std::cerr);
        exit(EXIT_FAILURE);
    }

    std::optional<std::regex> rgx;
    if (!regex_string.empty()) {
        try {
            rgx = std::regex(regex_string);
            if (rgx->mark_count() > 1) {
                std::cerr << "Regexes are limited to only one capture group" << std::endl;
                exit(EXIT_FAILURE);
            }
        } catch (std::regex_error & ex) {
            (void)ex;
            std::cerr << "-r " << regex_string << " could not be interpreted as a regular expression" << std::endl << std::endl;
            print_usage(std::cerr);
            exit(EXIT_FAILURE);
        }
    }

    std::unordered_set<Cdns> unique_cdns;
    for (const std::string & filename : input_filenames) {
        std::ifstream ifs(filename);
        if ( !ifs ) {
            std::cerr << "Could not open '" << filename << "' for reading" << std::endl << std::endl;
            exit(EXIT_FAILURE);
        }

        size_t line_no = 0;
        for ( std::string line; std::getline(ifs, line); ) {
            ++line_no;
            if ( line.find("#Unique Codons") != std::string::npos ) {
                std::getline(ifs, line); //skip headers
                break;
            }
        }

        for ( std::string line; std::getline(ifs, line); ) {
            ++line_no;
            if (line.find('#') != std::string::npos) break;
            while (!line.empty() && std::isspace(line.back())) line.pop_back();
            size_t tab = line.find('\t', 0);
            tab = line.find('\t', tab + 1);
            if (tab == std::string::npos) {
                std::cerr << "Bad formatting at '" << filename << "' line " << line_no << ":" << std::endl;
                std::cerr << "  Expected codon sequence in column 3 but got '" << line << "' instead" << std::endl;
                exit(EXIT_FAILURE);
            }
            Cdns cdns(line.c_str() + tab + 1);
            if (line.size() - (tab + 1) != cdns.size()) {
                std::cerr << "Bad formatting at '" << filename << "' line " << line_no << ":" << std::endl;
                std::cerr << "  Codon sequence contained invalid characters" << std::endl;
                exit(EXIT_FAILURE);
            }
            if (rgx) {
                Aas aas = cdns;
                std::string_view cdns_view = cdns.as_string_view();
                std::string_view  aas_view = aas.as_string_view();

                std::match_results<std::string_view::const_iterator> match;
                if (!std::regex_search(aas_view.begin(), aas_view.end(), match, *rgx)) continue;
                if ( rgx->mark_count() != 0 ) {
                    auto sm = match [ 1 ];
                    auto left = sm.first - aas_view.begin(), right = aas_view.end() - sm.second;
                    cdns.exo(left, right);
                }
            }
            unique_cdns.insert(std::move(cdns));
        }
        ifs.close();
    }

    std::ofstream ofs;
    if ( !output_filename.empty() ) {
        ofs.open(output_filename);
        if ( !ofs ) {
            std::cerr << "Could not open '" << output_filename << "' for writing" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    std::ostream & os = output_filename.empty() ? std::cout : ofs;
    if (format_string == "ascii") {
        for (const Cdns & cdns : unique_cdns) {
            os << cdns << std::endl;
        }
    } else if (format_string == "horizontal") {
        for (const Cdns & cdns : unique_cdns) {
            for (bio::Cdn cdn : cdns) os << cdn.p1() << cdn.p2() << cdn.p3();
            os << std::endl;
        }
    }

    if (ofs.is_open()) ofs.close();
}

void
print_usage(std::ostream & os) {
    os << "dsa-utils COMMAND [OPTIONS]" << std::endl;
    os << "\tRecognized commands: ";
    for (const auto &[cmd, _] : COMMAND_RUNNERS) os << cmd << " ";
    os << std::endl;
}

int
main(int argc, char *argv[]) {
    if (argc < 2) {
        print_usage(std::cerr);
        return EXIT_FAILURE;
    }

    try {
        COMMAND_RUNNERS.at(argv[1])(argc-1, argv+1);
    } catch (std::out_of_range &ex) {
        (void)ex;
        std::cerr << "Unrecognized command: '" << argv[1] << "'" << std::endl << std::endl;
        print_usage(std::cerr);
        return EXIT_FAILURE;
    }

    return 0;
}