#include "help.h"

#include "cdn.h"
#include "dna.h"
#include "defines.h"

#include <regex>

namespace help {

using bio::Cdn;
using bio::Nt;

void
print_option(std::ostream &os, const char *option, const char *helptext) {
    const int OPTW = 18;
    static const char *INDENT = "  ";
    os << INDENT << std::setw(OPTW) << std::left << option << INDENT << helptext << std::endl;
}

void
print_help() {
    const struct OptHelp options[] = {
        { 0,  "help",           "print this message and exit, see also (--help codons) and (--help templates)"},
        { 0,  "version",        "print the program version number"},
        {'s', "no_header",      "suppress printing of headers in program output"},
        {'f', "fw_ref",         "nucleotide sequence(s) used to determine UMI and reading frame for forward read"},
        {'r', "rv_ref",         "nucleotide sequence(s) used to determine UMI and reading frame for reverse read"},
        {'t', "template",       "amino acid sequence to which translated paired-end reads will be aligned (or 'none' for no alignment)"},
        {'d', "template_dna",   "dna sequence which will be traslated for alignmet with translated paired-end reads"},
        {'q', "min_qual",       "bases with quality scores of < min_qual will be removed from 3' ends of reads (default=A)"},
        {'x', "skip_assembly",  "skip paired read assemly and align forward and reverse reads to template independently (off by default)"},
        {'v', "min_overlap",    "minimum 3' overlap required for assembly of paired ends (default=9)"},
        {'m', "max_mismatch",   "maximum allowable nucleotide mismatches in paired 3' ends (default=0)"},
        {'g', "min_umi_grp",    "during umi collapse, sequences with < min_umi_grp members will be discarded (default=1)"},
        {'a', "min_aln",        "reads where (alignment score / max possible alignment score) < min_aln will be discarded (default=0.8)"},
        {'n', "number_from",    "number template amino acids starting from 'n' in the #Substitutions# section (default=1)"},
        {'c', "show_codons",    "output format for codons; can be ascii, horizontal, vertical, or none (none by default)"},
        { 0 , "split",          "regular expression to split translated ORFs into multiple pieces for alignment to separate templates (see --help templates)"},
        { 0 , "template_db",    ".fasta file containing a list of possible nucleotide templates for split sequences (see --help templates)"},
        { 0 , "trim"            "trim the N- and/or C-terminal ends of a template or template database to match the deep-sequenced region (default=0,0)"}
    }; 

    std::cout << "Deep Sequencing Analysis version " << VERSION_STRING << "\n\n"
                 "Program usage: dsa [options] [-f forward_reference] [-r reverse_reference]\n"
                 "  [-t template] forward_reads.fastq reverse_reads.fastq > output.csv\n" << std::endl;
    std::cout << "Aligns paired reads in fastq files forward_reads.fastq and reverse_reads.fastq,\n"
                 "  extracts UMI barcodes, translates, and aligns the translated sequence to the\n"
                 "  supplied amino acid or dna template seuqnece(s)." << std::endl;
    std::cout << "Reference sequences (-f, -r) are use to identify UMI barcodes and the open reading frame.\n"
                 "  References are composed of capital ATGCN and lowercase n characters. Capital ATGCN are\n"
                 "  used to match bases in the reads (N is a wild-card) and lowercase n characters define\n"
                 "  the UMI barcode. The end of the reference sequence defines the reading frame used for\n"
                 "  translation. More than one reference sequence can be defined by using the -f or -r\n"
                 "  arguments multipe times. Reference sequences will be tested in the order given and the\n"
                 "  first match will be accepted.\n"
                 "Example:\n"
                 "  Reference:    GAAnnCGnnNNN\n"
                 "  Fw Read:   AACGAAGACGAGGTTCTGCAGCCGCGGCTGGAGGCGGGGGTGTAGT\n"
                 "  Barcode:         GA  AG\n"
                 "  ORF:                      CTGCAGCCG...\n"
                 "                            LeuGlnPro..." << std::endl;
    std::cout << "\nOUTPUT:\n"
              << "Output is printed as tab-delimited text to the terminal stanard output stream.\n"
              << "To write to a file, use output redirection (e.g. \"dsa ... > output.csv\").\n"
              << "Program output is divided into several sections:\n"
              << "\n#Settings# lists the values of the input parameters.\n"
              << "\n#Parse# shows the numbers of sequences that were removed by quality control.\n"
              << "\n#Templates# lists the amino acid template sequence(s) used for alignments.\n"
              << "  Column 1 contains the template ID number\n"
              << "  Column 2 shows the name of the template\n"
              << "  Column 3 shows the amino acid sequence of the template\n"
              << "  Alignment to multiple templates is supported. See (--help templates) for details\n"
              << "\n#Template Usage# gives frequency statistics for which templates were used for alignments.\n"
              << "   This is used, for example, with template databases (see -- help templates) to determine\n"
              << "   frequencies of particular V segments, J segments, constant regions, etc.\n"
              << "  Column 1 contains the portion of the split (see --split)\n"
              << "  Column 2 contains the identifier for the template matched to that split\n"
              << "  Column 3 contains the number of matching UMI groups\n"
              << "  Column 4 contains the frequency of matching UMI groups\n"
              << "\n#Alignments# lists the alignments of the translated paired reads to the templates.\n"
              << "  One Alignments section will be produced per template.  If --show_codons was requested\n"
              << "  (see --help codons), every amino acid alignment is followed by its corresponding codon\n"
              << "  sequence.\n"
              << "  Column 1 contains the ID number of the amino template used for alignment\n"
              << "           or blank for the codon sequence.\n"
              << "  Column 2 contains the number of sequences merged during UMI collapse\n"
              << "           or blank for the codon sequence.\n"
              << "  Column 3 contains the UMI barcode\n"
              << "           or blank for the codon sequence.\n"
              << "  Column 4 contains the aligned amino acid sequence (see below for format)\n"
              << "           or the codon sequence if --show_codons was requested.\n"
              << "  Amino acid alignments are formatted as follows:\n" 
              << "    1. Capital letters indicate a match or mismatch with the template.\n"
              << "    2. A '-' character represents a deletion in the read relative to the template.\n"
              << "    3. Lower case characters show insertions relative to the template.\n"
              << "    Example:\n"
              << "       Template   MATIHKA\n"
              << "      Alignment   asML-VHqKA\n"
              << "    Should be interpreted as:\n"
              << "       Template --MATIH-KA\n"
              << "                  |: :| ||\n"
              << "      Alignment ASML-VHQKA" << std::endl;
    std::cout << "\n#Substitutions# contains a grid of amino acid mutation frequencies relative to the template.\n"
              << "  Column headers show the numbered residues of the template (see --number_from).\n"
              << "  Row headers show the different possible mutations.\n"
              << "  Each cell shows the frequency with with the amino acid indicated by the row header\n"
              << "    was found in the position of the amino acid indicated by the column header.\n"
              << "  The frequency of the wild type amino acid is always set to 0.0 to aid\n"
              << "    construction of stacked bar charts. The true frequency of the wild type\n"
              << "    residue is therefore 1.0-(sum of frequencies in a given column)." << std::endl;
    std::cout << "\n#Mutation Counts# provides the count of synonymous vs nonsynonymous mutations for each\n"
              << "  residue in the amino acid template. This section requires a template dna sequence as input\n"
              << "  (see -d, --template_dna). Indels in the aligned reads are not counted. Columns\n"
              << "  correspond to the amino acid residues in the translation of the dna template.\n"
              << "  Total:      the number of UMI groups with a non-indel at this position\n"
              << "  Non-Coding: the number of synonymous mutations\n"
              << "  Coding:     the number of non-synonymous muations\n" << std::endl;
    std::cout << "\n#Unique# shows a list of unique amino acid sequences and the corresponding number of unique\n"
              << "PCR events (UMI groups) and total reads (sum of UMI group sizes) for each.\n"
              << "Requires assembly of the paired ends (i.e. cannot be output when -x is set).\n"
              << "  Column 1 contains the number of UMI groups encoding this sequence\n"
              << "  Column 2 contains the number of PCR reads encoding this sequence\n"
              << "  Column 3 contains the amino unique acid sequence" << std::endl;
    std::cout << "\nOPTIONS:" << std::endl;
    print_opthelp(options);
}

void
print_help_codons() {
    std::cout << "Three codon output formats (-c, --show_codons) are available:\n"
              << "  1. 'horizontal' writes the nucleotides in left-to-right order"
              << "  2. 'vertical' writes each nucleotide triplet on the three lines\n"
              << "      beneath each amino acid in top-to-bottom order\n"
              << "  3. 'ascii' writes each codon as a single character according to\n"
              << "      the scheme outlined below\n" << std::endl;
    std::cout << "To convert an ASCII character codon to nucleotides:" << std::endl
              << "  1. subtract 48 from the decimal value (see www.asciitable.com)" << std::endl
              << "  2. deconstruct the binary value of resulting byte as follows:" << std::endl
              << "     bits 0 and 1 are ignored" << std::endl
              << "     bits 2 and 3 encode nucleotide #1" << std::endl
              << "     bits 4 and 5 encode nucleotide #2" << std::endl
              << "     bits 6 and 7 encode nucleotide #3" << std::endl
              << "     According to the following chart:" << std::endl
              << "       Binary DNA" << std::endl
              << "         00    A " << std::endl
              << "         01    C " << std::endl
              << "         10    T " << std::endl
              << "         11    G " << std::endl
              << "Example:" << std::endl
              << "  ASCII codon ';' has a decimal value of 59" << std::endl
              << "  59 - 48 = 11" << std::endl
              << "  11 in binary is 00001011" << std::endl
              << "                  ^^        ignored" << std::endl
              << "                    ^^      nucleotide #1 is 00 = A" << std::endl
              << "                      ^^    nucleotide #2 is 10 = T" << std::endl
              << "                        ^^  nucleotide #3 is 11 = G" << std::endl
              << "  so ';' is ATG, the start codon." << std::endl;
    std::cout << std::endl;

    std::cout << std::setw(6) << std::left  << "ASCII" 
              << std::setw(4) << std::left  << "DNA"
              << std::setw(4) << std::right << "AA"
              << std::endl;

    const Nt nts[4] = {Nt::A, Nt::C, Nt::G, Nt::T};
    char dna[] = "AAA";
    for (Nt n1 : nts)
    for (Nt n2 : nts)
    for (Nt n3 : nts) {
        Cdn nnn(n1, n2, n3);
        dna[0] = n1;
        dna[1] = n2;
        dna[2] = n3;
        bio::Aa res = bio::StandardTranslationTable.translate(nnn);
        std::cout << std::setw(6) << std::left  << nnn 
                  << std::setw(4) << std::left  << dna
                  << std::setw(4) << std::right << res
                  << std::endl;
    }
}

void
print_help_templates() {
    //            1                                                                             80
    std::cout << "dsa can handle multi-template alignments in different ways.\n"
              << "'Splitting' (--split) refers to dividing up each translated read into sections\n"
              << "  for alignment to different templates, or template databases. Splitting uses\n"
              << "  a regular expression based on the amino acid sequence.\n"
              << "  Examples of regular expressions for --split:\n"
              << "    1. --split=\"(.+[YF][YF]C..)(.+WG.G).+\" will divide an antibody VH\n"
              << "       segment into a V region (...YYCAR) and HCDR3 (...WGXG)\n"
              << "    2. --split=\"(.+[YF][YF]C..)(.+WG.G)(.+)\" will divide an antibody VH\n"
              << "       segment into a V region, HCDR3, and a portion of CH1\n"
              << "    3. --split=\".+([YF][YF]C.+WG.G).+\" will extract just the HCDR3 including\n"
              << "       the full YYCAR and WGXG motifs\n"
              << "    4. --split=\"(.{50})(.+)\" divides each read into one section containing\n"
              << "       the first 50 amino acids and a second section for the remainder of the\n"
              <<"        amino acids\n" << std::endl;
    std::cout << "A 'template database' (--template_db) is a .fasta file with the following\n"
              << "  properties:\n"
              << "  1. the sequences are nucleotides\n"
              << "  2. the nucleotides define an open reading frame\n"
              << "  3. headers are in IMGT format ('|' delimited with identifier in column 2)\n"
              << "     or a single identifier token following the '>' character.\n"
              << "     Example IMGT-format header for IGHV1-12*01:\n"
              << "       >AC090843|IGHV1-12*01|Mus musculus_C57BL/6|F|V-REGION|...\n"
              << "     Example alternatively formatted, single-identifier header for IGHV1-12*01\n"
              << "       >IGHV1-12*01\n" << std::endl;
    std::cout << "Templates can also be 'trimmed' (--trim) to limit the alignment to a sub-region\n"
              << "  of the original template. This is particularly useful for IMGT template\n"
              << "  databases (--template_db) because PCR primers will often be designed against\n"
              << "  interior regions of V segments, CH1, etc. --trim expects two integers\n"
              << "  separated by a comma, e.g. --trim=10,5 will remove 10 residues or codons\n"
              << "  from the left of each template and 5 residues from the right of each\n"
              << "  template prior to alignment. One instance of --trim must be supplied for\n"
              << "  every occurrence of --template, --template_dna, or --template_db given on\n"
              << "  the command line. Use --trim=0,0 to skip trimming a particular template.\n" << std::endl;
    std::cout << "For example, antibodies with diverse VH sequences, engineered to express a\n"
              << "  common HCDR3 can be divied into V region and HCDR3 using --split. dsa can\n"
              << "  determine the closest matching V region from a list of human or mouse V\n"
              << "  regions (--template_db) for each read, align the read to that V region\n"
              << "  while aligning the HCDR3 of the read to a different template (--template,\n"
              << "  --template_dna, or --template_db), concatenate the two alignments, and\n"
              << "  report the results.\n\n"
              << "A command to perform the multi-alignment proceedure outlined above might look\n"
              << "  like the following:\n\n"
              << "  $dsa --split=\"(.+[YF][YF]C..)(.+)WG.G.*\" \\\n"
              << "  $ --template_db=mouse_v_regions_imgt.fasta \\\n"
              << "  $ --trim=27,0 \\\n"
              << "  $ --template=RSEFYYYGNTYYYSAMDY \\\n"
              << "  $ --trim=0,0 \\ \n"
              << "  $ -f XXXXXXXXXXXXX -r XXXXXXXXXXXXX fw_reads.fastq rv_reads.fastq\n\n"
              << "  where -f and -r are appropriate reference sequences and the amino acid\n"
              << "  sequence supplied to --template is that of the engineered HCDR3. 27 amino acids\n"
              << "  will be removed from the N-terminus of each mouse V region in\n"
              << "  mouse_v_regions_imgt.fasta prior to alignment. The HCDR3 will be aligned to the\n"
              << "  amino acid sequence RSEFYYY..." << std::endl;
    std::cout << "Note, not everything needs to be aligned to a template or template database. To\n"
              << " skip alignment of a sequence or part of a split, use --template=none. For\n"
              << " example, to modify the command above to identify and align to V regions\n"
              << " from an IMGT database but not try to align the HCDR3 to any particular\n"
              << " template sequence, use the following:\n\n"
              << "  $dsa --split=\"(.+[YF][YF]C)(.+)WG.G.*\" \\\n"
              << "  $ --template_db=mouse_v_regions_imgt.fasta \\\n"
              << "  $ --trim=27,0 \\\n"
              << "  $ --template=none \\\n"
              << "  $ -f XXXXXXXXXXXXX -r XXXXXXXXXXXXX fw_reads.fastq rv_reads.fastq\n\n"
              << std::endl;
}

Params
parse_argv(int argc, char **argv) {
    Params p;

    static struct option long_options[] = {
        //flags
        {"no_header",      no_argument, &p.no_header_flag,      1},
        {"skip_assembly",  no_argument, &p.skip_assembly_flag,  1},
        //options  
        {"min_aln",        required_argument, 0, 'a'}, //minimum alignment score (fraction of max)
        {"fw_ref",         required_argument, 0, 'f'}, //forward UMI/reference DNA sequence
        {"min_umi_grp",    required_argument, 0, 'g'}, //minimum UMI group size
        {"max_mismatch",   required_argument, 0, 'm'}, //maximum mismatches in overlapped read ends
        {"min_qual",       required_argument, 0, 'q'}, //minimum quality score for 3' bases
        {"rv_ref",         required_argument, 0, 'r'}, //reverse UMI/reference DNA sequence
        {"template",       required_argument, 0, 't'}, //template AA sequence
        {"template_dna",   required_argument, 0, 'd'}, //template DNA sequence
        {"min_overlap",    required_argument, 0, 'v'}, //mimimum overlap in nucleotide for paired ends
        {"number_from",    optional_argument, 0, 'n'}, //number template residues in subs matrix starting from n
        {"show_codons",    required_argument, 0, 'c'}, //how to output codons
        {"split",          required_argument, 0,  0 }, //for a split template, e.g. V region CDR3, J/CH
        {"template_db",    required_argument, 0,  0 }, //file containing multiple templates
        {"trim",           required_argument, 0,  0 },
        //commands  
        {"version",        no_argument,     0,  0 }, //print version number
        {"help",           optional_argument, 0,  0 },
        //end-of-list sentinel
        {nullptr,        0,                        0,  0 }  //end-of-list sentinel
    };

    const char *opt_chars = "f:g:r:t:d:a:b:u:q:v:m:n:c:svx";

    std::regex trim_regex(R"-(([0-9]+),([0-9]+))-");
    std::smatch match;
    std::string optstring;
    std::optional<CodonOutput> co;

    for (;;) {
        int option_index = 0;
        int c = getopt_long(argc, argv, opt_chars, long_options, &option_index);
        if (c == -1) break;

        switch (c) {
            case 0: //case 0 indicates a long option with no corresponding 1-letter name
                if (long_options[option_index].flag != 0) break; //skip options that just set flags

                if (std::strcmp(long_options[option_index].name, "help") == 0) {
                    //this bizzare line gets the optional argument to --help if one is available
                    if (optarg == nullptr && optind < argc && argv[optind][0] != '-') optarg = argv[optind++];

                    if (optarg == nullptr) {
                        print_help();
                    } else if (std::strcmp(optarg, "codons") == 0) {
                        print_help_codons();
                    } else if (std::strcmp(optarg, "templates") == 0) {
                        print_help_templates();
                    } else {
                        std::cerr << "Unrecognized help topic: '" << optarg << "'" << std::endl;
                        exit (EXIT_FAILURE);
                    }
                    exit (EXIT_SUCCESS);
                } else if (std::strcmp(long_options[option_index].name, "version") == 0) {
                    std::cout << "Deep Sequencing Analysis version " << VERSION_STRING << std::endl;
                    exit (EXIT_SUCCESS);
                } else if (std::strcmp(long_options[option_index].name, "split") == 0) {
                    try {
                        p.split_template_regex = std::regex(optarg);
                        p.split_template_string = optarg;
                    } catch (std::regex_error &) {
                        std::cerr << "--split requires a valid regular expression; '"
                                  << optarg << "' could not be interpreted as one."
                                  << " Ensure the regular expression is correct and properly escaped for your shell."
                                  <<  std::endl;
                        exit (EXIT_FAILURE);
                    }
                } else if (std::strcmp(long_options[option_index].name, "template_db") == 0) {
                    p.template_sources.push_back(fs::path(optarg));
                } else if (std::strcmp(long_options[option_index].name, "trim") == 0) {
                    optstring = optarg;
                    if (!std::regex_match(optstring, match, trim_regex)) {
                        std::cerr << "--trim takes two comma-separated integers (e.g. --trim=5,0)" << std::endl;
                        exit (EXIT_FAILURE);
                    }
                    p.trims.push_back({std::stol(match.str(1)), std::stol(match.str(2))});
                }
                break;
            case 'a':
                p.min_alignment_score = std::strtof(optarg, nullptr);
                if (errno != 0 || p.min_alignment_score < 0.0 || p.min_alignment_score > 1.0) {
                    std::cerr << "min_alignment_score must be a number in the interval [0.0, 1.0]" << std::endl;
                    exit (EXIT_FAILURE);
                }
                break;
            case 'c':
                co = codon_output_from_string(optarg);
                if (!co) {
                    std::cerr << "show_codons must be one of 'none', 'ascii', 'horizontal', or 'vertical'" << std::endl;
                    exit (EXIT_FAILURE);
                }
                p.codon_output = *co;
                break;
            case 'd':
                p.dna_template = optarg;
                if (p.dna_template.size() % 3 != 0) {
                    std::cerr << "template_dna must encode a valid orf with length a multiple of 3" << std::endl;
                    exit (EXIT_FAILURE);
                }
                p.aa_template = bio::Nts(p.dna_template);
                p.template_sources.push_back(bio::Cdns(p.dna_template));
                break;
            case 'f':
                p.fw_refs.push_back(optarg);
                break;
            case 'g':
                p.min_umi_group_size = std::strtol(optarg, nullptr, 10);
                if (errno != 0 || p.min_umi_group_size < 1) {
                    std::cerr << "min_umi_grp must be an integer >= 1" << std::endl;
                    exit (EXIT_FAILURE);
                }
                break;
            case 'm':
                p.max_mismatches = std::strtol(optarg, nullptr, 10);
                if (errno != 0 || p.min_overlap < 0) {
                    std::cerr << "max_mismatches must be an integer >= 0" << std::endl;
                    exit (EXIT_FAILURE);
                }
                break;
            case 'n':
                p.number_from = std::strtol(optarg, nullptr, 10);
                if (errno != 0 || p.number_from < 0) {
                    std::cerr << "number_from must be an integer >= 0" << std::endl;
                    exit (EXIT_FAILURE);
                }
            case 'q':
                p.tp_qual_min = optarg[0];
                if (optarg[1] != 0 || p.tp_qual_min < '!' || p.tp_qual_min > '~') {
                    std::cerr << "qual min must be a single ascii character in the interval ['!', '~']" << std::endl;
                    exit (EXIT_FAILURE);
                }
                break;
            case 'r':
                p.rv_refs.push_back(optarg);
                break;
            case 's':
                p.no_header_flag = 1;
                break;
            case 't':
                if (std::strcmp(optarg, "none") == 0) {
                    p.aa_template = "";
                } else {
                    p.aa_template = optarg;
                }
                p.template_sources.push_back(p.aa_template);
                break;
            case 'v':
                p.min_overlap = std::strtol(optarg, nullptr, 10);
                if (errno != 0 || p.min_overlap < 1) {
                    std::cerr << "min_overlap must be an integer >= 1" << std::endl;
                    exit (EXIT_FAILURE);
                }
                break;
            case 'x':
                p.skip_assembly_flag = 1;
                break;
            case '?':
                break;
            default:
                exit (EXIT_FAILURE);
        }
    }

    if (optind ==  argc) {
        std::cerr << "missing positional argument: forward_reads.fastq" << std::endl;
        exit (EXIT_FAILURE);
    }
    p.fw_filename = argv[optind++];

    if (optind == argc) {
        std::cerr << "missing positional argument: reverse_reads.fastq" << std::endl;
        exit (EXIT_FAILURE);
    } 
    p.rv_filename = argv[optind++];

    if (optind != argc) {
        std::cerr << "unexpected positional argument: '" << argv[optind] << "'" << std::endl;
        exit (EXIT_FAILURE);
    }

    if (p.fw_refs.empty()) {
        std::cerr << "at least one reference sequence is required for the forward read (-f, --fw_ref)" << std::endl;
        exit (EXIT_FAILURE);
    }

    if (p.rv_refs.empty()) {
        std::cerr << "at least one reference sequence is required for the reverse read (-r, --rv_ref)" << std::endl;
        exit (EXIT_FAILURE);
    }

    /*
    if (p.template_sources.empty()) {
        std::cerr << "missing template sources (-t, --template) (-d, --template_dna) (--template_db)" << std::endl;
        exit (EXIT_FAILURE);
    }
    */

    if (p.trims.empty()) {
        for (const auto &src : p.template_sources) p.trims.push_back({0,0});
    }

    if (p.trims.size() != p.template_sources.size()) {
        std::cerr << "using -trim requires a separate --trim=x,y for each template source (--template, --template_dna, --template_db)" << std::endl;
        exit (EXIT_FAILURE);
    }

    return p;
}

};