/*
Copyright 2024, The Broad Institute of MIT and Harvard

Original Author: Charles C Bailey

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef HELP_PARAMS_H_
#define HELP_PARAMS_H_

#include <string>
#include <optional>

#include "aa.h"
#include "dna.h"

#include <filesystem>
#include <optional>
#include <regex>
#include <variant>
#include <vector>

namespace fs = std::filesystem;

namespace help {

/** Supported means of printing nucleotide data. */
enum class CodonOutput {
    None,        //< Don't print nucleotide/codon data
    Ascii,       //< Print single char representation of codons
    Horizontal,  //< Print dna single-line sequences rather than codon sequences
    Vertical     //< Print three-row vertically-aligned dna sequences 
};

/** Maybe get a CodonOuput enum value from a string. */
std::optional<CodonOutput>
codon_output_from_string(const char *s);

/** Alignment templates can be dna sequences (packed as Cdns),
* amino acid sequences (Aas), or special files containing lists
* of sequences (std::path to a .fasta file)
*/
using TemplateSource = std::variant<fs::path, bio::Cdns, bio::Aas>;

/** Params holds the various user-defined analysis settings determined by the command line
* arguments.
*/
struct Params {
    bool split_template_requested() const { return split_template_regex.mark_count() > 0; }

    std::string fw_filename;
    std::string rv_filename; 
    std::vector<std::string> fw_refs;
    std::vector<std::string> rv_refs;

    bio::Aas aa_template;
    bio::Nts dna_template;

    std::string split_template_string;
    std::regex  split_template_regex;
    std::vector<TemplateSource> template_sources;
    std::vector<std::pair<size_t, size_t>> trims;

    int no_header_flag     = 0;
    int skip_assembly_flag = 0;
    int allow_ptcs_flag    = 0;
    int separate_cdr3_flag = 0;

    float min_alignment_score = 0.8f;
    char  tp_qual_min         = 'A';
    long  min_umi_group_size  = 1;
    long  min_overlap         = 9;
    long  max_mismatches      = 0;
    long  number_from         = 1;

    CodonOutput codon_output = CodonOutput::None;
};

}; //namespace help

#endif //HELP_PARAMS_H_