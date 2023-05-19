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

enum class CodonOutput {
    None,
    Ascii,
    Horizontal,
    Vertical 
};

std::optional<CodonOutput> codon_output_from_string(const char *s);

using TemplateSource = std::variant<fs::path, bio::Cdns, bio::Aas>;

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