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

#include <fstream>
#include <regex>
#include <stdlib.h>

#include "abs.h"

#ifdef DSA_TARGET_WIN64
namespace detail {
std::string
    dsa_get_env(const char* what) {
    char* home = nullptr;
    size_t buf_size = 0;
    _dupenv_s(&home, &buf_size, what);
    std::string result(home);
    free(home);
    return result;
}
}; //naemspace detail
#elif defined(DSA_TARGET_LINUX)
namespace detail {
std::string
    dsa_get_env(const char* what) {
    return std::string(std::getenv(what));
}
};
#endif
namespace bio {

std::vector<std::string>
split(const std::string &str, const std::string &delim) {
    std::vector<std::string> tokens;
    if (str.empty()) return tokens;

    std::string::const_iterator lo = str.cbegin(), hi;
    while ((hi = std::search(lo, str.cend(), delim.cbegin(), delim.cend())) != str.cend()) {
        tokens.push_back(std::string(lo, hi));
        lo = hi + delim.size();
    }
    tokens.push_back(std::string(lo, hi));
    return tokens;
}

std::shared_ptr<TemplateDatabase>
TemplateDatabase::create_empty() {
    return std::shared_ptr<TemplateDatabase>(new TemplateDatabase);
}

std::shared_ptr<TemplateDatabase>
TemplateDatabase::from_imgt_fasta(const fs::path &path) {
    fs::path resolved = path;

    if (path.begin() != path.end() && *path.begin() == "~") {
        std::string home = detail::dsa_get_env("HOME");
        if (!home.empty()) {
            resolved = fs::path(home);
            auto ii = path.begin(); ++ii;
            for (; ii != path.end(); ++ii) resolved /= *ii;
        }
    }

    std::ifstream ifs(resolved);
    if (!ifs) {
        throw BadTemplateDatabaseParse(std::string("File '") + path.string() + "' could not be opened for reading.");
    }

    std::shared_ptr<TemplateDatabase> ptr;
    static const std::regex allelic_variant(R"-(\*0[2-9]$)-"); //when parsing IMGT-style .fasta we will skip allelic variants of things (e.g. FOO*02, FOO*03, etc.)
    std::smatch match;

    std::string line;
    std::vector<TemplateDatabaseEntry> records;

    size_t line_no = 0;

    std::string label;
    Nts nts;

    while (std::getline(ifs, line)) {
        ++line_no;
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r')) line.pop_back(); //strip LF/CRLF
        if (line.empty()) continue; //skip blank lines

        if (line[0] == '>') {
            if (!label.empty()) {
                if (!std::regex_search(label, match, allelic_variant)) {
                    records.push_back({std::move(label), Cdns(std::move(nts)), Aas("")});
                }
                label.clear();
                nts.clear();
            }

            std::vector<std::string> tokens = split(line, "|");
            if (tokens.size() == 1) { //a single token means we use the whole header (minus the '>') as the label
                if (tokens[0].size() < 2) {
                    throw BadTemplateDatabaseParse(
                        std::string("Bad header '") + line + "' on line " + std::to_string(line_no) + ": identifier field is empty."
                    );                    
                }
                label = tokens[0].substr(1); //remove the '>' from the front
            } else { //otherwise we assume we're dealing with an IMGT-style header; tokens[1] should contain the identifier
                if (tokens.size() < 2) {
                    throw BadTemplateDatabaseParse(
                        std::string("Bad header '") + line + "' on line " + std::to_string(line_no) + ": not enough fields."
                    );
                }

                if (tokens[1].empty()) {
                    throw BadTemplateDatabaseParse(
                        std::string("Bad header '") + line + "' on line " + std::to_string(line_no) + ": identifier field is empty."
                    );
                }
                label = std::move(tokens[1]);
            }
        } else if (label.empty()) {
            throw BadTemplateDatabaseParse(
                std::string("Unexpected sequence data '") + line + "' on line " + std::to_string(line_no)
            );
        } else {
            nts += Nts(line);
        }
    }

    if (!label.empty()) {
        if (!std::regex_search(label, match, allelic_variant)) {
            records.push_back({std::move(label), Cdns(std::move(nts)), Aas("")});
        }
        label.clear();
        nts.clear();
    }

    if (records.empty()) {
        throw BadTemplateDatabaseParse("No fasta records founnd");
    }

    for (TemplateDatabaseEntry &record : records) record.aas = Aas(record.cdns);

    return std::shared_ptr<TemplateDatabase>(new TemplateDatabase(std::move(records)));
}

void
TemplateDatabase::add_entry(const std::string &label, const Cdns &cdns, const Aas &aas) {
    targets_.push_back({label, cdns, aas});
}

void
TemplateDatabase::trim(const std::pair<size_t, size_t> &how_much) {
    size_t total = how_much.first + how_much.second;
    for (TemplateDatabaseEntry &entry : targets_) {
        if (total >= entry.aas.size()) throw ExcessiveTrimmingError(
            std::string("Cannot trim ") + 
            std::to_string(total) + 
            " amino acids from '" + 
            entry.label + 
            "', a template of only " + 
            std::to_string(entry.aas.size()) +
            "amino acids."
        );
        entry.aas.exo(how_much.first, how_much.second);
        if (!entry.cdns.empty()) entry.cdns.exo(how_much.first, how_much.second);
    }
}


size_t
TemplateDatabase::query_and_align(Cdns::const_iterator lo, Cdns::const_iterator hi, Alignment &result) const {
    size_t best_tpl = NOT_FOUND;
    result.clear();

    result.score = std::numeric_limits<decltype(result.score)>::min();

    Alignment current_alignment;
    for (size_t i=0; i<targets_.size(); ++i) {
        current_alignment.clear();

        nw_align<Cdn>(lo, hi, 
            targets_[i].cdns.cbegin(), 
            targets_[i].cdns.cend(), 
            CDNSUBS, gap_penalty_, current_alignment, true);

        if (current_alignment.score > result.score) {
            best_tpl = i+1;
            result = std::move(current_alignment);
        }
    }

    return best_tpl;
}

size_t
TemplateDatabase::query_and_align(Aas::const_iterator lo, Aas::const_iterator hi, Alignment &result) const {
    size_t best_tpl = NOT_FOUND;
    result.clear();

    result.score = std::numeric_limits<decltype(result.score)>::min();

    Alignment current_alignment;
    for (size_t i=0; i<targets_.size(); ++i) {
        current_alignment.clear();

        nw_align<Aa>(lo, hi, 
            targets_[i].aas.cbegin(), 
            targets_[i].aas.cend(), 
            BLOSUM62, gap_penalty_, current_alignment, true);

        if (current_alignment.score > result.score) {
            best_tpl = i+1;
            result = std::move(current_alignment);
        }
    }

    return best_tpl;
}

size_t
TemplateDatabase::query(Cdns::const_iterator lo, Cdns::const_iterator hi) const {
    Alignment alignment;
    size_t max_tpl = NOT_FOUND;
    float max_score = static_cast<float>(std::numeric_limits<decltype(alignment.score)>::min());
    for (size_t i=0; i<targets_.size(); ++i) {
        nw_align<Cdn>(lo, hi, targets_[i].cdns.cbegin(), targets_[i].cdns.cend(), CDNSUBS, gap_penalty_, alignment, true);
        if (alignment.score > max_score) {
            max_tpl = i+1;
            max_score = static_cast<float>(alignment.score);
        }
    }
    return max_tpl;
}

};