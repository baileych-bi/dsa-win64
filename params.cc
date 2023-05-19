#include <unordered_map>

#include "params.h"

namespace help {

std::optional<CodonOutput>
codon_output_from_string(const char *s) {
    static const std::unordered_map<std::string, CodonOutput> lookup = {
        {"none",       CodonOutput::None},
        {"ascii",      CodonOutput::Ascii},
        {"horizontal", CodonOutput::Horizontal},
        {"vertical",   CodonOutput::Vertical}
    };

    std::optional<CodonOutput> co;
    std::string lc; for (; *s; ++s) lc.push_back(std::tolower(*s));
    auto ii = lookup.find(lc);
    if (ii != lookup.end()) co = ii->second;
    return co;
}

};