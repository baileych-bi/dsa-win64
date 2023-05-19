#include "umi.h"

namespace bio {

UMIExtractor::UMIExtractor(const std::string &sequence) 
    : sequence_(sequence)
    , regex_() {
    static const std::string valid_chars = "ACGTNn";
    for (char &c : sequence_) {
        if (c != 'n' && c != 'N') c = std::toupper(c);
        if (valid_chars.find(c) == std::string::npos) throw std::runtime_error("Invalid UMI.");
    }
    bool capture = false;
    for (char c : sequence_) {
        if (capture) {
            if (c == 'n') {
                pattern_.push_back('.');
            } else {
                pattern_.push_back(')');
                pattern_.push_back(c == 'N' ? '.' : std::toupper(c));
                capture = false;
            }
        } else {
            if (c == 'n') {
                pattern_.push_back('(');
                pattern_.push_back('.');
                capture = true;
            } else {
                pattern_.push_back(c == 'N' ? '.' : std::toupper(c));
            }
        }
    }
    if (capture) pattern_.push_back(')');
    regex_ = std::regex(pattern_, std::regex::icase);
}

ExtractedUMI
UMIExtractor::operator()(const Nt *from_nt, const Nt *to_nt) const {
    ExtractedUMI result;
    std::cmatch match;
    const char *from = reinterpret_cast<const char *>(from_nt);
    const char *to   = reinterpret_cast<const char *>(  to_nt);
    if (!std::regex_search(from, to, match, regex_)) return result;

    for (size_t i=1; i<match.size();++i) {
        result.barcode += match[i];
    }

    result.from   = match.position(0);
    result.length = match.length(0);
    return result;
}

ExtractedUMI
UMIExtractor::operator()(Nts::const_iterator from, Nts::const_iterator to) const {
    return (*this)(from.operator->(), to.operator->());
}

}; //namespace bio