#ifndef BIO_ABS_H_
#define BIO_ABS_H_

#include <filesystem>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>

#include "align.h"

namespace fs = std::filesystem;

namespace bio {

class Aas;
class Cdns;
struct Orf;

class BadTemplateDatabaseParse : public std::runtime_error {
    using runtime_error::runtime_error;
};

class ExcessiveTrimmingError : public std::runtime_error {
    using runtime_error::runtime_error;
};


struct TemplateDatabaseEntry {
    std::string label;
    Cdns        cdns;
    Aas         aas;
};

class TemplateDatabase {
    TemplateDatabase() = default;
    TemplateDatabase(std::vector<TemplateDatabaseEntry> &&targets) : targets_(targets) {}
    std::vector<TemplateDatabaseEntry> targets_;

    int32_t gap_penalty_ = 4;

public:
    static std::shared_ptr<TemplateDatabase> create_empty();
    static std::shared_ptr<TemplateDatabase> from_imgt_fasta(const fs::path &path);

    bool codon_data_available() const { return !targets_.empty() && !targets_.front().cdns.empty(); }
    int32_t gap_penalty() const { return gap_penalty_; }

    void add_entry(const std::string &label, const Cdns &cdns, const Aas &aas);

    void trim(const std::pair<size_t, size_t> &how_much); 

    size_t query(const Cdns &cdns) const;
    size_t query(Cdns::const_iterator lo, Cdns::const_iterator hi) const;

    size_t query_and_align(Cdns::const_iterator lo, Cdns::const_iterator hi, Alignment &result) const;
    size_t query_and_align(const Cdns &cdns, Alignment &result) const {
        return query_and_align(cdns.cbegin(), cdns.cend(), result);
    };

    size_t query_and_align(Aas::const_iterator lo, Aas::const_iterator hi, Alignment &result) const;
    size_t query_and_align(const Aas &aas, Alignment &result) const {
        return query_and_align(aas.cbegin(), aas.cend(), result);
    };

    size_t size() const { return targets_.size(); }

    const TemplateDatabaseEntry &operator[](size_t i) const { assert(i != 0); return targets_[i-1]; }

    const std::string &get_label(size_t i) const { assert(i != 0); return targets_[i-1].label; }
    const Cdns &get_codons(size_t i) const { assert(i != 0); return targets_[i-1].cdns; }
    const Aas  &get_aas   (size_t i) const { assert(i != 0); return targets_[i-1].aas ; }

    using const_iterator = decltype(targets_.cbegin());

    const_iterator begin()  const { return targets_.cbegin(); }
    const_iterator end()    const { return targets_.cend();   }
    const_iterator cbegin() const { return targets_.cbegin(); }
    const_iterator cend()   const { return targets_.cend();   }

    static constexpr const size_t NOT_FOUND = 0;
};

};

#endif