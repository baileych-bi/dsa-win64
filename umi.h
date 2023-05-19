#ifndef BIO_UMIEXTRACTOR_
#define BIO_UMIEXTRACTOR_

#include <iostream>
#include <regex>
#include <vector>

#include "dna.h"

namespace bio {

/** Stores the result of an attempt to match a reference sequence to an ASCII string */
struct ExtractedUMI {
    std::string barcode;                 ///< the barcode (if the reference specifices UMI capture)
    size_t from   = std::string::npos;   ///< the index of base 1 of the reference match
    size_t length = 0;                   ///< the length of the match

    bool valid()   const { return length != 0; } ///< true if reference was found
    bool invalid() const { return length == 0; } ///< true if reference not found
};

/** Class for extracting UMI barcodes from ascii sequence data.
  * Takes input as a string containing ATGC bases and Ns.  Contiguous
  * strecthes of Ns are converted to regex capture groups.  Running
  * the extractor by means of the () operator creates a concatenated string
  * of nucleotides found at the N positions or an empty string if UMI was
  * not found.
  */
class UMIExtractor {
    std::regex regex_;
    std::string pattern_;
    std::string sequence_;
public:
    UMIExtractor() = default;
    UMIExtractor(UMIExtractor &&) = default;
    /** Construct a UMIExtractor that recognizes a given ASCII nucleotide sequence.
      *
      * The sequence has certain formatting requirements. Capital ATGC perform literal
      * matching. Capital 'N' is a wildcard. Lowercase 'n' is a wildcard but also
      * records the base it matches as part of the UMI barcode. If sequence contains
      * invalid characters (i.e. anything other than ATGCNn) this constructor will throw
      * std::runtime_error
      *
      * @param sequence an ASCII DNA reference sequence formatted as described above.
      */
    explicit UMIExtractor(const std::string &sequence); //can throw std::regex_error
    UMIExtractor &operator=(const UMIExtractor &) = default;
    UMIExtractor &operator=(UMIExtractor &&) = default;

    /** Look for reference sequence in [from, to)
      * 
      * The given interval [to, from) will be searched for the reference sequence
      * supplied at construction. If found, the ExtractedUMI returned will be valid
      * and contain the UMI barcode (if any).
      *
      * @param from the first char of the interval to search.
      * @param to   the last+1 char of the interval to search.
      *
      * @return an ExtractedUMI representing the result of the search
    */
    ExtractedUMI operator()(const Nt *from, const Nt *to) const;

    /** Look for reference sequennce in [begin, end) */ 
    ExtractedUMI operator()(Nts::const_iterator begin, Nts::const_iterator end) const;

    bool empty() const { return pattern_.empty(); }

    const std::string &sequence() const { return sequence_; } ///< The original reference sequence as supplied to the constructor.
    const std::string & pattern() const { return  pattern_; } ///< The regular expression string created from the reference sequence.
    const std::regex  &   regex() const { return    regex_; } ///< The regular exrpession object used for matching and extracting UMIs.
};

}; //namespace bio

#endif //#ifndef BIO_UMIEXTRACTOR_
