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

#ifndef BIO_MAINFUNCTIONS_H_
#define BIO_MAINFUNCTIONS_H_
#include <ostream>
#include <unordered_map>

#include "align.h"
#include "defines.h"
#include "help.h"
#include "io.h"
#include "params.h"
#include "umi.h"

namespace bio {

class TemplateDatabase;

/** Class similar to Python's collections.Counter */
template<typename T, typename Hash=std::hash<T>, typename KeyEqual=std::equal_to<T>>
struct Counter {
  using iterator       = typename std::unordered_map<T, size_t, Hash, KeyEqual>::const_iterator;
  using const_iterator = typename std::unordered_map<T, size_t, Hash, KeyEqual>::const_iterator;

  iterator begin()  const { return counts_.cbegin(); }
  iterator end()    const { return counts_.cend();   }
  iterator cbegin() const { return counts_.cbegin(); }
  iterator cend()   const { return counts_.cend();   }

  void push_back(const T &k) { counts_[k] += 1; total_ += 1; }
  size_t operator[](const T &k) const { auto ii = counts_.find(k); return (ii == counts_.end()) ? 0 : ii->second; }
  size_t total() const { return total_; }

private:
  size_t total_ = 0;
  std::unordered_map<T, size_t, Hash, KeyEqual> counts_;
};

/** Maintains a count of sequences dropped from analysis for various reasons. */
struct ParseLog {
    size_t filter_invalid_chars            = 0;
    size_t filter_no_fw_umi                = 0;
    size_t filter_no_rv_umi                = 0;
    size_t filter_could_not_assemble       = 0;
    size_t filter_umi_group_size_too_small = 0;
    size_t filter_duplicate_umi            = 0;
    size_t filter_premature_stop_codon     = 0;
    size_t filter_split_failed             = 0;
    size_t filter_no_matching_template     = 0;
    size_t filter_bad_alignment            = 0;

    ParseLog operator+(const ParseLog &l) {
        ParseLog sum = *this;
        sum.filter_invalid_chars += l.filter_invalid_chars;
        sum.filter_no_fw_umi += l.filter_no_fw_umi;
        sum.filter_no_rv_umi += l.filter_no_rv_umi;
        sum.filter_could_not_assemble += l.filter_could_not_assemble;
        sum.filter_umi_group_size_too_small += l.filter_umi_group_size_too_small;
        sum.filter_duplicate_umi += l.filter_duplicate_umi;
        sum.filter_premature_stop_codon += l.filter_premature_stop_codon;
        sum.filter_split_failed += l.filter_split_failed;
        sum.filter_no_matching_template += l.filter_no_matching_template;
        sum.filter_bad_alignment += l.filter_bad_alignment;
        return sum;
    }
};

/** Enumerated template used for alignment. Multiple templates are used with --split
    that breaks into different parts (e.g. of an antibody) to align to different
    templates. */
struct AlignmentTemplate {
    size_t id = 0;
    std::vector<std::string> labels;
    Aas aas;
    Cdns cdns;

    std::string label(const std::string &delim=" / ") const;
};

/** GroupAlignment represents an alignment of a UMI group to a particular template.
* Populating a vector of GroupAlignments is basically the goal of dsa.
*/
struct GroupAlignment {
    size_t umi_group_size = 0;                //< Number of PCR reads in this umi group
    std::shared_ptr<AlignmentTemplate> templ; //< The AlignmentTemplate
    std::string barcode;                      //< The UMI nucleotide sequence
    std::string alignment;                    //< The alignment itself (as amino acids)
    std::string cdns;                         //< The alignment itself (as codons)

    /** GroupAlignments can be concatenated. This is useful if, for example,
    * The left and right halves of the reads are aligned to different templates
    * but we want to reconstruct the entire alignment in a single string.
    */
    GroupAlignment &operator+=(const GroupAlignment &g) {
      alignment += g.alignment;
      cdns      += g.cdns;
      return *this;
    }
};

/** Store the counts of coding and silet mutations at each codon position. */
struct MutationCount {
    MutationCount() = default;
    MutationCount(size_t cols);

    std::vector<unsigned> synonymous;    //< Position-wise counts of silent mutations
    std::vector<unsigned> nonsynonymous; //< Position-wise counts of coding mutations
    std::vector<unsigned> total;         //< Total number of sequences examined at each position

    MutationCount &operator=(const MutationCount &) = default;
    MutationCount  operator+(const MutationCount &) const;
};

/**
  * Parse a memory mapped .fastq file.
  *
  * Reads will be filtered if the contain non ATGC characters (i.e. Ns)
  * or if there is a mismatch between the sequence length and the fastq
  * quality length.
  *
  * @params mapping a memory mapped fastq file
  * @return the unpaired reads and quality data
  */
std::vector<Read>
extract_read_data(const ConstMapping &mapping);

/**
  * Remove poor quality sequeces from data.
  *
  * Bases at 3' read ends are removed if they fall below params.tq_qual_min.
  * The reference sequence itself will be trimmed from the read and the UMI barcode extracted. 
  * Reads fail QC if the fw or rv reference sequence/UMI cannot be identified.
  *
  * @param fw the unpaired forward reads
  * @param rv the unpaired reverse reads
  * @param fwex a UMIExtractor initialized with the forawrd reference sequence
  * @param rvex a UMIExtractor initialized with the reverse reference sequence
  * @param params run options from command line arguments
  * @param log ParseLog to store counts of reads that fail QC for one reason or another
  *
  * @return pairs of reads (not yet assembled) for which both fw and rv passed QC
  */
std::vector<ReadPair>
qc_reads(
    std::vector<Read> &&fw,
    std::vector<Read> &&rv,
    const std::vector<UMIExtractor> &fwexs,
    const std::vector<UMIExtractor> &rvexs,
    const help::Params &params,
    ParseLog &log);

/**
  * Assemble paired-end reads.
  *
  * Will assemble reads by aligning the 3' ends of the fw read and the reverse complement
  * of the rv read. Requires a minimum number of overlapping bases (params.min_overlap)
  * but permits some number of mismatces (params.max_mismatches) within the region of
  * overlap.
  *
  * @param pairs the unpaired reads
  * @param params run options from command line arguments
  * @param log ParseLog to store counts of reads that to align
  *
  * @return pairs of reads (not yet assembled) for which assembly was successful
  */
std::vector<Read>
assemble_reads(
    std::vector<ReadPair> &&pairs,
    const help::Params &params,
    ParseLog &log);

/**
  * Build a conensus nucleotide sequence from a group of reads.
  *
  * Two means of building the consensus are available based on whether or
  * not all reads are expected to be the same length (paired reads) or
  * different lengths varying in their 3' tails (unpaired reads).
  * The ragged_ends flag should be set to false for paired reads or true for
  * unpaired reads.
  * <br/>
  * UMI groups with fewer than params.min_umi_group_size members will be discarded.
  * <br/>
  * For unpaired reads a consensus will be created with length such that a 
  * the consensus nucleotide at any given position represents at least
  * params.min_umi_group_size sequences.
  * <br/>
  * For paired reads, the consensus sequence length will be the modal length
  * of the sequences in the group.
  * <br/>
  * For each position in the consensus, the consensus nucleotide will be the modal
  * nucleotide at that position in the inputs. If there is no single mode (e.g. for a group
  * of two inputs), fastq quality will be used to break the tie.
  *
  * @param reads the reads, paired or unpaired
  * @param params run options from command line arguments
  * @param log ParseLog to store counts of reads "collapsed"
  * @param ragged_ends set true when reads are expected to vary in length (i.e. unpaired reads)
  *
  * @return the consensus sequences for each UMI group
  */
std::vector<Read>
umi_collapse(
    std::vector<Read> &&reads,
    const help::Params &params,
    ParseLog &log,
    bool ragged_ends);

/**
  * Translate ORFs from nucleotide data.
  *
  * Reads should be in frame such that the first nucleotide of a
  * forward read is the first nucleotide of the first codon in the ORF.
  * The first nucleotide of a reverse read should be the last
  * nucleotide of the last codon in the ORF.
  *
  * @param reads the nucleotide read data (paired or unpaired)
  * @param params run options from command line arguments
  * @param log ParseLog to store counts of ORFs containing stop codons
  * @param reverse_complement set true to translate the reverse complement strand (i.e. for reverse reads)
  *
  * @return ORFs without stop codons
  */
std::vector<Orf>
translate_and_filter_ptcs(std::vector<Read> &&reads, 
                          const help::Params &params,
                          ParseLog &log,
                          bool reverse_complement);


/**
  * Align ORFs to an amino acid template.
  *
  * Output is a gapped string of ASCII amino acids.
  *
  * Low scoring alignments (those scoring below params.min_alignment_score) 
  * will be filtered.
  * @param orfs the ORFs to align to the template params.aa_template
  * @param params run options from command line arguments
  * @param log ParseLog to store counts of ORFs containing stop codons
  *
  * @return the alignments that pass QC
  */
std::vector<GroupAlignment>
align_to_template(std::vector<Orf> &&orfs,
                  const std::vector<AlignmentTemplate> &templates,
                  const help::Params &params,
                  ParseLog &log,
                  bool ragged_ends);

enum class AlignmentMethod {
    ByCdns,
    ByAas
};

std::vector<GroupAlignment>
align_to_multiple_templates(vecvec<Orf> &&orfs,
                   const std::vector<std::shared_ptr<const TemplateDatabase>> &dbs,
                   const help::Params &params,
                   ParseLog &log,
                   bool ragged_ends=false);

/**
  * Split ORFs according to params.split_template_regex.
  *
  * @param orfs the orfs to split
  * @param params run options from command line arguments
  * @param log ParseLog to store counts of ORFs that could not be split
  *
  * @return a vector of vector of Orf where each element of the inner vectors is one of the split orfs
  */
vecvec<Orf>
split_orfs(std::vector<Orf> &&orfs,
           const help::Params &params,
           ParseLog &log);

/**
  * Set the template_id of each Orf to the index of the closest template in the corresponding database.
  *
  * Note, this function assumes split Orfs and multiple template databases.
  * For example, orfs[i][1] will be matched against the templates in dbs[1]
  * and orfs[j][2] will be matched against the templates in dbs[2].
  *
  * @param orfs 2D array of Orf; 2nd dimension must be dbs.size()
  * @param dbs the template databases
  * @param params run options from command line arguments
  * @param log ParseLog to store counts of ORFs where no matching template could be found
  *
  * @return 2D array of Orf with template ids assigned, less any where template could not be matched
  */
vecvec<Orf>
assign_templates(std::vector<std::vector<Orf>> &&orfs,
                 std::vector<std::shared_ptr<const TemplateDatabase>> &dbs,
                 const help::Params &params,
                 ParseLog &log);

}; //namespace bio

#endif