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

/** \mainpage Deep Sequencing Analysis
  * \section Overview
  * dsa analyzes paired-end reads that encode a protein sequence, aligns the reads to user-supplied
  * template(s), and outputs those alignments and various statistics about the mutations.<br/>
  * The current implementation relies on AVX2 instructions and requires somewhat recent x86 CPUs.
  */
#include <immintrin.h>

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "defines.h"

#include "aa.h"
#include "abs.h"
#include "align.h"
#include "cdn.h"
#include "dna.h"
#include "help.h"
#include "io.h"
#include "mainfunctions.h"
#include "parallelism.h"
#include "params.h"
#include "polymer.h"
#include "umi.h"
#include "tests.h"

namespace fs = std::filesystem;
using namespace bio;
using help::Params;

int
main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Deep Sequencing Analysis version " << VERSION_STRING << ": "
                  << "run dsa --help for instructions." << std::endl;
        exit (EXIT_FAILURE);
    }

    if (argc == 2 && !std::strcmp("test", argv[1])) {
        try {
            test::run_all();
        } catch (test::test_failed_error &ex) {
            std::cerr << "test failed:" << std::endl;
            std::cerr << ex.what() << std::endl;
            throw ex;
        }
        std::cout << "All tests successful." << std::endl;
        exit (EXIT_SUCCESS);
    }

    const help::Params p = help::parse_argv(argc, argv);

    const std::string VERSION = VERSION_STRING;

    //check the command line arguments for bad values and mutually exclusive options
    if (p.min_overlap < p.max_mismatches) {
        std::cerr << "max_mismatches must be less than min_overlap" << std::endl;
        exit (EXIT_FAILURE);
    }

    //std::cerr << "p.template_sources.size()==" << p.template_sources.size() << std::endl;
    if (p.skip_assembly_flag && p.template_sources.size() > 1) {
        std::cerr << "skipping assembly (i.e. -x, --skip_assembly) is incompatible with "
                  << "split templates and multiple template alignment" << std::endl;
        exit (EXIT_FAILURE);
    }

    std::vector<std::shared_ptr<const TemplateDatabase>> template_dbs;

    if (p.split_template_regex.mark_count() != 0 &&
        p.split_template_regex.mark_count() != p.template_sources.size()) {
        std::cerr << "when splitting reads for multi-template alignment (--split), "
                  << "a template source (--template, --template_dna, --template_db) must be provided for each capturing "
                  << "subgroup of the regular expression (see --help_split)" << std::endl;
        exit (EXIT_FAILURE);
    } else {
        for (size_t i=0; i<p.template_sources.size(); ++i) {
            const help::TemplateSource &source = p.template_sources[i];
            std::shared_ptr<TemplateDatabase> db;
            if (std::holds_alternative<fs::path>(source)) {
                const fs::path &filename = std::get<fs::path>(source);
                try {
                    db = TemplateDatabase::from_imgt_fasta(filename);
                } catch (const BadTemplateDatabaseParse &ex) {
                    std::cerr << "could not parse '" << filename.string() << "' as a template database: " << std::endl
                            << "Error: " << ex.what() << std::endl
                            << "databases should be .fasta files of in-frame nucleotides with IGMT-style headers (see --help_split)" << std::endl;
                    exit (EXIT_FAILURE);
                }
            } else if (std::holds_alternative<Cdns>(source)) {
                db = TemplateDatabase::create_empty();
                db->add_entry("user_defined_cdns", std::get<Cdns>(source), Aas(std::get<Cdns>(source)));
            } else if (std::holds_alternative<Aas>(source)) {
                Aas aas = std::get<Aas>(source);
                if (aas.empty()) {
                    db = nullptr;
                } else {
                    db = TemplateDatabase::create_empty();
                    db->add_entry("user_defined_aas", Cdns(), std::get<Aas>(source));
                }
            } else {
                std::cerr << "Unkown template source." << std::endl;
                exit (EXIT_FAILURE);
            }

            try {
                if (db) db->trim(p.trims[i]);
            } catch (ExcessiveTrimmingError &ex) {
                std::cerr << ex.what() << std::endl;
                exit (EXIT_FAILURE);
            }

            template_dbs.push_back(db);
        }
    }

    if (p.skip_assembly_flag && (template_dbs.size() * template_dbs.front()->size() > 1)) {
        std::cerr << "skipping assembly (i.e. -x, --skip_assembly) is incompatible with "
                  << "split templates and multiple template alignment" << std::endl;
        exit (EXIT_FAILURE);
    }

    std::vector<UMIExtractor> fwexs;
    for (const std::string &ref : p.fw_refs) {
        try {
            fwexs.push_back(UMIExtractor(ref));
        } catch (std::exception &) {
            std::cerr << "fw_ref '" << ref << "' is not a valid reference sequence (see --help)" << std::endl;
            exit (EXIT_FAILURE);
        }
    }

    std::vector<UMIExtractor> rvexs;
    for (const std::string &ref : p.rv_refs) {
        try {
            rvexs.push_back(UMIExtractor(ref));
        } catch (std::exception &) {
            std::cerr << "rv_ref '" << ref << "' is not a valid reference sequence (see --help)" << std::endl;
            exit (EXIT_FAILURE);
        }
    }

    //Filling out 'alignments' is the ultimate goal of our program.
    //These GroupAlignments represent the Needleman-Wunsch alignments
    //of translated paired or unpaired read data to the user-supplied
    //template(s)
    //If no templates are given, dsa will just return UMI-grouped
    //lists of sequences found between the references
    std::vector<GroupAlignment> alignments;

    auto clock_start = std::chrono::high_resolution_clock::now();

    ParseLog log;
    std::vector<Read> fwreads, rvreads; //at first we hold the reads from the two fastq files separately

    //parse the fastq files into Read data structures
    try {
        ConstMapping fwmap = ConstMapping::map(p.fw_filename);
        fwreads = extract_read_data(fwmap);
        fwmap.unmap();
    } catch (std::exception &) {
        std::cerr << "error parsing '" << p.fw_filename << "'" << std::endl;
        exit (EXIT_FAILURE);
    }

    try {
        ConstMapping rvmap = ConstMapping::map(p.rv_filename);
        rvreads = extract_read_data(rvmap);
        rvmap.unmap();
    } catch (std::exception &) {
        std::cerr << "error parsing '" << p.rv_filename << "'" << std::endl;
        exit (EXIT_FAILURE);
    }

    //make sure we got the same number of forward and reverse reads
    if (fwreads.size() != rvreads.size()) {
        std::cerr << "read count disagreement between " << p.fw_filename << " and " << p.rv_filename << std::endl;
        exit (EXIT_FAILURE);
    }

    const size_t total_reads = fwreads.size();

    /* uncomment to limit reads for testing purposes
    const size_t max_reads = 10000;
    fwreads.resize(max_reads);
    rvreads.resize(max_reads);
    */

    //Perform qc and get back read pairs
    //QC includes locating the primers, extracting the UMI,
    //and trimming low-quality bases from the 3' ends of the reads.
    std::vector<ReadPair> qcd_pairs = qc_reads(
        std::move(fwreads),
        std::move(rvreads),
        fwexs, rvexs, p, log);

    fwreads.clear();
    rvreads.clear();

    //Sometimes data are low enough quality that the 3' ends are too hard to
    //align or the PCR template may be too long to sequence. In these cases,
    //we can skip assembling the read pairs and process them anyway.
    if (p.skip_assembly_flag) {
        for (ReadPair &rp : qcd_pairs) {
            rp.rv.barcode = rp.fw.barcode;
            fwreads.push_back(std::move(rp.fw));
            rvreads.push_back(std::move(rp.rv));
        }
        qcd_pairs.clear(); qcd_pairs.shrink_to_fit();
    
        //UMI collapse gives us consensus sequences for the UMI groups
        fwreads = umi_collapse(std::move(fwreads), p, log, true);
        rvreads = umi_collapse(std::move(rvreads), p, log, true);

        //translate
        std::vector<Orf> nterm = translate_and_filter_ptcs(std::move(fwreads), p, log, false);
        fwreads.clear(); fwreads.shrink_to_fit();

        //We don't support splitting unpaired reads so this step just reorganizes
        //the data structures so they can be passed to the template alignment functions
        vecvec<Orf> nsplits = split_orfs(std::move(nterm), p, log);
        nterm.clear(); nterm.shrink_to_fit();

        std::vector<Orf> cterm = translate_and_filter_ptcs(std::move(rvreads), p, log, true);
        rvreads.clear(); rvreads.shrink_to_fit();

        vecvec<Orf> csplits = split_orfs(std::move(cterm), p, log);
        cterm.clear(); cterm.shrink_to_fit();

        //Align our reads to the template; 5' and 3' are aligned separately
        std::vector<GroupAlignment> fwaln = align_to_multiple_templates(
            std::move(nsplits),
            template_dbs,
            p,
            log,
            true);

        std::vector<GroupAlignment> rvaln = align_to_multiple_templates(
            std::move(csplits),
            template_dbs,
            p,
            log,
            true);

        alignments.reserve(fwaln.size() + rvaln.size());

        //Now we collate alignments such that the fw alignment of a read pair always
        //preceeds the rv alignment. Some pairs might only have one viable alignment
        //at this point (for example, one Orf is fine but the other has a PTC)
        //so those will be listed at the end of the alignments section. We do this by
        //sorting fw and rv alignments by barcode then interleaving them into the alignments
        //vector while storing the unpaired alignments in the unpaired vector.
        std::vector<GroupAlignment> unpaired; //to store the unpaired alignments

        //reverse-sort by barcode
        struct by_barcode { 
            bool operator()(const GroupAlignment &a, const GroupAlignment &b) const {
                return a.barcode > b.barcode;
            }
        };
        std::sort(fwaln.begin(), fwaln.end(), by_barcode());
        std::sort(rvaln.begin(), rvaln.end(), by_barcode());

        //collate
        auto ff = fwaln.cbegin(), rr = rvaln.cbegin();
        while (!fwaln.empty() && !rvaln.empty()) {
            const int cmp = fwaln.back().barcode.compare(rvaln.back().barcode);
            if (cmp == 0) {
                alignments.push_back(std::move(fwaln.back()));
                fwaln.pop_back();
                alignments.push_back(std::move(rvaln.back()));
                rvaln.pop_back();
            } else if (cmp < 0) {
                unpaired.push_back(std::move(fwaln.back()));
                fwaln.pop_back();
            } else {
                unpaired.push_back(std::move(rvaln.back()));
                rvaln.pop_back();
            }
        }

        //combine with unpaired and left-over alignments
        alignments.insert(alignments.cend(),
                          std::make_move_iterator(unpaired.begin()),
                          std::make_move_iterator(unpaired.end()));

        alignments.insert(alignments.cend(),
                          std::make_move_iterator(fwaln.rbegin()),
                          std::make_move_iterator(fwaln.rend()));

        alignments.insert(alignments.cend(),
                          std::make_move_iterator(rvaln.rbegin()),
                          std::make_move_iterator(rvaln.rend()));
    } else { //assembling the read ends makes life much easier
        fwreads.shrink_to_fit(); rvreads.shrink_to_fit();

        std::vector<Read> reads = assemble_reads(std::move(qcd_pairs), p, log);
                          reads = umi_collapse(std::move(reads), p, log, false);
        std::vector<Orf>  orfs  = translate_and_filter_ptcs(std::move(reads), p, log, false);

        //Note that the split/multitemplate code path and the sigle template code
        //path are the same. If there is no regex for splitting, split_orfs just turns
        //the 1D orfs vector, shape=(orfs.size(), ) into a 2D vector of shape=(orfs.size(), 1)
        vecvec<Orf> splits = split_orfs(
            std::move(orfs),
            p,
            log
        );
        orfs.clear(); orfs.shrink_to_fit();

        //Again, the single- and multi-template code paths are the same.
        //Single templates are just folded into single-entry template
        //databases.
        const size_t n_splits = splits.size();
        alignments = align_to_multiple_templates(
            std::move(splits),
            template_dbs,
            p,
            log
        );
    }

    std::vector<std::shared_ptr<AlignmentTemplate>> templates;
    std::vector<Matrix<float>> substitution_matrices;

    //If we have more than one template, we sort the alignments by template id
    //so that the output has similar sequences adjacent to one another
    std::sort(alignments.begin(), 
        alignments.end(), 
        [](const GroupAlignment &a, const GroupAlignment &b)->bool{
            if ( a.templ ==  b.templ) return false;
            if (!a.templ &&  b.templ) return true;
            if ( a.templ && !b.templ) return false;
            return a.templ->id < b.templ->id;
        }
    );

    //Alignments are sorted by tempate_id where template_id is the index into templates
    //we now iterate over the templates and process only alignments with the corresponding template_id
    //the range of alignments will be stored in [lo, hi)
    std::vector<GroupAlignment>::iterator lo=alignments.begin(), hi=alignments.begin();
    while (hi != alignments.end()) {
        //find the batch of alignments whose template_id matches the current template, i;
        //these will be in the range [lo, hi)
        lo = hi;

        //substitutions only make sense in the context of a template
        //so we skip untemplated stuff.
        if (lo->templ == nullptr) {
            ++hi;
            continue;
        }

        size_t i = lo->templ->id;
        const Aas &templ = lo->templ->aas;
        templates.push_back(lo->templ);

        hi = std::find_if_not(lo, alignments.end(),
            [=](const GroupAlignment &g)->bool{ return i == g.templ->id; }
        );

        //ignore indels and count substitutions at each position in the
        //template; output is a matrix whose columns correspond to the
        //amino acid positions in the template and whose rows correspond
        //to the amino acids found at those positions 
        typedef decltype (alignments.cbegin()) SubsIt;
        auto count_substitutions = [&templ](SubsIt first, SubsIt last)->Matrix<float> {
            const size_t tpl_size = templ.size();
            Matrix<float> out(Aa::valid_chars.size(), tpl_size);
            for (; first != last; ++first) {
                const std::string &query = first->alignment;
                assert(tpl_size <= query.size());
                for (size_t q=0, t=0; t != tpl_size; ++q) {       //q and t are indices into query and template respectively
                    const char c = query[q];
                    if (c == '-')        { ++t; continue; }       //skip insertions
                    if (std::islower(c)) {      continue; }       //skip deletions
                    out.elem(Aa::from_char(c)->index(), t) += 1.; //increment the count in out[residues, position]
                    ++t;
                }
            }
            return out;
        };

        //accumulate mutation counts
        Matrix<float> substitutions = parallel_reduce(
            lo,
            hi,
            count_substitutions
        );

        //calculate the column totals
        std::vector<float> column_totals(substitutions.cols(), 0.);
        for (size_t r=0; r<substitutions.rows(); ++r)
        for (size_t c=0; c<substitutions.cols(); ++c) {
            column_totals[c] += substitutions.elem(r, c);
        }

        //convert counts to frequencies
        for (size_t c=0; c<substitutions.cols(); ++c) {
            if (column_totals[c] == 0.) continue; //treat 0/0 as 0
            for (size_t r=0; r<substitutions.rows(); ++r) substitutions.elem(r, c) /= column_totals[c];
        }

        //zero out the wild type frequencies
        for (size_t c=0; c<substitutions.cols(); ++c) {
            substitutions.elem(templ[c].index(), c) = 0.;
        }

        substitution_matrices.push_back(std::move(substitutions));
    }

    auto clock_stop = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(clock_stop-clock_start).count();
    size_t ss = static_cast<size_t>(ms / 1000);
    size_t mm = ss / 60;
    size_t hh = ss / 60;
    ms -= ss * 1000;
    ss -= mm * 60;
    mm -= hh * 60;

    std::time_t end_t = std::time(nullptr);
    std::tm end_tm = *std::localtime(&end_t);
    
    if (!p.no_header_flag) {
        std::cout << "#Settings#" << std::endl;
        std::cout << "#program version\t" << VERSION << std::endl;
        std::cout << "#run complete\t" << std::put_time(&end_tm, "%Y-%m-%d %H:%M:%S") << std::endl;
        std::cout << "#wall clock time\t" << std::setw(2) << std::setfill('0') << hh << ":"
                                          << std::setw(2) << std::setfill('0') << mm << ":"
                                          << std::setw(2) << std::setfill('0') << ss << "."
                                          << std::setw(3) << std::setfill('0') << static_cast<int>(ms) << std::endl; 
        std::cout << "#forward reads fastq file\t" << p.fw_filename << std::endl;
        std::cout << "#reverse reads fastq file\t" << p.rv_filename << std::endl;
        for (const UMIExtractor &fwex : fwexs) std::cout << "#forward nucleotide reference sequence (-f, --fw_ref)\t" << fwex.sequence() << std::endl;
        for (const UMIExtractor &rvex : rvexs) std::cout << "#reverse nucleotide reference sequence (-r, --rv_ref)\t" << rvex.sequence() << std::endl;
        if (!p.split_template_string.empty()) {
            std::cout << "#split template regular expression (--split)\t" << p.split_template_string << std::endl;
        }
        for (const help::TemplateSource &source : p.template_sources) {
            if  (std::holds_alternative<Aas>(source)) {
                std::cout << "#amino acid template sequence (-t, --template)\t" << std::get<Aas>(source) << std::endl;
            } else if (std::holds_alternative<Cdns>(source)) {
                std::cout << "#dna template sequence (-d, --template_dna)\t" << std::get<Cdns>(source).to_nts() << std::endl;
            } else if (std::holds_alternative<fs::path>(source)) {
                std::cout << "#template database (--template_db)\t" << std::get<fs::path>(source) << std::endl;
            }
        }
        std::cout << "#minimum 3 prime quality (-q, --min_qual)\t" << p.tp_qual_min << std::endl;
        std::cout << "#minimum umi group size (-g, --min_umi_grp)\t" << p.min_umi_group_size << std::endl;
        std::cout << "#reads aligned to template separately (-x, --skip_assembly)\t" << p.skip_assembly_flag << std::endl;
        std::cout << "#minimum nucleotide alignment overlap (-v, --min_overlap)\t" << p.min_overlap << std::endl;
        std::cout << "#maximum nucleotide mismatches allowed (-m, --max_mismatch)\t" << p.max_mismatches << std::endl;
        std::cout << "#minimum template alignment score (-a, --min_aln)\t" << p.min_alignment_score << std::endl;
        std::cout << "#Parse#" << std::endl; 
        std::cout << "#paired end reads parsed\t" << total_reads << std::endl;
        std::cout << "#reads filtered because of non-ATGC characters\t" << log.filter_invalid_chars << std::endl;
        std::cout << "#reads filtered because reference could not be identified in forward sequence\t" << log.filter_no_fw_umi << std::endl;
        std::cout << "#reads filtered because reference could not be identified in reverse sequence\t" << log.filter_no_rv_umi << std::endl;
        std::cout << "#reads filtered because they could not be assembled\t" << log.filter_could_not_assemble << std::endl;
        std::cout << "#reads filtered because of small umi group size\t" << log.filter_umi_group_size_too_small << std::endl;
        std::cout << "#reads merged during umi collapse\t" << log.filter_duplicate_umi << std::endl;
        std::cout << "#reads filtered because of premature stop codons\t" << log.filter_premature_stop_codon << std::endl;
        std::cout << "#reads filtered because no matching template was identified\t" << log.filter_no_matching_template << std::endl;
        std::cout << "#reads filtered because of poor alignment to template\t" << log.filter_bad_alignment << std::endl;
        std::cout << "#alignments calculated after qc and umi collapse\t" << alignments.size() << std::endl;
    }

    if (template_dbs.size()) {
        std::cout << "#Templates#" << std::endl;
        std::cout << "Template Id\tTemplate Name\tSequence" << std::endl;
        for (const auto& tpl : templates) {
            std::cout << tpl->id << '\t'
                << tpl->label() << '\t'
                << tpl->aas << std::endl;
        }

        //get frequency of template usage
        std::vector<Counter<std::string>> template_counters(template_dbs.size());
        for (const GroupAlignment& aln : alignments) {
            const AlignmentTemplate& tpl = *aln.templ;
            for (size_t i = 0; i < tpl.labels.size(); ++i) template_counters[i].push_back(tpl.labels[i]);
        }

        std::cout << "#Template Usage#" << std::endl;
        std::cout << "Split\tTemplate\tCount\tFrequency" << std::endl;
        for (size_t i = 0; i < template_counters.size(); ++i) {
            for (const auto& [label, count] : template_counters[i]) {
                std::cout << (i + 1) << '\t'
                    << label << '\t'
                    << count << '\t'
                    << count / static_cast<double>(template_counters[i].total()) << std::endl;
            }
        }
    }

    std::cout << "#Alignments#" << std::endl;
    std::cout << "Template\tUMI Group Size\tBarcode\tSequence" << std::endl;
    Cdns cdns; Nts nts; std::vector<std::optional<Cdn>> ocdns;
    for (const GroupAlignment &al : alignments) {
        std::cout << (al.templ ? std::to_string(al.templ->id) : std::string()) << '\t'
                  << al.umi_group_size << '\t'
                  << al.barcode << '\t'
                  << al.alignment << std::endl;
        switch (p.codon_output) {
        case help::CodonOutput::Ascii:
            std::cout << "\t\t\t" << al.cdns << std::endl;
            break;
        case help::CodonOutput::Horizontal:
            ocdns.clear();
            for (char c : al.cdns) ocdns.push_back(Cdn::from_char(c));
            std::cout << "\t\t\t"; //<< nts << std::endl;
            for (const auto &oc : ocdns) {
                if (oc) std::cout << oc->p1() << oc->p2() << oc->p3();
            }
            std::cout << std::endl;
            break;
        case help::CodonOutput::Vertical:
            ocdns.clear();
            for (char c : al.cdns) ocdns.push_back(Cdn::from_char(c));
            for (size_t i=0; i<3; ++i) {
                std::cout << "\t\t\t";
                for (size_t j=0; j<ocdns.size(); ++j) {
                    std::cout << (ocdns[j] ? static_cast<char>(ocdns[j]->at(i)) : ' ');
                }
                std::cout << std::endl;
            }
            break;
        case help::CodonOutput::None:
            break;
        };
    }

    if (template_dbs.size()) {
        lo = alignments.begin(); hi = alignments.begin();
        for (size_t i = 0; i < substitution_matrices.size(); ++i) {
            Matrix<float>& substitutions = substitution_matrices[i];

            lo = hi;
            hi = std::find_if_not(lo, alignments.end(),
                [&](const GroupAlignment& g)->bool { return templates[i]->id == g.templ->id; });

            std::cout << "#Substitutions (" << templates[i]->label() << ")#" << std::endl;
            //print the matrix
            for (size_t c = 0; c < substitutions.cols(); ++c) std::cout << '\t' << templates[i]->aas[c] << (c + p.number_from);
            std::cout << std::endl;
            for (size_t r = 0; r < substitutions.rows(); ++r) {
                std::cout << Aa::valid_chars[r];
                for (size_t c = 0; c < substitutions.cols(); ++c) std::cout << '\t' << substitutions.elem(r, c);
                std::cout << std::endl;
            }

            if (!templates[i]->cdns.empty()) {
                const Aas& aa_template = templates[i]->aas;
                const Cdns& cdn_template = templates[i]->cdns;

                //compare an alignment string/codons with an amino acid template/codons
                //and count the coding vs noncoding mutations
                typedef decltype (alignments.cbegin()) SubsIt;
                auto categorize_mutations = [&p, &aa_template, &cdn_template](SubsIt first, SubsIt last)->MutationCount {
                    MutationCount out(cdn_template.size());
                    assert(aa_template.size() == cdn_template.size());
                    const char* ta = aa_template.c_str();
                    const char* tc = cdn_template.c_str();
                    const size_t t_size = aa_template.size();

                    for (; first != last; ++first) {
                        assert(first->alignment.size() == first->cdns.size());
                        assert(aa_template.size() <= first->alignment.size());
                        const char* qa = first->alignment.c_str();
                        const char* qc = first->cdns.c_str();

                        for (size_t q = 0, t = 0; t != t_size; ++q) {       //q and t are indices into query and template
                            if (qa[q] == '-') { ++t; continue; } //skip deletions
                            if (std::islower(qa[q])) { continue; } //skip insertions
                            out.total[t] += 1;
                            if (qc[q] != tc[t]) {                       //codon mismatch means a mutation
                                if (qa[q] == ta[t]) {                   //mutation is synonymous if residues match
                                    out.synonymous[t] += 1;
                                }
                                else {
                                    out.nonsynonymous[t] += 1;
                                }
                            }
                            ++t;
                        }
                    }
                    return out;
                };

                MutationCount mutation_count = parallel_reduce(
                    lo,
                    hi,
                    categorize_mutations
                );

                std::cout << "#Mutation Counts (" << templates[i]->label() << ")#" << std::endl;
                for (size_t c = 0; c < aa_template.size(); ++c) std::cout << '\t' << aa_template[c] << (c + p.number_from);
                std::cout << std::endl;

                std::cout << "Total";
                for (size_t c = 0; c < aa_template.size(); ++c) std::cout << '\t' << mutation_count.total[c];
                std::cout << std::endl;

                std::cout << "Non-Coding";
                for (size_t c = 0; c < aa_template.size(); ++c) std::cout << '\t' << mutation_count.synonymous[c];
                std::cout << std::endl;

                std::cout << "Coding";
                for (size_t c = 0; c < aa_template.size(); ++c) std::cout << '\t' << mutation_count.nonsynonymous[c];
                std::cout << std::endl;
            }
        }
    }

    //output lists of unique amino acid and codon sequences
    //FIXME: de-duplicate this stuff
    if (!p.skip_assembly_flag) {
        struct Counts {
            size_t template_id = 0;
            const char* seq = "";
            unsigned groups = 0;
            unsigned reads = 0;
        };

        std::unordered_map<std::string, Counts> uniq;
        std::unordered_map<std::string, Counts> uniq_cdns;

        for (size_t i = 0; i < alignments.size(); ++i) {
            GroupAlignment& aln = alignments[i];
            std::erase(aln.alignment, '-'); //remove gaps from alignment string before calculating uniqueness
            std::erase(aln.cdns, ' '); //remove gaps from codons
            {
                auto [ii, success] = uniq.insert(std::make_pair(aln.alignment, Counts{}));
                ii->second.seq = ii->first.c_str();
                ii->second.groups += 1;
                ii->second.reads += static_cast<unsigned int>(aln.umi_group_size);
            }
            {
                auto [ii, success] = uniq_cdns.insert(std::make_pair(aln.cdns, Counts{}));
                ii->second.seq = ii->first.c_str();
                ii->second.groups += 1;
                ii->second.reads += static_cast<unsigned int>(aln.umi_group_size);
            }
        }

        //print Unique ORFs
        {
            std::cout << "#Unique Amino Acids (" << /*templates[i]->label()*/ "" << ")#" << std::endl;
            std::cout << "Num UMI Groups\tNum PCR Reads\tSequence" << std::endl;
            std::vector<Counts> flat;
            flat.reserve(uniq.size());
            for (const auto& [aln, count] : uniq) flat.push_back(count);
            std::sort(flat.begin(),
                flat.end(),
                [](const Counts& a, const Counts& b)->bool {return a.groups > b.groups; }
            );

            for (const Counts& c : flat) {
                std::cout << c.groups << '\t'
                    << c.reads << '\t'
                    << c.seq << std::endl;
            }
        }

        //print Unique Codons
        {
            std::cout << "#Unique Codons (" << /*templates[i]->label()*/ ""  << ")#" << std::endl;
            std::cout << "Num UMI Groups\tNum PCR Reads\tSequence" << std::endl;
            std::vector<Counts> flat;
            flat.reserve(uniq_cdns.size());
            for (const auto& [aln, count] : uniq_cdns) flat.push_back(count);
            std::sort(flat.begin(),
                flat.end(),
                [](const Counts& a, const Counts& b)->bool {return a.groups > b.groups; }
            );

            for (const Counts& c : flat) {
                std::cout << c.groups << '\t'
                    << c.reads << '\t'
                    << c.seq << std::endl;
            }
        }
    }

    return EXIT_SUCCESS;
}