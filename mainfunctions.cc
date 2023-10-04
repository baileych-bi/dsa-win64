#include "mainfunctions.h"

#include <algorithm>
#include <numeric>
#include <thread>
#include <unordered_map>

#include "abs.h"
#include "align.h"
#include "io.h"
#include "parallelism.h"
#include "umi.h"

#undef min
#undef max

namespace bio {

using help::Params;

void
vector_accumulate(std::vector<unsigned> &a, const std::vector<unsigned> &b) {
    assert(a.size() == b.size());
    auto ii=a.begin(); auto jj=b.cbegin();
    for (; jj != b.cend(); ++ii, ++jj) *ii += *jj;
}

MutationCount::MutationCount(size_t cols)
    : synonymous(cols, 0)
    , nonsynonymous(cols, 0)
    , total(cols, 0) { }

MutationCount
MutationCount::operator+(const MutationCount &c) const {
    assert(synonymous.size()    == c.synonymous.size()   );
    assert(nonsynonymous.size() == c.nonsynonymous.size());
    assert(total.size()         == c.total.size()       );

    MutationCount sum(synonymous.size());
    vector_accumulate(sum.synonymous,      synonymous);
    vector_accumulate(sum.synonymous,    c.synonymous);
    vector_accumulate(sum.nonsynonymous,   nonsynonymous);
    vector_accumulate(sum.nonsynonymous, c.nonsynonymous);
    vector_accumulate(sum.total,           total);
    vector_accumulate(sum.total,         c.total);

    return sum;
}

std::string
AlignmentTemplate::label(const std::string &delim) const {
    std::string s;
    for (size_t i=1; i<labels.size(); ++i) s += labels[i-1] + delim;
    if (!labels.empty()) s += labels.back();
    return s;
}


std::vector<Read>
extract_read_data(const ConstMapping &mapping) {
    const unsigned int thread_count = std::thread::hardware_concurrency();

    //divide the memory up into evenly sized chunks
    size_t chunk = mapping.size() / thread_count;
    std::vector<const char *> breakpoints(thread_count+1, nullptr);
    for (size_t i=0; i<thread_count; ++i) breakpoints[i] = mapping.begin() + i * chunk;
    breakpoints.back() = mapping.end();

    //move each chunk pointer to the beginning of the next record
    for (size_t i=1; i<breakpoints.size()-1; ++i) {
        breakpoints[i] = seek_next(breakpoints[i], mapping.begin(), breakpoints.back());
    }

    auto process_fastq = [](const char *begin,
                            const char *end,
                            std::vector<Read> &result)->void {
        result.clear();
        Read rd;
        size_t stripped = 0;
        while (begin != end) {
            begin = bio::skipline(begin, end, '\n');            //skip header
            begin = bio::getline(begin, end, rd.dna, stripped); //copy dna
            begin = bio::skipline(begin, end, '\n');            //skip '+'
            begin = bio::getline(begin, end, rd.qual);          //copy quality

            if (stripped != 0 || rd.dna.size() != rd.qual.size()) {
                result.push_back(Read());
            } else {
                result.push_back(rd);
            }
        }
    };

    //spawn threads to perform the extraction
    std::vector<std::thread> threads(thread_count-1);
    vecvec<Read>             partial_results(thread_count);

    size_t i=0;
    for (; i<threads.size(); ++i) {
        const char *begin = breakpoints[i  ];
        const char *end   = breakpoints[i+1];
        threads[i] = std::thread(process_fastq,
                                 begin,
                                 end,
                                 std::ref(partial_results[i]));
    }
    process_fastq(breakpoints[i], breakpoints[i+1], partial_results[i]);

    for (auto &th : threads) th.join();

    //concat partial results from each thread in a linked list
    std::vector<Read> result;
    size_t total_size = 0;
    for (const auto &pr : partial_results) total_size += pr.size();
    result.reserve(total_size);

    for (auto &pr : partial_results) {
        result.insert(
            result.cend(),
            std::make_move_iterator(pr.begin()),
            std::make_move_iterator(pr.end()  ));
    }

    return result;
}

std::vector<ReadPair>
qc_reads(
    std::vector<Read> &&fw,
    std::vector<Read> &&rv,
    const std::vector<UMIExtractor> &fwexs,
    const std::vector<UMIExtractor> &rvexs,
    const Params &params,
    ParseLog &log)
{
    assert(fw.size() == rv.size());
    std::vector<ReadPair> result;

    const unsigned int thread_count = std::thread::hardware_concurrency();
    const size_t chunk = fw.size() / thread_count;

    std::vector<std::thread> threads(thread_count-1);
    vecvec<ReadPair>         partial_results(thread_count);
    std::vector<ParseLog>    partial_logs(thread_count);

    typedef std::vector<Read>::iterator IterT;

    auto perform_qc = [&](IterT ff, 
                        IterT rr,
                        size_t n,
                        std::vector<ReadPair> &output,
                        ParseLog &log)->void {
        output.clear();
        for (const IterT last = ff + n; ff != last; ++ff, ++rr) {
            if (ff->empty()) {
                ++log.filter_invalid_chars;
                continue;
            }

            if (rr->empty()) {
                ++log.filter_invalid_chars;
                continue;
            }

            while (!ff->empty() && ff->qual.back() < params.tp_qual_min) ff->pop_back();
            while (!rr->empty() && rr->qual.back() < params.tp_qual_min) rr->pop_back();

            ExtractedUMI fwumi;
            for (const UMIExtractor &fwex : fwexs) {
                fwumi = fwex(ff->dna.begin(), ff->dna.end());
                if (fwumi.valid()) break;
            }
            if (fwumi.invalid()) {
                ++log.filter_no_fw_umi;
                continue;
            }

            ExtractedUMI rvumi;
            for (const UMIExtractor &rvex : rvexs) {
                rvumi = rvex(rr->dna.begin(), rr->dna.end());
                if (rvumi.valid()) break;
            }
            if (rvumi.invalid()) {
                ++log.filter_no_rv_umi;
                continue;
            }

            ff->dna.exo(fwumi.from + fwumi.length, 0);
            ff->qual = ff->qual.substr(fwumi.from + fwumi.length);

            rr->dna.exo(rvumi.from + rvumi.length, 0);
            rr->qual = rr->qual.substr(rvumi.from + rvumi.length);

            ff->barcode.reserve(fwumi.barcode.size() + rvumi.barcode.size());
            ff->barcode += fwumi.barcode;
            ff->barcode += rvumi.barcode;

            ReadPair rd;
            rd.fw = std::move(*ff);
            rd.rv = std::move(*rr);

            output.push_back(std::move(rd));
        }
    };

    size_t i=0;
    auto ff=fw.begin(), rr=rv.begin();
    for (; i<thread_count-1; ++i, ff += chunk, rr += chunk) {
        threads[i] = std::thread(
            perform_qc, ff, rr, chunk, std::ref(partial_results[i]), std::ref(partial_logs[i])
        );
    }
    perform_qc(ff, rr, fw.end()-ff, partial_results[i], partial_logs[i]);

    for (auto &th : threads) th.join();

    size_t result_size = 0;
    for (const auto &pr : partial_results) result_size += pr.size();
    result.reserve(result_size);

    for (auto &pr : partial_results) {
        result.insert(
            result.cend(),
            std::make_move_iterator(pr.begin()),
            std::make_move_iterator(pr.end()  )
        );
    }

    log = std::accumulate(partial_logs.begin(), partial_logs.end(), log);

    fw.clear(); fw.shrink_to_fit();
    rv.clear(); rv.shrink_to_fit();

    return result;
}

std::vector<Read>
assemble_reads(
    std::vector<ReadPair> &&pairs,
    const Params &params,
    ParseLog &log) {
    std::vector<Read> result;

    const unsigned int thread_count = 1; //std::thread::hardware_concurrency();
    const size_t chunk = pairs.size() / thread_count;

    std::vector<std::thread> threads(thread_count-1);
    vecvec<Read>             partial_results(thread_count);
    std::vector<ParseLog>    partial_logs(thread_count);

    typedef std::vector<ReadPair>::iterator IterT;

    auto perform_assembly = [&](IterT ii, 
                        size_t n,
                        std::vector<Read> &output,
                        ParseLog &log)->void {

        output.clear();
        for (const IterT last = ii + n; ii != last; ++ii) {
            Read rd = Read::assemble(
                std::move(ii->fw),
                std::move(ii->rv),
                params.min_overlap,
                params.max_mismatches);
            
            if (rd.empty()) {
                ++log.filter_could_not_assemble;
                continue;
            }

            output.push_back(std::move(rd));
        }
    };

    size_t i=0;
    auto ii=pairs.begin();
    for (; i<thread_count-1; ++i, ii += chunk) {
        threads[i] = std::thread(
            perform_assembly, ii, chunk, std::ref(partial_results[i]), std::ref(partial_logs[i])
        );
    }
    perform_assembly(ii, pairs.end()-ii, partial_results[i], partial_logs[i]);

    for (auto &th : threads) th.join();

    size_t result_size = 0;
    for (const auto &pr : partial_results) result_size += pr.size();
    result.reserve(result_size);

    for (auto &pr : partial_results) {
        result.insert(
            result.cend(),
            std::make_move_iterator(pr.begin()),
            std::make_move_iterator(pr.end()  )
        );
    }

    log = std::accumulate(partial_logs.begin(), partial_logs.end(), log);

    pairs.clear();
    pairs.shrink_to_fit();

    return result;
}

struct Choice {
    Nt nt = Nt::A;
    unsigned occurs   = 0;
    char     max_qual = 0;

    bool operator<(const Choice &c) const {
        if (occurs   < c.occurs  ) return true;
        if (occurs   > c.occurs  ) return false;
        if (max_qual < c.max_qual) return true;
        return false;
    }
};

std::array<Choice, 5>
make_default_choices() {
    std::array<Choice, 5> defaults;
    defaults[Nt::A.index()] = {Nt::A, 0, 0};
    defaults[Nt::C.index()] = {Nt::C, 0, 0};
    defaults[Nt::G.index()] = {Nt::G, 0, 0};
    defaults[Nt::T.index()] = {Nt::T, 0, 0};
    defaults[Nt::N.index()] = {Nt::N, 0, 0};

    return defaults;
}

void
build_consensus_sequence(std::vector<Read> &reads, const Params &params, bool ragged_ends) {
    assert (reads.size() >= params.min_umi_group_size);

    static const std::array<Choice, 5> default_choices = make_default_choices();

    thread_local std::vector<std::array<Choice, 5>> choices;
    choices.clear();

    //for ragged ends we left-justify all sequences and make sure that every nucleotide in
    //the consensus output draws from least params.min_umi_group_size sequences
    if (ragged_ends) {
        struct sort_by_size_desc {
            bool operator()(const Read &a, const Read &b) const { return a.size() > b.size(); }
        };
        std::sort(reads.begin(), reads.end(), sort_by_size_desc());
        choices.resize(reads[params.min_umi_group_size-1].size(), default_choices);
        reads.front().umi_group_size = reads.size();

        for (const Read &rd : reads) {
            const size_t limit = std::min(choices.size(), rd.size());
            for (size_t i=0; i<limit; ++i) {
                Choice &ch = choices[i][rd.dna[i].index()];
                ch.occurs += 1;
                if (rd.qual[i] > ch.max_qual) ch.max_qual = rd.qual[i];
            }
        }
    //!ragged_ends means sequences are expected to be the same length, any indels
    //should be considered PCR artifacts; therefore we determine the modal length
    //of the UMI group and discard sequences that are too long or too short
    } else {
        std::unordered_map<size_t, size_t> size_counts;
        for (const Read &rd : reads) size_counts[rd.size()] += 1;

        size_t modal_size = std::max_element(
            size_counts.begin(),
            size_counts.end(),
            [](const std::pair<size_t, size_t> &a, 
            const std::pair<size_t, size_t> &b)->bool{
                return a.second < b.second;
            })->first;
        
        choices.resize(modal_size, default_choices);

        reads.front().umi_group_size = 0;
        for (const Read &rd : reads) {
            if (rd.size() != modal_size) continue;
            reads.front().umi_group_size += 1;

            const size_t limit = std::min(rd.size(), choices.size());
            for (size_t i=0; i<limit; ++i) {
                Choice &ch = choices[i][rd.dna[i].index()];
                ch.occurs += 1;
                if (rd.qual[i] > ch.max_qual) ch.max_qual = rd.qual[i];
            }
        }
    }

    reads.front().resize(choices.size());
    for (size_t i=0; i<choices.size(); ++i) {
        Choice &ch = *std::max_element(choices[i].begin(), choices[i].end());
        reads.front().dna[i]  = ch.nt;
        reads.front().qual[i] = ch.max_qual;
    }

    reads.resize(1);
}

std::vector<Read>
umi_collapse(
    std::vector<Read> &&reads,
    const Params &params,
    ParseLog &log,
    bool ragged_ends) {
    std::vector<Read> result;

    //gather up reads by umi into a hash table
    std::unordered_map<std::string, std::vector<Read>> groups;
    for (Read &rd : reads) groups[rd.barcode].push_back(std::move(rd));
    reads.clear();

    const unsigned int thread_count = std::thread::hardware_concurrency();
    const size_t chunk = groups.size() / thread_count;

    std::vector<std::thread> threads(thread_count-1);
    vecvec<Read>             partial_results(thread_count);
    std::vector<ParseLog>    partial_logs(thread_count);

    typedef decltype(groups.begin()) IterT;
 
    auto build_consensus = [&](IterT from, IterT to, std::vector<Read> &output, ParseLog &log)->void {
        output.clear();
        for (; from != to; ++from) {
            std::vector<Read> &group = from->second;
            size_t pre_consensus_group_size = group.size();

            if (group.size() < params.min_umi_group_size) {
                log.filter_umi_group_size_too_small += pre_consensus_group_size;
                continue;
            }

            if (group.size() > 1) {
                build_consensus_sequence(group, params, ragged_ends);
            }

            if (group.front().umi_group_size < params.min_umi_group_size) {
                log.filter_umi_group_size_too_small += pre_consensus_group_size;
                continue;
            }
            
            if (std::find(group.front().dna.begin(), group.front().dna.end(), Nt::N) != group.front().dna.end()) {
                ++log.filter_invalid_chars;
                continue;
            }

            log.filter_duplicate_umi += pre_consensus_group_size - group.size();
            output.push_back(std::move(group.front()));
        }
    };

    size_t i=0;
    auto from=groups.begin(), to=from;
    for (; i<thread_count-1; ++i) {
        std::advance(to, chunk);
        threads[i] = std::thread(
            build_consensus, 
            from, 
            to,
            std::ref(partial_results[i]),
            std::ref(partial_logs[i]));
        from = to;
    }
    build_consensus(from, groups.end(), partial_results[i], partial_logs[i]);

    for (auto &th : threads) th.join();

    size_t result_size = 0;
    for (const auto &pr : partial_results) result_size += pr.size();
    result.reserve(result_size);

    for (auto &pr : partial_results) {
        result.insert(
            result.cend(),
            std::make_move_iterator(pr.begin()),
            std::make_move_iterator(pr.end()  ));
    }

    log = std::accumulate(partial_logs.begin(), partial_logs.end(), log);

    reads.clear(); reads.shrink_to_fit();

    return result;
}

std::vector<Orf>
translate_and_filter_ptcs(std::vector<Read> &&preads, 
                          const help::Params &p,
                          ParseLog &log,
                          bool reverse_complement)
{
    //translate and discard all reads with PTCs
    auto translate_and_filter_ptcs = [&p, &reverse_complement](Read &&rd, ParseLog &log)->std::optional<Orf> {
        std::optional<Orf> opt_orf;

        if (reverse_complement) {
            rd.resize(rd.size() / 3 * 3);
            rd.reverse_complement();
        }

        Orf orf = std::move(rd);
        if (orf.contains_ptc()) {
            ++log.filter_premature_stop_codon;
        } else {
            opt_orf = std::move(orf);
        }
        return opt_orf;
    };

    std::vector<Orf> orfs;
    orfs.reserve(preads.size());

    parallel_transform_filter(
        std::make_move_iterator(preads.begin()),
        std::make_move_iterator(preads.end()  ),
        std::back_inserter(orfs),
        translate_and_filter_ptcs,
        log
    );

    preads.clear(); preads.shrink_to_fit();

    return orfs;
};

std::vector<GroupAlignment>
align_to_multiple_templates(vecvec<Orf> &&orfs,
                   const std::vector<std::shared_ptr<const TemplateDatabase>> &dbs,
                   const help::Params &params,
                   ParseLog &log,
                   bool ragged_ends) {
    assert (!dbs.empty());
    std::vector<GroupAlignment> alignments;
    if (orfs.empty()) {
        return alignments;
    }

    struct WorkerOutput {
        GroupAlignment alignment;
        std::vector<size_t> template_ids;
    };

    auto align_to_multiple_templates_worker = [&](
        std::vector<Orf> &&orfs,
        ParseLog &log)->std::optional<WorkerOutput> {
        assert (orfs.size() == dbs.size());

        std::optional<WorkerOutput> output;

        GroupAlignment           alignment;
        std::vector<size_t>   template_ids;
        std::string          full_template;

        for (size_t i=0; i<orfs.size(); ++i) {
            if (dbs[i] == nullptr) {
                template_ids.push_back(0);
                alignment.alignment += orfs[i].aas.c_str();
                alignment.cdns += orfs[i].cdns.c_str();
                continue;
            }

            Alignment aln;

            size_t template_id = dbs[i]->codon_data_available()
                               ? dbs[i]->query_and_align(orfs[i].cdns, aln)
                               : dbs[i]->query_and_align(orfs[i].aas,  aln);

            if (template_id == TemplateDatabase::NOT_FOUND) {
                ++log.filter_no_matching_template;
                break;
            }

            const Aas  &template_aas  = dbs[i]->get_aas(template_id);
            const Cdns &template_cdns = dbs[i]->get_codons(template_id);

            float max_score = static_cast<float>(
                dbs[i]->codon_data_available()
                    ? nw_self_align_score<Cdn>(template_cdns, CDNSUBS )
                    : nw_self_align_score<Aa> (template_aas , BLOSUM62)
            );

            if (ragged_ends) max_score -= dbs[i]->gap_penalty()
                                        * std::abs(static_cast<int32_t>( orfs[i].aas.size()) - 
                                                   static_cast<int32_t>(template_aas.size()));

            if (aln.score / max_score < params.min_alignment_score) {
                ++log.filter_bad_alignment;
                break;
            }

            template_ids.push_back(template_id);
            alignment.alignment += aln.build_string(orfs[i].aas);
            alignment.cdns      += aln.build_string(orfs[i].cdns);
        }

        if (template_ids.size() == orfs.size()) {
            alignment.umi_group_size = orfs.front().umi_group_size;
            alignment.barcode        = orfs.front().barcode;
            output = WorkerOutput{std::move(alignment), std::move(template_ids)};
        }

        return output;
    };

    std::vector<WorkerOutput> worker_outputs;
    worker_outputs.reserve(orfs.size());

    parallel_transform_filter(
        std::make_move_iterator(orfs.begin()),
        std::make_move_iterator(orfs.end()  ),
        std::back_inserter(worker_outputs),
        align_to_multiple_templates_worker,
        log
    );

    orfs.clear();

    //hash function from Thomas Mueller, Stack Overflow #664014
    struct vector_size_t_hasher {
        size_t operator()(const std::vector<size_t> &ids) const {
            size_t h = ids.size();
            for (size_t x : ids) {
                x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9u;
                x = (x ^ (x >> 27)) * 0x94d049bb133111ebu;
                x =  x ^ (x >> 31);
                h ^= x;
            }
            return h;
        }
    };

    std::unordered_map<std::vector<size_t>,
                       std::shared_ptr<AlignmentTemplate>,
                       vector_size_t_hasher,
                       std::equal_to<std::vector<size_t>>> template_lookup;

    size_t next_id = 0;
    for (WorkerOutput &wo : worker_outputs) {
        auto [ii, inserted] = template_lookup.insert(
            {std::move(wo.template_ids), std::shared_ptr<AlignmentTemplate>()}
        );
        if (inserted) {
            std::shared_ptr<AlignmentTemplate> tpl = std::shared_ptr<AlignmentTemplate>(new AlignmentTemplate);
            tpl->id = ++next_id;
            const std::vector<size_t> &ids = ii->first;
            const std::string nil;
            for (size_t i=0; i<ids.size(); ++i) {
                size_t id = ids[i];
                if (dbs[i]) {
                    tpl->labels.push_back(dbs[i]->get_label(id));
                    tpl->aas += dbs[i]->get_aas(id);
                    tpl->cdns += dbs[i]->get_codons(id);
                } else {
                    tpl->labels.push_back("none");
                }
            }
            ii->second = tpl;
        }
        wo.alignment.templ = ii->second;
    }

    alignments.reserve(worker_outputs.size());
    for (WorkerOutput &wo : worker_outputs) alignments.push_back(std::move(wo.alignment));

    return alignments;
}

vecvec<Orf>
split_orfs(std::vector<Orf> &&orfs,
           const help::Params &params,
           ParseLog &log) {
    vecvec<Orf> result;
    result.reserve(orfs.size());

    //if there's nothing to split we just return a matrix of shape (orfs.size(), 1)
    if (params.split_template_regex.mark_count() == 0) {
        result.resize(orfs.size());
        for (size_t i=0; i<orfs.size(); ++i) result[i].push_back(std::move(orfs[i]));
        return result;
    }

    //otherwise the matrix shape is (orfs.size() - unsplittable, params.split_template_regex.mark_count())
    std::cmatch match;
    for (Orf &orf : orfs) {
        std::vector<Orf> splits;
        splits.resize(params.split_template_regex.mark_count());
        if (std::regex_match(orf.aas.c_str(), match, params.split_template_regex)) {
            for (size_t i=1; i<match.size(); ++i) {
                Orf &sub = splits[i-1];
                sub.umi_group_size = orf.umi_group_size;
                sub.template_id    = orf.template_id;
                sub.barcode        = orf.barcode;
                sub.aas            = orf.aas .subclone(match.position(i), match.length(i));
                sub.cdns           = orf.cdns.subclone(match.position(i), match.length(i));
            }
            result.push_back(std::move(splits));
        } else {
            ++log.filter_split_failed;
        }
    }

    return result;
}

}; //namespace bio