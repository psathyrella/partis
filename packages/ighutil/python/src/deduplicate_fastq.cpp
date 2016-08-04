#include "deduplicate_fastq.hpp"

#include "kseq.h"
#include "zlib.h"

#include <cassert>

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using std::string;
using std::unordered_map;
using std::vector;

KSEQ_INIT(gzFile, gzread);

struct ReadWithCount
{
    ReadWithCount(const string& name, const string& seq, const size_t count=0) :
        name(name), seq(seq), count(count) {};
    string name, seq;
    size_t count;
};

struct GzFile
{
    GzFile(const char *path, const char *mode) :
        fp(gzopen(path, mode))
    {
        assert(fp != nullptr);
    }

    ~GzFile() { gzclose(fp); }

    gzFile fp;
};

void deduplicate_fastq(const char* in_path,
                       const char* out_path,
                       const int est_size)
{
    GzFile in_fp(in_path, "r");
    GzFile out_fp(out_path, "w");

    kseq_t *seq = kseq_init(in_fp.fp);

    unordered_map<string, ReadWithCount> result;
    result.reserve(est_size);

    size_t n_processed = 0, n_unique = 0;

    while(kseq_read(seq) > 0) {
        n_processed++;

        string s(seq->seq.s);

        unordered_map<string, ReadWithCount>::iterator it = result.find(s);
        if(it == result.end()) {
            const string name(seq->name.s);
            it = result.insert(std::make_pair(s, ReadWithCount(name, s))).first;
            n_unique++;
        }
        it->second.count++;

        if(n_processed % 50000 == 0) {
            std::cout << "[dedup_fq] "
                << std::right << std::setw(10) << n_unique << "/"
                << std::right << std::setw(10) << n_processed << " unique\r";
        }
    }
    kseq_destroy(seq);

    vector<ReadWithCount> reads;
    reads.reserve(n_unique);

    for(auto p : result)
        reads.push_back(std::move(p.second));
    result.clear();

    /* Sort by descending frequency */
    auto dec_count = [](const ReadWithCount& a, const ReadWithCount& b) {
        return a.count > b.count;
    };

    std::sort(reads.begin(), reads.end(), dec_count);

    for(const ReadWithCount& r : reads) {
        gzprintf(out_fp.fp, ">%s_%lu\n%s\n", r.name.c_str(), r.count, r.seq.c_str());
    }
}
