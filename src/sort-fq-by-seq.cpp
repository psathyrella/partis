#include "kseq.h"

#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <functional>
#include <vector>
#include <string>

#include <zlib.h>

using std::string;
using std::vector;

/* KLIB */
KSEQ_INIT(gzFile, gzread);

struct GzFile
{
    GzFile(const char *path, const char *mode) :
        fp(gzopen(path, mode))
    {
        assert(fp != nullptr && "Failed to open gz file.");
    }

    ~GzFile() { gzclose(fp); }

    gzFile fp;
};

string blank_if_null(char *s) {
    return string(s == nullptr ? "" : s);
}

struct Sequence
{
    string name;
    string description;
    string seq;
    string qual;
};

int main(int argc, char *argv[])
{
    if(argc != 2) {
        std::cerr << "usage: " << argv[0] << " " << "<fast[aq]>\n";
        return EXIT_FAILURE;
    }
    assert(argc == 2);
    GzFile in_fp(argv[1], "r");
    kseq_t *seq = kseq_init(in_fp.fp);

    vector<Sequence> sequences;
    while(kseq_read(seq) > 0) {
        Sequence s { blank_if_null(seq->name.s),
                     blank_if_null(seq->comment.s),
                     blank_if_null(seq->seq.s),
                     blank_if_null(seq->qual.s) };
        sequences.push_back(std::move(s));
    }
    kseq_destroy(seq);

    auto cmp = [](const Sequence& s1, const Sequence& s2) { return s1.seq < s2.seq; };
    std::sort(sequences.begin(), sequences.end(), cmp);

    for(const Sequence& s : sequences) {
        if(!s.qual.empty()) { /* write fastq */
            std::cout << '@' << s.name;
            if(!s.description.empty())
                std::cout << ' ' << s.description;
            std::cout << '\n' << s.seq << "\n+\n" << s.qual << '\n';
        } else {
            std::cout << '>' << s.name;
            if(!s.description.empty())
                std::cout << ' ' << s.description;
            std::cout << '\n' << s.seq << '\n';
        }
    }

    return EXIT_SUCCESS;
}
