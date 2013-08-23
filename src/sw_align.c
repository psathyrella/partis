/* Maybe:
 * - read a batch of sequences
 * - Create some threads
 * - process / write
 * - read more
 */

/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 07/31/12.
 */

#include "sw_align.h"

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <zlib.h>
#include <pthread.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "kseq.h"
#include "ksort.h"
#include "ksw.h"
#include "kvec.h"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place)
  @discussion x will be modified.
  */
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

KSEQ_INIT(gzFile, gzread);

typedef kvec_t(kseq_t) kseq_v;

/* Destroy a vector of pointers, calling f on each item */
#define kvp_destroy(f, v) \
    for(size_t __kpd = 0; __kpd < kv_size(v); ++__kpd) \
        f(kv_A(v, __kpd)); \
    kv_destroy(v)

#define kvi_destroy(f, v) \
    for(size_t __kpd = 0; __kpd < kv_size(v); ++__kpd) \
        f(&kv_A(v, __kpd)); \
    kv_destroy(v)

static void kstring_copy(kstring_t* dest, const kstring_t* other) {
    dest->l = other->l;
    dest->m = dest->l + 1;
    if(dest->l > 0) {
        assert(other->s != NULL);
        dest->s = malloc(dest->l + 1);
        memcpy(dest->s, other->s, dest->l);
        dest->s[dest->l] = 0; // null-terminate
    } else {
        dest->s = NULL;
    }
}

static void kseq_stack_destroy(kseq_t* seq)
{
    free(seq->name.s);
    free(seq->comment.s);
    free(seq->seq.s);
    free(seq->qual.s);
}

static void kseq_copy(kseq_t* dest, const kseq_t* seq)
{
    dest->f = NULL;
    kstring_copy(&dest->name, &seq->name);
    kstring_copy(&dest->comment, &seq->comment);
    kstring_copy(&dest->seq, &seq->seq);
    kstring_copy(&dest->qual, &seq->qual);
}

kseq_v read_seqs(kseq_t* seq,
                 size_t n_wanted) {
    kseq_v result;
    kv_init(result);
    for(size_t i = 0; i < n_wanted || n_wanted == 0; i++) {
        if(kseq_read(seq) <= 0)
            break;
        kseq_t s;
        kseq_copy(&s, seq);
        kv_push(kseq_t, result, s);
    }
    return result;
}

typedef struct {
    size_t target_idx;
    kswr_t loc;
    uint32_t *cigar;
    int n_cigar;
} aln_t;
typedef kvec_t(aln_t) aln_v;

#define __aln_score_lt(a, b) ((a).loc.score > (b).loc.score)
KSORT_INIT(dec_score, aln_t, __aln_score_lt)

typedef struct {
    aln_v alignments;
} aln_task_t;
typedef kvec_t(aln_task_t) aln_task_v;

typedef struct {
    int32_t gap_o; /* 3 */
    int32_t gap_e; /* 1 */
    uint8_t *table;
    int m; /* Number of residue tyes */
    int8_t *mat; /* Scoring matrix */
} align_config_t;

static aln_v align_read(const kseq_t* read,
                        const kseq_v targets,
                        const align_config_t *conf)
{
    kseq_t *r;
    const int32_t read_len = read->seq.l;

    aln_v result;
    kv_init(result);
    kv_resize(aln_t, result, kv_size(targets));

    uint8_t *read_num = calloc(read_len, sizeof(uint8_t));

    for (size_t k = 0; k < read_len; ++k)
        read_num[k] = conf->table[(int)read->seq.s[k]];

    // Align to each target
    kswq_t *qry = NULL;
    for(size_t j = 0; j < kv_size(targets); j++) {
        // Encode target
        r = &kv_A(targets, j);
        uint8_t *ref_num = calloc(r->seq.l, sizeof(uint8_t));
        for (size_t k = 0; k < r->seq.l; ++k)
            ref_num[k] = conf->table[(int)r->seq.s[k]];

        aln_t aln;
        aln.target_idx = j;
        aln.loc = ksw_align(read_len, read_num,
                            r->seq.l, ref_num,
                            conf->m,
                            conf->mat,
                            conf->gap_o,
                            conf->gap_e,
                            KSW_XSTART,
                            &qry);
        ksw_global(aln.loc.qe - aln.loc.qb + 1,
                   &read_num[aln.loc.qb],
                   aln.loc.te - aln.loc.tb + 1,
                   &ref_num[aln.loc.tb],
                   conf->m,
                   conf->mat,
                   conf->gap_o,
                   conf->gap_e,
                   40, /* TODO: Magic number */
                   &aln.n_cigar,
                   &aln.cigar);
        kv_push(aln_t, result, aln);
        free(ref_num);
    }
    free(qry);
    free(read_num);
    ks_introsort(dec_score, kv_size(result), result.a);
    return result;
}

static void write_sam_records(gzFile dest,
                              const kseq_t *read,
                              const aln_v result,
                              const kseq_v ref_seqs)
{
    for(size_t i = 0; i < kv_size(result); i++) {
        aln_t a = kv_A(result, i);
        gzprintf(dest, "%s\t%d\t", read->name.s,
                i == 0 ? 0 : 256); // Secondary
        gzprintf(dest, "%s\t%d\t%d\t",
                kv_A(ref_seqs, a.target_idx).name.s, /* Reference */
                a.loc.tb + 1,                        /* POS */
                40);                                 /* MAPQ */
        if (a.loc.qb)
            gzprintf(dest, "%dS", a.loc.qb);
        for (size_t c = 0; c < a.n_cigar; c++) {
            int32_t letter = 0xf&*(a.cigar + c);
            int32_t length = (0xfffffff0&*(a.cigar + c))>>4;
            gzprintf(dest, "%d", length);
            if (letter == 0) gzprintf(dest, "M");
            else if (letter == 1) gzprintf(dest, "I");
            else gzprintf(dest, "D");
        }

        if (a.loc.qe + 1 != read->seq.l)
            gzprintf(dest, "%luS", read->seq.l - a.loc.qe - 1);

        gzprintf(dest, "\t*\t0\t0\t");
        gzprintf(dest, "%s\t", i > 0 ? "*" : read->seq.s);
        if (read->qual.s && i == 0)
            gzprintf(dest, "%s", read->qual.s);
        else
            gzprintf(dest, "*");

        gzprintf(dest, "\tAS:i:%d\n", a.loc.score);
    }
}

void align_reads (const char* ref_path,
                  const char* qry_path,
                  const char* output_path,
                  const int32_t match,       /* 2 */
                  const int32_t mismatch,    /* 2 */
                  const int32_t gap_o,       /* 3 */
                  const int32_t gap_e,       /* 1 */
                  const uint8_t n_threads) { /* 1 */
    gzFile read_fp, ref_fp, out_fp;
    int32_t j, k, l;
    const int m = 5;
    kseq_t *seq;
    int8_t* mat = (int8_t*)calloc(25, sizeof(int8_t));

    /* This table is used to transform nucleotide letters into numbers. */
    uint8_t table[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };

    // initialize scoring matrix for genome sequences
    for (l = k = 0; LIKELY(l < 4); ++l) {
        for (j = 0; LIKELY(j < 4); ++j) mat[k++] = l == j ? match : -mismatch;	/* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base
    }
    for (j = 0; LIKELY(j < 5); ++j) mat[k++] = 0;


    // Read reference sequences
    ref_fp = gzopen(ref_path, "r");
    seq = kseq_init(ref_fp);
    kseq_v ref_seqs;
    ref_seqs = read_seqs(seq, 0);
    kseq_destroy(seq);
    gzclose(ref_fp);

    fprintf(stderr, "Read %lu references\n",
            kv_size(ref_seqs));

    // Print SAM header
    out_fp = gzopen(output_path, "w");
    gzprintf(out_fp, "@HD\tVN:1.4\tSO:queryname\n");
    for(size_t i = 0; i < kv_size(ref_seqs); i++) {
        seq = &kv_A(ref_seqs, i);
        gzprintf(out_fp, "@SQ\tSN:%s\tLN:%d\n",
                seq->name.s, (int32_t)seq->seq.l);
    }

    align_config_t conf;
    conf.gap_o = gap_o;
    conf.gap_e = gap_e;
    conf.m = m;
    conf.table = table;
    conf.mat = mat;

    read_fp = gzopen(qry_path, "r");
    size_t count = 0;
    seq = kseq_init(read_fp);
    while(true) {
        kseq_v reads = read_seqs(seq, 5); // TODO: Magic number
        count += kv_size(reads);
        if(!kv_size(reads)) {
            break;
        }

        const size_t n_reads = kv_size(reads);
        for(size_t i = 0; i < n_reads; i++) {
            kseq_t *s = &kv_A(reads, i);
            aln_v result = align_read(s,
                                      ref_seqs,
                                      &conf);

            write_sam_records(out_fp,
                              s,
                              result,
                              ref_seqs);

            for(size_t j = 0; j < kv_size(result); j++)
                free(kv_A(result, j).cigar);
            kv_destroy(result);
            kseq_stack_destroy(s);
        }
        kv_destroy(reads);
    }
    kseq_destroy(seq);

    // Clean up reference sequences
    kvi_destroy(kseq_stack_destroy, ref_seqs);

    gzclose(read_fp);
    gzclose(out_fp);
    free(mat);
}
