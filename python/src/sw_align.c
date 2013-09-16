#include "sw_align.h"

#include <assert.h>
#include <emmintrin.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <zlib.h>
#include "kseq.h"
#include "ksort.h"
#include "kstring.h"
#include "ksw.h"
#include "kvec.h"

/* Cigar operations */
/**
 * Describing how CIGAR operation/length is packed in a 32-bit integer.
 */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

#define BAM_CIGAR_STR  "MIDNSHP=XB"
#define BAM_CIGAR_TYPE 0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c)>>BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l)<<BAM_CIGAR_SHIFT|(o))
#define bam_cigar_type(o) (BAM_CIGAR_TYPE>>((o)<<1)&3) // bit 1: consume query; bit 2: consume reference
/* end */


#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place) @discussion x will be modified.
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

static void kstring_copy(kstring_t *dest, const kstring_t *other)
{
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

static void kseq_stack_destroy(kseq_t *seq)
{
    free(seq->name.s);
    free(seq->comment.s);
    free(seq->seq.s);
    free(seq->qual.s);
}

static void kseq_copy(kseq_t *dest, const kseq_t *seq)
{
    dest->f = NULL;
    kstring_copy(&dest->name, &seq->name);
    kstring_copy(&dest->comment, &seq->comment);
    kstring_copy(&dest->seq, &seq->seq);
    kstring_copy(&dest->qual, &seq->qual);
}

kseq_v read_seqs(kseq_t *seq,
                 size_t n_wanted)
{
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
    uint32_t nm;
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
    int m;            /* Number of residue tyes */
    int8_t *mat;      /* Scoring matrix */
    int32_t n_keep;   /* Number of alignments to keep */
    int32_t max_drop; /* Maximum drop from best alignment score */
} align_config_t;

static aln_v align_read(const kseq_t *read,
                        const kseq_v targets,
                        const align_config_t *conf)
{
    kseq_t *r;
    const int32_t read_len = read->seq.l;

    aln_v result;
    kv_init(result);
    kv_resize(aln_t, result, kv_size(targets));

    uint8_t *read_num = calloc(read_len, sizeof(uint8_t));

    for(size_t k = 0; k < read_len; ++k)
        read_num[k] = conf->table[(int)read->seq.s[k]];

    // Align to each target
    kswq_t *qry = NULL;
    for(size_t j = 0; j < kv_size(targets); j++) {
        // Encode target
        r = &kv_A(targets, j);
        uint8_t *ref_num = calloc(r->seq.l, sizeof(uint8_t));
        for(size_t k = 0; k < r->seq.l; ++k)
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
                   50, /* TODO: Magic number - band width */
                   &aln.n_cigar,
                   &aln.cigar);

        aln.nm = 0;
        size_t qi = aln.loc.qb, ri = aln.loc.tb;
        for(size_t k = 0; k < aln.n_cigar; k++) {
            const int32_t oplen = bam_cigar_oplen(aln.cigar[k]),
                          optype = bam_cigar_type(aln.cigar[k]);

            if(optype & 3) { // consumes both - check for mismatches
                for(size_t j = 0; j < oplen; j++) {
                    if(UNLIKELY(read_num[qi + j] != ref_num[ri + j]))
                        aln.nm++;
                }
            } else {
                aln.nm += oplen;
            }
            if(optype & 1) qi += oplen;
            if(optype & 2) ri += oplen;
        }

        kv_push(aln_t, result, aln);
        free(ref_num);
    }
    free(qry);
    free(read_num);
    ks_introsort(dec_score, kv_size(result), result.a);

    return result;
}

static void write_sam_records(kstring_t *str,
                              const kseq_t *read,
                              const aln_v result,
                              const kseq_v ref_seqs,
                              const char *read_group_id,
                              const int32_t n_keep,
                              const int32_t max_drop)
{
    if(kv_size(result) == 0)
        return;

    const int32_t min_score = kv_A(result, 0).loc.score - max_drop;

    /* Alignments are sorted by decreasing score */
    for(size_t i = 0; i < kv_size(result) && (n_keep <= 0 || i < n_keep); i++) {
        aln_t a = kv_A(result, i);

        if(a.loc.score < min_score)
            break;

        ksprintf(str, "%s\t%d\t", read->name.s,
                 i == 0 ? 0 : 256); // Secondary
        ksprintf(str, "%s\t%d\t%d\t",
                 kv_A(ref_seqs, a.target_idx).name.s, /* Reference */
                 a.loc.tb + 1,                        /* POS */
                 40);                                 /* MAPQ */
        if(a.loc.qb)
            ksprintf(str, "%dS", a.loc.qb);
        for(size_t c = 0; c < a.n_cigar; c++) {
            int32_t letter = 0xf&*(a.cigar + c);
            int32_t length = (0xfffffff0&*(a.cigar + c))>>4;
            ksprintf(str, "%d", length);
            if(letter == 0) ksprintf(str, "M");
            else if(letter == 1) ksprintf(str, "I");
            else ksprintf(str, "D");
        }

        if(a.loc.qe + 1 != read->seq.l)
            ksprintf(str, "%luS", read->seq.l - a.loc.qe - 1);

        ksprintf(str, "\t*\t0\t0\t");
        ksprintf(str, "%s\t", i > 0 ? "*" : read->seq.s);
        if(read->qual.s && i == 0)
            ksprintf(str, "%s", read->qual.s);
        else
            ksprintf(str, "*");

        ksprintf(str, "\tAS:i:%d\tNM:i:%d", a.loc.score, a.nm);
        if(read_group_id)
            ksprintf(str, "\tRG:Z:%s", read_group_id);
        kputs("\n", str);
    }
}


typedef struct {
    size_t start;
    size_t step;
    size_t n;
    kseq_v ref_seqs;
    kseq_v reads;
    kstring_t *sams;
    align_config_t *config;
    const char *read_group_id;
} worker_t;


static void *worker(void *data)
{
    worker_t *w = (worker_t *)data;
    for(size_t i = w->start; i < w->n; i+= w->step) {
        kseq_t *s = &kv_A(w->reads, i);
        aln_v result = align_read(s,
                                  w->ref_seqs,
                                  w->config);

        kstring_t str = { 0, 0, NULL };

        write_sam_records(&str,
                          s,
                          result,
                          w->ref_seqs,
                          w->read_group_id,
                          w->config->n_keep,
                          w->config->max_drop);

        w->sams[i] = str;

        for(size_t j = 0; j < kv_size(result); j++)
            free(kv_A(result, j).cigar);
        kv_destroy(result);
        kseq_stack_destroy(s);
    }

    return 0;
}

void align_reads(const char *ref_path,
                 const char *qry_path,
                 const char *output_path,
                 const int32_t match,       /* 2 */
                 const int32_t mismatch,    /* 2 */
                 const int32_t gap_o,       /* 3 */
                 const int32_t gap_e,       /* 1 */
                 const uint8_t n_threads,   /* 1 */
                 const int32_t n_keep,
                 const int32_t max_drop,
                 const char *read_group,
                 const char *read_group_id)
{
    gzFile read_fp, ref_fp;
    FILE *out_fp;
    int32_t j, k, l;
    const int m = 5;
    kseq_t *seq;
    int8_t *mat = (int8_t *)calloc(25, sizeof(int8_t));

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
    for(l = k = 0; LIKELY(l < 4); ++l) {
        for(j = 0; LIKELY(j < 4); ++j) mat[k++] = l == j ? match : -mismatch;	/* weight_match : -weight_mismatch */
        mat[k++] = 0; // ambiguous base
    }
    for(j = 0; LIKELY(j < 5); ++j) mat[k++] = 0;


    // Read reference sequences
    ref_fp = gzopen(ref_path, "r");
    assert(ref_fp != NULL && "Failed to open reference");
    seq = kseq_init(ref_fp);
    kseq_v ref_seqs;
    ref_seqs = read_seqs(seq, 0);
    kseq_destroy(seq);
    gzclose(ref_fp);

    fprintf(stderr, "[sw_align] Read %lu references\n",
            kv_size(ref_seqs));

    // Print SAM header
    out_fp = fopen(output_path, "w");
    fprintf(out_fp, "@HD\tVN:1.4\tSO:unsorted\n");
    for(size_t i = 0; i < kv_size(ref_seqs); i++) {
        seq = &kv_A(ref_seqs, i);
        fprintf(out_fp, "@SQ\tSN:%s\tLN:%d\n",
                seq->name.s, (int32_t)seq->seq.l);
    }
    if(read_group) {
        fputs(read_group, out_fp);
        fputc('\n', out_fp);
    }

    align_config_t conf;
    conf.gap_o = gap_o;
    conf.gap_e = gap_e;
    conf.m = m;
    conf.table = table;
    conf.mat = mat;
    conf.n_keep = n_keep;
    conf.max_drop = max_drop;

    read_fp = gzopen(qry_path, "r");
    assert(read_fp != NULL && "Failed to open query");
    size_t count = 0;
    seq = kseq_init(read_fp);
    while(true) {
        kseq_v reads = read_seqs(seq, 5000 * n_threads);
        const size_t n_reads = kv_size(reads);
        if(!n_reads) {
            break;
        }

        worker_t *w = calloc(n_threads, sizeof(worker_t));
        kstring_t *sams  = calloc(n_reads, sizeof(kstring_t));
        for(size_t i = 0; i < n_threads; i++) {
            w[i].start = i;
            w[i].n = n_reads;
            w[i].step = n_threads;
            w[i].ref_seqs = ref_seqs;
            w[i].reads = reads;
            w[i].sams = sams;
            w[i].config = &conf;
            w[i].read_group_id = read_group_id;
        }

        if(n_threads == 1) {
            worker(w);
        } else {
            pthread_t *tid = calloc(n_threads, sizeof(pthread_t));
            for(size_t i = 0; i < n_threads; ++i)
                pthread_create(&tid[i], 0, worker, &w[i]);
            for(size_t i = 0; i < n_threads; ++i)
                pthread_join(tid[i], 0);
        }
        free(w);

        for(size_t i = 0; i < n_reads; i++) {
            fputs(sams[i].s, out_fp);
            free(sams[i].s);
        }
        free(sams);
        count += n_reads;
        kv_destroy(reads);
    }
    kseq_destroy(seq);
    fprintf(stderr, "[sw_align] Aligned %lu reads\n", count);

    // Clean up reference sequences
    kvi_destroy(kseq_stack_destroy, ref_seqs);

    gzclose(read_fp);
    fclose(out_fp);
    free(mat);
}
