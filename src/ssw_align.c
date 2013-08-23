/*  main.c
 *  Created by Mengyao Zhao on 06/23/11.
 *	Version 0.1.4
 *  Last revision by Mengyao Zhao on 07/31/12.
 */

#include <stdlib.h>
#include <stdint.h>
#include <emmintrin.h>
#include <zlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "ssw.h"
#include "kseq.h"
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

int ascore_lt_dec(const void* va, const void* vb) {
    const s_align* a = *((const s_align**) va);
    const s_align* b = *((const s_align**) vb);
    if(a->score1 == b->score1) return 0;
    if(a->score1 < b->score1) return 1;
    return -1;
}

void ssw_write_sam (gzFile* dest,
                    const s_align* a,
                    const kseq_t* ref_seq,
                    const kseq_t* read,
                    const char* read_seq,
                    const int8_t* table,
                    const int32_t flag) {

        gzprintf(dest, "%s\t", read->name.s);
        if (a->score1 == 0) gzprintf(dest, "4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n");
        else {
            int32_t c, l = a->read_end1 - a->read_begin1 + 1, qb = a->ref_begin1, pb = a->read_begin1, p;
            uint32_t mapq = -4.343 * log(1 - (double)abs(a->score1 - a->score2)/(double)a->score1);
            mapq = (uint32_t) (mapq + 4.99);
            mapq = mapq < 254 ? mapq : 254;
            gzprintf(dest, "%d\t", flag);
            gzprintf(dest, "%s\t%d\t%d\t", ref_seq->name.s, a->ref_begin1 + 1, mapq);
            if (a->read_begin1)
                gzprintf(dest, "%dS", a->read_begin1);

            for (c = 0; c < a->cigarLen; ++c) {
                int32_t letter = 0xf&*(a->cigar + c);
                int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
                gzprintf(dest, "%d", length);
                if (letter == 0) gzprintf(dest, "M");
                else if (letter == 1) gzprintf(dest, "I");
                else gzprintf(dest, "D");
            }

            if (a->read_end1 + 1 != read->seq.l)
                gzprintf(dest, "%dS", read->seq.l - a->read_end1 - 1);

            gzprintf(dest, "\t*\t0\t0\t");
            if (!(flag & 256)) {
                gzprintf(dest, read->seq.s);
            } else {
                gzprintf(dest, "*");
            }
            gzprintf(dest, "\t");
            if (read->qual.s && !(flag & 256)) {
                gzprintf(dest, read->qual.s);
            } else gzprintf(dest, "*");
            gzprintf(dest, "\tAS:i:%d", a->score1);
            mapq = 0;	// counter of difference
            for (c = 0; c < a->cigarLen; ++c) {
                int32_t letter = 0xf&*(a->cigar + c);
                int32_t length = (0xfffffff0&*(a->cigar + c))>>4;
                if (letter == 0) {
                    for (p = 0; p < length; ++p){
                        if (table[(int)*(ref_seq->seq.s + qb)] != table[(int)*(read_seq + pb)]) ++mapq;
                        ++qb;
                        ++pb;
                    }
                } else if (letter == 1) {
                    pb += length;
                    mapq += length;
                } else {
                    qb += length;
                    mapq += length;
                }
            }
            gzprintf(dest,"\tNM:i:%d\t", mapq);
            if (a->score2 > 0) gzprintf(dest, "ZS:i:%d\n", a->score2);
            else gzprintf(dest, "\n");
        }
}

void align_reads (const char* ref_path,
                  const char* qry_path,
                  const char* output_path,
                  const int32_t match, /* 2 */
                  const int32_t mismatch, /* 2 */
                  const int32_t gap_open, /* 3 */
                  const int32_t gap_extension) { /* 1 */
    clock_t start, end;
    float cpu_time;
    gzFile read_fp, ref_fp, out_fp;
    kseq_t *read_seq, *ref_seq;
    int32_t l, m, k, n = 5, s1 = 1048575, s2 = 128, filter = 0;
    kroundup32(s1);
    int8_t* mata = (int8_t*)calloc(25, sizeof(int8_t)), *mat = mata;
    int8_t* ref_num = (int8_t*)malloc(s1);
    int8_t* num = (int8_t*)malloc(s2);

    /* This table is used to transform nucleotide letters into numbers. */
    int8_t table[128] = {
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
        for (m = 0; LIKELY(m < 4); ++m) mata[k++] = l == m ? match : -mismatch;	/* weight_match : -weight_mismatch */
        mata[k++] = 0; // ambiguous base
    }
    for (m = 0; LIKELY(m < 5); ++m) mata[k++] = 0;

    read_fp = gzopen(qry_path, "r");

    // Read reference sequences
    kvec_t(kseq_t*) ref_seqs;
    kv_init(ref_seqs);

    ref_fp = gzopen(ref_path, "r");
    ref_seq = kseq_init(ref_fp);
    while ((l = kseq_read(ref_seq)) >= 0) {
        kv_push(kseq_t*, ref_seqs, ref_seq);
        ref_seq = kseq_init(ref_fp);
    }
    kseq_destroy(ref_seq); // One extra allocated
    gzclose(ref_fp);

    // Print SAM header
    out_fp = gzopen(output_path, "w0");
    gzprintf(out_fp, "@HD\tVN:1.4\tSO:queryname\n");
    for(l = 0; l < kv_size(ref_seqs); l++) {
        ref_seq = kv_A(ref_seqs, l);
        gzprintf(out_fp, "@SQ\tSN:%s\tLN:%d\n", ref_seq->name.s, (int32_t)ref_seq->seq.l);
    }

    read_seq = kseq_init(read_fp);

    // alignment
    start = clock();
    while (kseq_read(read_seq) >= 0) {
        s_profile* p;
        const int32_t readLen = read_seq->seq.l;
        const int32_t maskLen = readLen / 2;

        while (readLen >= s2) {
            ++s2;
            kroundup32(s2);
            num = (int8_t*)realloc(num, s2);
        }
        for (m = 0; m < readLen; ++m) num[m] = table[(int)read_seq->seq.s[m]];
        p = ssw_init(num, readLen, mat, n, 2);

        kvec_t(s_align*) alignments;
        kv_init(alignments);
        kv_resize(s_align*, alignments, kv_size(ref_seqs));

        for(l = 0; l < kv_size(ref_seqs); ++l) {
            ref_seq = kv_A(ref_seqs, l);
            int32_t refLen = ref_seq->seq.l;
            int8_t flag = 0;
            while (refLen > s1) {
                ++s1;
                kroundup32(s1);
                ref_num = (int8_t*)realloc(ref_num, s1);
            }
            for (m = 0; m < refLen; ++m) ref_num[m] = table[(int)ref_seq->seq.s[m]];
            flag = 2;
            s_align* result = ssw_align (p, ref_num, refLen, gap_open, gap_extension, flag, filter, 0, maskLen);
            kv_push(s_align*, alignments, result);
        }

        // Descending alignment score
        qsort(alignments.a,
              kv_size(alignments),
              sizeof(s_align*),
              ascore_lt_dec);

        for(l = 0; l < kv_size(alignments); l++) {
            ssw_write_sam(out_fp, kv_A(alignments, l),
                          ref_seq, read_seq, read_seq->seq.s, table,
                          l == 0 ? 0 : 256);
            align_destroy(kv_A(alignments, l));
        }
        kv_destroy(alignments);

        init_destroy(p);
    }


    end = clock();
    cpu_time = ((float) (end - start)) / CLOCKS_PER_SEC;
    fprintf(stderr, "CPU time: %f seconds\n", cpu_time);

    kseq_destroy(read_seq);
    // Clean up reference sequences
    for(l = 0; l < kv_size(ref_seqs); ++l)
        free(kv_A(ref_seqs, l));
    kv_destroy(ref_seqs);

    gzclose(read_fp);
    gzclose(out_fp);
    free(num);
    free(ref_num);
    free(mata);
}
