#include "ig_align.h"

#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "kseq.h"
#include "ksort.h"
#include "kstring.h"
#include "ksw.h"
#include "kvec.h"
#include <zlib.h>

#define xstr(a) str(a)
#define str(a) #a

/* Cigar operations */
/**
 * Describing how CIGAR operation/length is packed in a 32-bit integer.
 */
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK ((1 << BAM_CIGAR_SHIFT) - 1)

#define BAM_CIGAR_STR "MIDNSHP=XB"
#define BAM_CIGAR_TYPE 0x3C1A7

#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c) >> BAM_CIGAR_SHIFT)
#define bam_cigar_opchr(c) (BAM_CIGAR_STR[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l) << BAM_CIGAR_SHIFT | (o))
#define bam_cigar_type(o)                                                      \
  (BAM_CIGAR_TYPE >> ((o) << 1) &                                              \
   3) // bit 1: consume query; bit 2: consume reference
/* end */

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/*! @function
  @abstract  Round an integer to the next closest power-2 integer.
  @param  x  integer to be rounded (in place) @discussion x will be modified.
  */
#ifndef kroundup32
#define kroundup32(x)                                                          \
  (--(x), (x) |= (x) >> 1, (x) |= (x) >> 2, (x) |= (x) >> 4, (x) |= (x) >> 8,  \
   (x) |= (x) >> 16, ++(x))
#endif

KSEQ_INIT(gzFile, gzread)

typedef kvec_t(kseq_t) kseq_v;

/* Destroy a vector of pointers, calling f on each item */
#define kvp_destroy(f, v)                                                      \
  for (size_t __kpd = 0; __kpd < kv_size(v); ++__kpd)                          \
    f(kv_A(v, __kpd));                                                         \
  kv_destroy(v)

#define kvi_destroy(f, v)                                                      \
  for (size_t __kpd = 0; __kpd < kv_size(v); ++__kpd)                          \
    f(&kv_A(v, __kpd));                                                        \
  kv_destroy(v)

static void kstring_copy(kstring_t *dest, const kstring_t *other) {
  dest->l = other->l;
  dest->m = dest->l + 1;
  if (dest->l > 0) {
    assert(other->s != NULL);
    dest->s = malloc(dest->l + 1);
    memcpy(dest->s, other->s, dest->l);
    dest->s[dest->l] = 0; // null-terminate
  } else {
    dest->s = NULL;
  }
}

static void kseq_stack_destroy(kseq_t *seq) {
  free(seq->name.s);
  free(seq->comment.s);
  free(seq->seq.s);
  free(seq->qual.s);
}

static void kseq_copy(kseq_t *dest, const kseq_t *seq) {
  dest->f = NULL;
  kstring_copy(&dest->name, &seq->name);
  kstring_copy(&dest->comment, &seq->comment);
  kstring_copy(&dest->seq, &seq->seq);
  kstring_copy(&dest->qual, &seq->qual);
}

static kseq_v read_seqs(kseq_t *seq, size_t n_wanted) {
  kseq_v result;
  kv_init(result);
  for (size_t i = 0; i < n_wanted || n_wanted == 0; i++) {
    if (kseq_read(seq) <= 0)
      break;
    kseq_t s;
    kseq_copy(&s, seq);
    kv_push(kseq_t, result, s);
  }
  return result;
}

typedef struct {
  char *target_name;
  kswr_t loc;
  uint32_t *cigar;
  int n_cigar;
  uint32_t nm;
} aln_t;
typedef kvec_t(aln_t) aln_v;

#define __aln_score_lt(a, b) ((a).loc.score > (b).loc.score)
KSORT_INIT(cdec_score, aln_t, __aln_score_lt)

typedef struct { aln_v alignments; } aln_task_t;
typedef kvec_t(aln_task_t) aln_task_v;

typedef struct {
  int32_t gap_o; /* 3 */
  int32_t gap_e; /* 1 */
  uint8_t *table;
  int m;       /* Number of residue tyes */
  int8_t *mat; /* Scoring matrix */
  int max_drop;
  int min_score;
  unsigned bandwidth;
} align_config_t;

static aln_t align_read_against_one(kseq_t *target, const int read_len,
                                    uint8_t *read_num, kswq_t **qry,
                                    const align_config_t *conf,
                                    const int min_score) {
  uint8_t *ref_num = calloc(target->seq.l, sizeof(uint8_t));
  for (size_t k = 0; k < target->seq.l; ++k)
    ref_num[k] = conf->table[(int)target->seq.s[k]];

  aln_t aln;
  aln.cigar = NULL;
  aln.loc = ksw_align(read_len, read_num, target->seq.l, ref_num, conf->m,
                      conf->mat, conf->gap_o, conf->gap_e, KSW_XSTART, qry);

  aln.target_name = target->name.s;

  if (aln.loc.score < min_score) {
    free(ref_num);
    return aln;
  }

  ksw_global(aln.loc.qe - aln.loc.qb + 1, &read_num[aln.loc.qb],
             aln.loc.te - aln.loc.tb + 1, &ref_num[aln.loc.tb], conf->m,
             conf->mat, conf->gap_o, conf->gap_e, conf->bandwidth, &aln.n_cigar,
             &aln.cigar);

  aln.nm = 0;
  size_t qi = aln.loc.qb, ri = aln.loc.tb;
  for (int k = 0; k < aln.n_cigar; k++) {
    const int32_t oplen = bam_cigar_oplen(aln.cigar[k]),
                  optype = bam_cigar_type(aln.cigar[k]);

    if (optype & 3) { // consumes both - check for mismatches
      for (int j = 0; j < oplen; j++) {
        if (UNLIKELY(read_num[qi + j] != ref_num[ri + j]))
          aln.nm++;
      }
    } else {
      aln.nm += oplen;
    }
    if (optype & 1)
      qi += oplen;
    if (optype & 2)
      ri += oplen;
  }

  free(ref_num);

  /* size_t cigar_len = aln.loc.qb; */
  /* for (int c = 0; c < aln.n_cigar; c++) { */
  /*   int32_t length = (0xfffffff0 & *(aln.cigar + c)) >> 4; */
  /*   cigar_len += length; */
  /* } */
  /* cigar_len += read_len - aln.loc.qe - 1; */
  /* if(cigar_len != (size_t)read_len) { */
  /*   /\* printf("[ig_align] Error: cigar length (score %d) not equal to read length for XXX (target %s): %zu vs %d\n", aln.loc.score, target->name.s, cigar_len, read_len); *\/ */
  /*   // NOTE: */
  /*   //   It is *really* *fucking* *scary* that it's spitting out cigars that are not the same length as the query sequence. */
  /*   //   Nonetheless, fixing it seems to involve delving into the depths of ksw_align() and ksw_global(), which would be very time consuming, and the length discrepancy seems to ony appear in very poor matches. */
  /*   //   I.e., poor enough that we will subsequently ignore them in partis/python/waterer.py, so it seems to not screw anything up downstream to just set the length-discrepant matches' scores to zero, such that ig-sw doesn't write them to its sam output. */
  /*   //   Note also that it is not always the lowest- or highest-scoring matches that have discrepant lengths (i.e. setting their scores to zero promotes matches swith poorer scores, but which do not have discrepant lengths. */
  /*   /\* aln.loc.score = 0; *\/ */
  /*   aln.cigar = NULL; */
  /* } */

  return aln;
}

/* Returns total size after dropping low scores */
void drop_low_scores(aln_v *vec, int offset, int max_drop) {
  const int size = kv_size(*vec);
  ks_introsort(cdec_score, size - offset, vec->a + offset);
  const int min_score = kv_A(*vec, offset).loc.score - max_drop;
  for (int i = offset; i < size; i++) {
    if (kv_A(*vec, i).loc.score < min_score) {
      vec->n = i;

      /* Free remaining */
      for (int j = i; j < size; j++)
        free(kv_A(*vec, j).cigar);
      return;
    }
  }
}

static aln_v align_read(const kseq_t *read, const kseq_v targets,
                        const size_t n_extra_targets,
                        const kseq_v *extra_targets,
                        const align_config_t *conf) {
  kseq_t *r;
  const int32_t read_len = read->seq.l;

  aln_v result;
  kv_init(result);
  kv_resize(aln_t, result, kv_size(targets));

  uint8_t *read_num = calloc(read_len, sizeof(uint8_t));  // vector of integers encoding the read sequence (i.e. one integer for each character in the read)

  for (int k = 0; k < read_len; ++k)
    read_num[k] = conf->table[(int)read->seq.s[k]];

  // first align entire read against first batch of targets (i.e. the v genes)
  kswq_t *qry = NULL;
  int min_score = -1000;
  int max_score = 0;
  for (size_t j = 0; j < kv_size(targets); j++) {
    // Encode target
    r = &kv_A(targets, j);
    aln_t aln =
      align_read_against_one(r, read_len, read_num, &qry, conf, min_score);
    if (aln.cigar != NULL) {
      max_score = aln.loc.score > max_score ? aln.loc.score : max_score;
      min_score = (aln.loc.score - conf->max_drop) > min_score
                      ? (aln.loc.score - conf->max_drop)
                      : min_score;
      kv_push(aln_t, result, aln);
    }
  }

  /* If no alignments to the first set of targets reached the minimum score,
   * abort.
   */
  if (max_score < conf->min_score) {
    // kv_size returns the n field of a kvec_t, which is a size_t.
    for (size_t i = 0; i < kv_size(result); i++)
      free(kv_A(result, i).cigar);
    kv_size(result) = 0;
    free(qry);
    free(read_num);
    return result;
  }

  drop_low_scores(&result, 0, conf->max_drop);

  // Extra references - qe points to the exact end of the sequence
  int qend = kv_A(result, 0).loc.qe + 1;  // + 1 converts to python slice conventions
  int read_len_trunc = read_len - qend;
  uint8_t *read_num_trunc = read_num + qend;  // truncate the encoded read sequence by incremeting the start position (well, pointer to)

  free(qry);
  qry = NULL;

  // then align any part of the read that remains after truncation against any targets in <extra_targets> (i.e., presumably, first against j [after which we truncate the j match] and then against d)
  if (read_len_trunc > 2) {
    for (size_t i = 0; i < n_extra_targets; i++) {
      const size_t idx = n_extra_targets - i - 1;
      min_score = -1000;
      const size_t init_count = kv_size(result);
      for (size_t j = 0; j < kv_size(extra_targets[idx]); j++) {
        r = &kv_A(extra_targets[idx], j);
        aln_t aln = align_read_against_one(r, read_len_trunc, read_num_trunc,
                                           &qry, conf, min_score);

        if (aln.cigar != NULL) {
          min_score = (aln.loc.score - conf->max_drop) > min_score
                          ? (aln.loc.score - conf->max_drop)
                          : min_score;
          aln.loc.qb += qend;
          aln.loc.qe += qend;
          kv_push(aln_t, result, aln);
        }
      }
      drop_low_scores(&result, init_count, conf->max_drop);

      // truncate read/query sequence by adjusting the length that we pass to the aligner (i.e. we first do j, then effectively remove the j portion of the read by decreasing <read_len_trunc>)
      read_len_trunc = kv_A(result, init_count).loc.qb - qend;
      free(qry);
      qry = NULL;
    }
  }

  free(qry);
  free(read_num);

  return result;
}

static void write_sam_records(kstring_t *str, const kseq_t *read,
                              const aln_v result, const char *read_group_id) {
  if (kv_size(result) == 0)
    return;

  /* Alignments are sorted by decreasing score */
  bool wrote_primary_stuff = false;  // can't any more use i==0 criterion, since the first match may be one with discordant cigar and read lengths
  aln_t *tmp_aln = NULL;  // pointer to first alignment, so we can write a sensible target name to the dummy line
  for (size_t i = 0; i < kv_size(result); i++) {
    aln_t a = kv_A(result, i);
    if (i == 0)
      tmp_aln = &a;

    kstring_t tmpstr = {0, 0, NULL};  // temporary, so we can avoid writing to <str> until we know if the cigar and read are the same length
    ksprintf(&tmpstr, "%s\t%d\t", read->name.s, !wrote_primary_stuff ? 0 : 256); /* Secondary */
    ksprintf(&tmpstr, "%s\t%d\t%d\t", a.target_name,               /* Reference */
             a.loc.tb + 1,                                     /* POS */
             40);                                              /* MAPQ */
    if (a.loc.qb)
      ksprintf(&tmpstr, "%dS", a.loc.qb);
    size_t match_length = 0;
    for (int c = 0; c < a.n_cigar; c++) {
      int32_t letter = 0xf & *(a.cigar + c);
      int32_t length = (0xfffffff0 & *(a.cigar + c)) >> 4;
      ksprintf(&tmpstr, "%d", length);
      if (letter == 0) {
        ksprintf(&tmpstr, "M");
	match_length += length;
      } else if (letter == 1) {
        ksprintf(&tmpstr, "I");
	match_length += length;
      } else {
        ksprintf(&tmpstr, "D");
      }
    }

    if(a.loc.qe - a.loc.qb != (int)match_length - 1) {
      /* fprintf(stderr, "[ig_align] Error: match length mismatch for query %s score %d target %s: qe - qb = %d - %d = %d != match_length - 1 = %zu\n", read->name.s, a.loc.score, a.target_name, a.loc.qe, a.loc.qb, a.loc.qe - a.loc.qb, match_length - 1); */
      /* fprintf(stderr, "        %s    %s\n", tmpstr.s, read->seq.s); */
      /* // assert(0); */
      continue;
    }

    ksprintf(str, "%s", tmpstr.s);

    if (a.loc.qe + 1 != (int)read->seq.l)
      ksprintf(str, "%luS", read->seq.l - a.loc.qe - 1);

    ksprintf(str, "\t*\t0\t0\t");
    ksprintf(str, "%s\t", wrote_primary_stuff ? "*" : read->seq.s);
    if (read->qual.s && !wrote_primary_stuff)
      ksprintf(str, "%s", read->qual.s);
    else
      ksprintf(str, "*");

    ksprintf(str, "\tAS:i:%d\tNM:i:%d", a.loc.score, a.nm);
    if (read_group_id)
      ksprintf(str, "\tRG:Z:%s", read_group_id);
    kputs("\n", str);
    wrote_primary_stuff = true;
  }

  if (!wrote_primary_stuff) {  // no matches whatsoever made it through (probably due to lots of cigar/read length discrepancies) -- write a dummy line so the code reading it knows we didn't just lose the query
    ksprintf(str, "%s\t%d\t%s\t%d\t%d\t%dS\t*\t0\t0\t%s\t*\tAS:i:%d\tNM:i:%d",
	     read->name.s, 0, tmp_aln ? tmp_aln->target_name : "xxx", 1, 999, (int)read->seq.l, read->seq.s, 0, 0);
    if (read_group_id)
      ksprintf(str, "\tRG:Z:%s", read_group_id);
    kputs("\n", str);
  }
}

typedef struct {
  size_t start;
  size_t step;
  size_t n;
  kseq_v ref_seqs;
  int n_extra_refs;
  kseq_v *extra_ref_seqs;
  kseq_v reads;
  kstring_t *sams;
  align_config_t *config;
  const char *read_group_id;
} worker_t;

static void *worker(void *data) {
  worker_t *w = (worker_t *)data;
  for (size_t i = w->start; i < w->n; i += w->step) {
    kseq_t *s = &kv_A(w->reads, i);
    aln_v result = align_read(s, w->ref_seqs, w->n_extra_refs,
                              w->extra_ref_seqs, w->config);

    kstring_t str = {0, 0, NULL};

    write_sam_records(&str, s, result, w->read_group_id);

    w->sams[i] = str;

    for (size_t j = 0; j < kv_size(result); j++)
      free(kv_A(result, j).cigar);
    kv_destroy(result);
    kseq_stack_destroy(s);
  }

  return 0;
}

void ig_align_reads(const char *ref_path, const uint8_t n_extra_refs,
                    const char **extra_ref_paths, const char *qry_path,
                    const char *output_path, const int32_t match, /* 2 */
                    const int32_t mismatch,                       /* 2 */
                    const int32_t gap_o,                          /* 3 */
                    const int32_t gap_e,                          /* 1 */
                    const unsigned max_drop,                      /* 1000 */
                    const int min_score,                          /* 0 */
                    const unsigned bandwidth,                     /* 150 */
                    const uint8_t n_threads,                      /* 1 */
                    const char *read_group, const char *read_group_id) {
  gzFile read_fp, ref_fp;
  FILE *out_fp;
  int32_t j, k, l;
  const int m = 5;
  kseq_t *seq;
  int8_t *mat = (int8_t *)calloc(25, sizeof(int8_t));

  /* This table is used to transform nucleotide letters into numbers. */
  uint8_t table[128] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                        4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4,
                        4, 4, 4, 4, 4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                        4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                        4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

  // initialize scoring matrix for genome sequences
  for (l = k = 0; LIKELY(l < 4); ++l) {
    for (j = 0; LIKELY(j < 4); ++j)
      mat[k++] =
          l == j ? match : -mismatch; /* weight_match : -weight_mismatch */
    mat[k++] = 0;                     // ambiguous base
  }
  for (j = 0; LIKELY(j < 5); ++j)
    mat[k++] = 0;

  // Read reference sequences
  ref_fp = gzopen(ref_path, "r");
  if(ref_fp == NULL) {
    fprintf(stderr, "Failed to open reference %s\n", ref_path);
    assert(0);
  }
  seq = kseq_init(ref_fp);
  kseq_v ref_seqs;
  ref_seqs = read_seqs(seq, 0);
  kseq_destroy(seq);
  gzclose(ref_fp);

  fprintf(stderr, "[ig_align] Read %lu references\n", kv_size(ref_seqs));

  kseq_v *extra_ref_seqs = malloc(sizeof(kseq_v) * n_extra_refs);

  for (size_t i = 0; i < n_extra_refs; i++) {
    ref_fp = gzopen(extra_ref_paths[i], "r");
    assert(ref_fp != NULL && "Failed to open reference");
    seq = kseq_init(ref_fp);
    extra_ref_seqs[i] = read_seqs(seq, 0);
    kseq_destroy(seq);
    gzclose(ref_fp);
    fprintf(stderr, "[ig_align] Read %lu extra references from %s\n",
            kv_size(extra_ref_seqs[i]), extra_ref_paths[i]);
  }

  // Print SAM header
  out_fp = fopen(output_path, "w");
  fprintf(out_fp, "@HD\tVN:1.4\tSO:unsorted\n");
  fprintf(out_fp, "@PG\tID:ig_align\tPN:ig_align\tCL:match=%d,mismatch=%d,go=%"
                  "d,ge=%d\tVN:%s\n",
          match, mismatch, gap_o, gap_e, xstr(VDJALIGN_VERSION));
  for (size_t i = 0; i < kv_size(ref_seqs); i++) {
    seq = &kv_A(ref_seqs, i);
    fprintf(out_fp, "@SQ\tSN:%s\tLN:%d\n", seq->name.s, (int32_t)seq->seq.l);
  }
  for (size_t i = 0; i < n_extra_refs; i++) {
    for (size_t j = 0; j < kv_size(extra_ref_seqs[i]); j++) {
      seq = &kv_A(extra_ref_seqs[i], j);
      fprintf(out_fp, "@SQ\tSN:%s\tLN:%d\n", seq->name.s, (int32_t)seq->seq.l);
    }
  }
  if (read_group) {
    fputs(read_group, out_fp);
    fputc('\n', out_fp);
  }

  align_config_t conf;
  conf.gap_o = gap_o;
  conf.gap_e = gap_e;
  conf.max_drop = max_drop;
  conf.min_score = min_score;
  conf.m = m;
  conf.table = table;
  conf.mat = mat;
  conf.bandwidth = bandwidth;

  read_fp = gzopen(qry_path, "r");
  assert(read_fp != NULL && "Failed to open query");
  size_t count = 0;
  seq = kseq_init(read_fp);
  while (true) {
    kseq_v reads = read_seqs(seq, 5000 * n_threads);
    const size_t n_reads = kv_size(reads);
    if (!n_reads) {
      break;
    }

    worker_t *w = calloc(n_threads, sizeof(worker_t));
    kstring_t *sams = calloc(n_reads, sizeof(kstring_t));
    for (size_t i = 0; i < n_threads; i++) {
      w[i].start = i;
      w[i].n = n_reads;
      w[i].step = n_threads;
      w[i].ref_seqs = ref_seqs;
      w[i].n_extra_refs = n_extra_refs;
      w[i].extra_ref_seqs = extra_ref_seqs;
      w[i].reads = reads;
      w[i].sams = sams;
      w[i].config = &conf;
      w[i].read_group_id = read_group_id;
    }

    if (n_threads == 1) {
      worker(w);
    } else {
      pthread_t *tid = calloc(n_threads, sizeof(pthread_t));
      for (size_t i = 0; i < n_threads; ++i)
        pthread_create(&tid[i], 0, worker, &w[i]);
      for (size_t i = 0; i < n_threads; ++i)
        pthread_join(tid[i], 0);
    }
    free(w);

    for (size_t i = 0; i < n_reads; i++) {
      if (sams[i].s) {
        fputs(sams[i].s, out_fp);
        free(sams[i].s);
      }
    }
    free(sams);
    count += n_reads;
    kv_destroy(reads);
  }
  kseq_destroy(seq);
  fprintf(stderr, "[ig_align] Aligned %lu reads\n", count);

  // Clean up reference sequences
  kvi_destroy(kseq_stack_destroy, ref_seqs);

  // And extra reference sequences
  for (size_t i = 0; i < n_extra_refs; i++) {
    kvi_destroy(kseq_stack_destroy, extra_ref_seqs[i]);
  }
  free(extra_ref_seqs);

  gzclose(read_fp);
  fclose(out_fp);
  free(mat);
}
