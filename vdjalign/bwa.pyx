# distutils: language = c
# distutils: include_dirs = BWA_INCLUDE
# distutils: extra_compile_args = BWA_FLAGS
# distutils: sources = BWA_SRC
# distutils: libraries = BWA_LIBS

from libc.stdint cimport uint64_t, uint8_t, int8_t, int64_t, uint32_t
from libc.stdlib cimport free

cdef extern from "bntseq.h":
    ctypedef struct bntseq_t:
        pass

cdef extern from "bwt.h":
    ctypedef uint64_t bwtint_t

    ctypedef struct bwt_t:
        pass

cdef extern from "bwa.h":
    ctypedef struct bwaidx_t:
        bwt_t *bwt
        bntseq_t *bns
        uint8_t *pac

    bwaidx_t* bwa_idx_load(const char *hint, int which)
    void bwa_idx_destroy(bwaidx_t *idx)


cdef extern from "bwamem.h":
    ctypedef struct mem_opt_t:
        pass

    mem_opt_t *mem_opt_init()
    void mem_fill_scmat(int a, int b, int8_t mat[25])

    ctypedef struct mem_alnreg_t:
        pass

    ctypedef struct mem_alnreg_v:
        size_t n
        size_t m
        mem_alnreg_t *a

    ctypedef struct mem_aln_t:
        int64_t pos
        int rid # Reference id (<0 = unmapped)
        int flag
        int n_cigar
        uint32_t *cigar # CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234

        int score

    mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq)

    mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar)


def align_sequence(index_path, bytes seq):
    cdef mem_opt_t *opt = mem_opt_init()
    cdef char* s = seq
    cdef int l_seq = len(seq)

    cdef bwaidx_t *idx = NULL
    cdef int BWA_IDX_ALL = 0x7

    cdef mem_alnreg_v alignments
    cdef mem_aln_t alignment
    try:
        idx = bwa_idx_load(index_path, BWA_IDX_ALL)
        if idx == NULL:
            raise IOError("Could not open " + seq)
        alignments = mem_align1(opt, idx.bwt, idx.bns, idx.pac, l_seq, s)
        for i in range(alignments.n):
            alignment = mem_reg2aln(opt, idx.bns, idx.pac, l_seq, s, &alignments.a[i])
            print alignment.rid, alignment.pos
    finally:
        if idx != NULL:
            bwa_idx_destroy(idx)
        free(opt)

