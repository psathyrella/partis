# distutils: language = c
# distutils: include_dirs = BWA_INCLUDE
# distutils: extra_compile_args = BWA_FLAGS
# distutils: sources = BWA_SRC
# distutils: libraries = BWA_LIBS

from libc.stdint cimport uint64_t, uint8_t, int8_t, int64_t, uint32_t

cdef extern from "bntseq.h":
    cdef struct bntseq_t:
        pass

cdef extern from "bwt.h":
    ctypedef uint64_t bwtint_t

    cdef struct bwt_t:
        pass

cdef extern from "bwa.h":
    cdef struct bwaidx_t:
        bwt_t *bwt
        bntseq_t *bns
        uint8_t *pac

    bwaidx_t* bwa_idx_load(const char *hint, int which)
    void bwa_idx_destroy(bwaidx_t *idx)


cdef extern from "bwamem.h":
    cdef struct mem_opt_t:
        pass

    mem_opt_t *mem_opt_init()
    void mem_fill_scmat(int a, int b, int8_t mat[25])

    cdef struct mem_alnreg_t:
        pass

    cdef struct mem_alnreg_v:
        size_t n
        size_t m
        mem_alnreg_t *a

    cdef struct mem_aln_t:
        int64_t pos
        int rid # Reference id (<0 = unmapped)
        int flag
        int n_cigar
        uint32_t *cigar # CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234

        int score

    mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq)

    mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar)



print "hello world!"
