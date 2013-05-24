# distutils: language = c
# distutils: include_dirs = BWA_INCLUDE
# distutils: extra_compile_args = BWA_FLAGS
# distutils: sources = BWA_SRC
# distutils: libraries = BWA_LIBS

from libc.stdint cimport uint64_t, uint8_t, int8_t, int64_t, uint32_t, int32_t
from libc.stdlib cimport free, malloc

cdef extern from "bntseq.h":
    ctypedef struct bntann1_t:
        char *name
    ctypedef struct bntseq_t:
        int64_t l_pac
        int32_t n_seqs
        uint32_t seed
        bntann1_t *anns

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
        uint32_t is_rev, mapq, NM
        int n_cigar
        uint32_t *cigar # CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234

        int score

    mem_alnreg_v mem_align1(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq) nogil

    mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_seq, const char *seq, const mem_alnreg_t *ar)

cdef list parse_cigar(uint32_t *cigar, int n_cigar):
    cdef list result = []
    for i in xrange(n_cigar):
        result.append((cigar[i]>>4, "MIDSH"[cigar[i] & 0xF]))
    return result

cdef class BwaIndex:
    cdef bwaidx_t *idx
    cdef bytes path

    def __init__(self):
        raise ValueError("This class cannot be instantiated from Python.")

    def __dealloc__(self):
        if self.idx != NULL:
            bwa_idx_destroy(self.idx)

    def __repr__(self):
        return '<BwaIndex {0}>'.format(self.path)

    def align(self, bytes seq):
        cdef mem_opt_t *opt = mem_opt_init()
        cdef char* s = seq
        cdef int l_seq = len(seq)

        cdef int BWA_IDX_ALL = 0x7

        cdef mem_alnreg_v ar
        cdef mem_aln_t a
        try:
            with nogil:
                ar = mem_align1(opt, self.idx.bwt, self.idx.bns, self.idx.pac, l_seq, s)
            result = []
            for i in xrange(ar.n):
                a = mem_reg2aln(opt, self.idx.bns, self.idx.pac, l_seq, s, &ar.a[i])
                result.append({'pos': a.pos, 'rid': a.rid, 'flag': a.flag,
                    'reference': self.idx.bns.anns[a.rid].name,
                    'mapq': a.mapq,
                    'strand': "+-"[a.is_rev],
                    'NM': a.NM,
                    'cigar': parse_cigar(a.cigar, a.n_cigar), 'score': a.score})
                free(a.cigar)
            return result
        finally:
            free(opt)
            free(ar.a)

def open_bwa_index(bytes index_path):
    cdef int BWA_IDX_ALL = 0x7
    cdef bwaidx_t *idx = bwa_idx_load(index_path, BWA_IDX_ALL)
    if idx == NULL:
        raise IOError("Could not open " + index_path)
    cdef BwaIndex result = BwaIndex.__new__(BwaIndex)
    result.idx = idx
    result.path = index_path
    return result
