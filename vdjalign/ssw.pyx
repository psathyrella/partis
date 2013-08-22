# distutils: language = c++
# distutils: include_dirs = ssw-src/src
# distutils: sources = ssw-src/src/ssw_cpp.cpp ssw-src/src/ssw.c
import re
from libcpp.string cimport string
from libc.stdlib cimport free, malloc
from libc.string cimport strdup

cdef extern from "inttypes.h":
    ctypedef int int32_t
    ctypedef unsigned short uint16_t
    ctypedef unsigned char uint8_t

cdef extern from "ssw_cpp.h" namespace "StripedSmithWaterman":
    cdef cppclass Alignment:
        int32_t sw_score
        int32_t sw_score_next_best
        int32_t ref_begin
        int32_t ref_end
        int32_t ref_end_next_best
        int32_t query_begin
        int32_t query_end
        int32_t mismatches
        string cigar_string
    cdef cppclass Filter:
        pass
    cdef cppclass Aligner:
        Aligner()
        Aligner(uint8_t& match_score,
              uint8_t& mismatch_penalty,
              uint8_t& gap_opening_penalty,
              uint8_t& gap_extending_penalty)

        int Align(char *query, char *rf, int& ref_len,
                Filter& filter, Alignment* alignment) nogil

_cigar_re = re.compile(r'(\d+)([A-Z])')

cdef align_to_pair(Alignment a, bytes query, bytes ref):
    """
    Convert an alignment (with cigar string) to a tuple containing
    (aligned_ref, aligned_query)
    """
    cdef bytes cigar = a.cigar_string
    cdef bytes c, clen_str
    cdef int i, clen, ri, qi, j
    qi = a.query_begin
    ri = a.ref_begin
    q_aln = []
    r_aln = []

    cigar_parts = _cigar_re.findall(cigar)
    assert cigar_parts

    for clen_str, c in cigar_parts:
        clen = int(clen_str)

        if c == 'M':
            for j in range(clen):
                q_aln.append(query[qi + j])
                r_aln.append(ref[ri + j])
            ri += clen
            qi += clen
        elif c == 'I':
            for j in range(clen):
                q_aln.append(query[qi + j])
                r_aln.append('-')
            qi += clen
        elif c == 'D':
            for j in range(clen):
                q_aln.append('-')
                r_aln.append(ref[ri + j])
            ri += clen

    return ''.join(r_aln), ''.join(q_aln)

def align(bytes query, bytes ref, parse_cigar=True, uint8_t match_score=2,
        uint8_t mismatch_penalty=2, uint8_t gap_open_penalty=3, uint8_t gap_extend_penalty=1):
    """
    ``align(query, ref, parse_cigar=True)``

    Aligns ``query`` to ``ref`` using banded smith waterman.

    Returns a dict.
    """
    cdef Filter f
    cdef Aligner *aligner = new Aligner(match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty)
    cdef Alignment alignment
    cdef int res, ref_len = len(ref)
    cdef char *qry_str = strdup(query), *ref_str = strdup(ref)

    try:
        with nogil:
            res = aligner.Align(qry_str, ref_str, ref_len, f, &alignment)
    finally:
        del aligner
        free(qry_str)
        free(ref_str)
    if not res:
        raise ValueError("Error in alignment")

    result = {'sw_score': alignment.sw_score,
              'sw_score_next_best': alignment.sw_score_next_best,
              'ref_begin': alignment.ref_begin,
              'ref_end': alignment.ref_end,
              'ref_end_next_best': alignment.ref_end_next_best,
              'query_begin': alignment.query_begin,
              'query_end': alignment.query_end,
              'mismatches': alignment.mismatches,
              'cigar_string': alignment.cigar_string}
    if parse_cigar:
        r_aln, q_aln = align_to_pair(alignment, query, ref)
        result['ref_align'] = r_aln
        result['query_align'] = q_aln
    return result
