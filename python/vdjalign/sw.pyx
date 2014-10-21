# distutils: language = c
# distutils: include_dirs = src
# distutils: sources = src/ksw.c src/sw_align.c src/ig_align.c src/kstring.c
# distutils: extra_compile_args = -std=c99
# distutils: extra_link_args = -lz
from libc.stdlib cimport free, malloc
from libc.stdint cimport int32_t, uint8_t

import re


cdef extern from "sw_align.h":
    void align_reads(const char*,
                     const char*,
                     const char*,
                     const int32_t,
                     const int32_t,
                     const int32_t,
                     const int32_t,
                     const uint8_t,
                     const int32_t,
                     const int32_t,
                     const char*,
                     const char*) nogil

cdef extern from "ig_align.h":
    void ig_align_reads(const char*,
                        const uint8_t,
                        const char**,
                        const char*,
                        const char*,
                        const int32_t,
                        const int32_t,
                        const int32_t,
                        const int32_t,
                        const unsigned,
                        const int,
                        const unsigned,
                        const uint8_t,
                        const char*,
                        const char*) nogil

class InvalidReadGroupError(ValueError):
    pass

def ig_align(bytes ref_path,
             bytes qry_path,
             bytes output_path,
             list extra_ref_paths = [],
             int match=2,
             int mismatch=2,
             int gap_open=3,
             int gap_extend=1,
             int max_drop=0,
             int min_score=0,
             int bandwidth=150,
             int n_threads=1,
             bytes read_group=None):
    """
    :param ref_path: Path to reference sequence file (optionally gzipped) fasta
    :param qry_path: Path to query sequence file (optionally gzipped) fasta
    :param output_path: Path to write output (SAM)
    :param match: positive match score
    :param mismatch: positive mismatch penalty
    :param extra_ref_paths: Extra references to align the tail of each read against.
        If more than 1, should be in D, J order.
        Alignment occurs first against the J, then D.
    """
    cdef char* ref = ref_path
    cdef char* qry = qry_path
    cdef char* out = output_path
    cdef list read_group_id
    cdef char* rg = NULL
    cdef char* rg_id = NULL

    if max_drop < 0:
        raise ValueError("Invalid max drop: {0}".format(max_drop))

    if read_group is not None:
        read_group = read_group.replace('\\t', '\t')
        rg = read_group
        if not read_group.startswith('@RG'):
            raise InvalidReadGroupError('Does not start with @RG: "{0}"'.format(read_group))

        read_group_id = re.findall('\tID:([^\t]+)', read_group)
        if len(read_group_id) != 1:
            raise InvalidReadGroupError('Could not find ID: "{0}" {1}'.format(read_group, read_group_id))
        rg_id = read_group_id[0]

    cdef int32_t m = match, p = mismatch, go = gap_open, ge = gap_extend, ne = len(extra_ref_paths)
    cdef uint8_t threads = n_threads
    cdef char** extra = <char**>malloc(sizeof(char*) * len(extra_ref_paths))
    cdef unsigned md = max_drop, bw = bandwidth
    cdef int mscore = min_score

    cdef bytes r
    try:
        for i in xrange(len(extra_ref_paths)):
            r = extra_ref_paths[i]
            extra[i] = r
        with nogil:
            ig_align_reads(ref, ne, extra, qry, out, m, p, go, ge, md, mscore, bandwidth, threads, rg, rg_id)
    finally:
        free(extra)

def align(bytes ref_path,
          bytes qry_path,
          bytes output_path,
          int match=2,
          int mismatch=2,
          int gap_open=3,
          int gap_extend=1,
          int n_threads=1,
          int n_keep=-1,
          int max_drop=100,
          bytes read_group=None):
    """
    :param ref_path: Path to reference sequence file (optionally gzipped) fasta
    :param qry_path: Path to query sequence file (optionally gzipped) fasta
    :param output_path: Path to write output (SAM)
    :param match: positive match score
    :param mismatch: positive mismatch penalty
    :param n_keep: Maximum number of alignments to output
    :param max_drop: Maximum drop from best alignment score
    :param read_group: read group string
    """
    cdef char* ref = ref_path
    cdef char* qry = qry_path
    cdef char* out = output_path
    cdef list read_group_id
    cdef char* rg = NULL
    cdef char* rg_id = NULL

    if read_group is not None:
        read_group = read_group.replace('\\t', '\t')
        rg = read_group
        if not read_group.startswith('@RG'):
            raise InvalidReadGroupError('Does not start with @RG: "{0}"'.format(read_group))

        read_group_id = re.findall('\tID:([^\t]+)', read_group)
        if len(read_group_id) != 1:
            raise InvalidReadGroupError('Could not find ID: "{0}" {1}'.format(read_group, read_group_id))
        rg_id = read_group_id[0]

    cdef int32_t m = match, p = mismatch, go = gap_open, ge = gap_extend, md = max_drop
    cdef uint8_t threads = n_threads

    with nogil:
        align_reads(ref, qry, out, m, p, go, ge, threads, n_keep, md, rg, rg_id)

