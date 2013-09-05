# distutils: language = c
# distutils: include_dirs = src
# distutils: sources = src/ksw.c src/sw_align.c src/kstring.c
# distutils: extra_compile_args = -std=c99
# distutils: extra_link_args = -lz
from libc.stdlib cimport free, malloc
from libc.stdint cimport int32_t, uint8_t

import re


cdef extern from "sw_align.h":
    void align_reads(const char*,
                     const char*,
                     const char*,
                     int32_t,
                     int32_t,
                     int32_t,
                     int32_t,
                     uint8_t,
                     const int32_t,
                     const char*,
                     const char*) nogil

class InvalidReadGroupError(ValueError):
    pass


def align(bytes ref_path,
          bytes qry_path,
          bytes output_path,
          int match=2,
          int mismatch=2,
          int gap_open=3,
          int gap_extend=1,
          int n_threads=1,
          int n_keep=-1,
          bytes read_group=None):
    """
    :param ref_path: Path to reference sequence file (optionally gzipped) fasta
    :param qry_path: Path to query sequence file (optionally gzipped) fasta
    :param output_path: Path to write output (SAM)
    :param match: positive match score
    :param mismatch: positive mismatch penalty
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

    cdef int32_t m = match, p = mismatch, go = gap_open, ge = gap_extend
    cdef uint8_t threads = n_threads

    with nogil:
        align_reads(ref, qry, out, m, p, go, ge, threads, n_keep, rg, rg_id)

