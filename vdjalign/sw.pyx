# distutils: language = c
# distutils: include_dirs = src
# distutils: sources = src/ksw.c src/sw_align.c src/kstring.c
# distutils: extra_compile_args = -std=c99
# distutils: extra_link_args = -lz
from libc.stdlib cimport free, malloc
from libc.stdint cimport int32_t, uint8_t


cdef extern from "sw_align.h":
    void align_reads(const char*,
                     const char*,
                     const char*,
                     int32_t,
                     int32_t,
                     int32_t,
                     int32_t,
                     uint8_t) nogil

def align_reads_to_ref(bytes ref_path,
                       bytes qry_path,
                       bytes output_path,
                       int match=2,
                       int mismatch=2,
                       int gap_open=3,
                       int gap_extend=1):
    cdef char* ref = ref_path
    cdef char* qry = qry_path
    cdef char* out = output_path

    cdef int32_t m = match, p = mismatch, go = gap_open, ge = gap_extend
    cdef uint8_t threads = 12

    with nogil:
        align_reads(ref, qry, out, m, p, go, ge, threads)

