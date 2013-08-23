# distutils: language = c
# distutils: include_dirs = lib/ssw/src src
# distutils: sources = lib/ssw/src/ssw.c src/ssw_align.c
# distutils: extra_compile_args = -fopenmp
# distutils: extra_link_args = -fopenmp -lz
from libc.stdlib cimport free, malloc
from libc.stdint cimport int32_t


cdef extern from "ssw_align.h":
    void align_reads(const char*,
                     const char*,
                     const char*,
                     int32_t,
                     int32_t,
                     int32_t,
                     int32_t) nogil

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

    with nogil:
        align_reads(ref, qry, out, m, p, go, ge)

