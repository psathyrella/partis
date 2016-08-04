# distutils: language = c++
# distutils: include_dirs = src
# distutils: sources = src/deduplicate_fastq.cpp
# distutils: extra_compile_args = -std=c++0x
# distutils: extra_link_args = -lz
from libc.stdlib cimport free, malloc
from libc.stdint cimport int32_t, uint8_t

cdef extern from "deduplicate_fastq.hpp":
    void deduplicate_fastq(const char*, const char*, const int) nogil

def deduplicate(bytes input_path, bytes output_path, int est_size=2097152): # 2**21
    cdef const char* inp = input_path
    cdef const char* out = output_path
    cdef int est = est_size

    with nogil:
        deduplicate_fastq(inp, out, est)
