# distutils: language = c++
# distutils: include_dirs = src
# distutils: extra_compile_args = -std=c++0x -Wall -Wextra
# distutils: sources = src/c_codons.cpp
"""
Codon translation
"""
from libcpp.string cimport string
from cpython cimport bool

cdef extern from "c_codons.hpp" namespace "codons":
    string translate_dna(string&) except +
    void toupper(string&) except +

def _toupper(bytes s not None):
    cdef string st = s
    toupper(st)
    return st

def translate(bytes s not None, bool is_upper=False):
    """
    translate(string s) -> string
    Translate string ``s`` to amino acids.

    Unknown codons are marked 'X'
    """
    cdef string st = s
    if not is_upper:
        toupper(st)
    return translate_dna(st)
