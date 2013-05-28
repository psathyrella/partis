import os.path
from setuptools import setup, find_packages

from Cython.Build import cythonize

compile_args = ['-g', '-Wall', '-O0', '-DUSE_MALLOC_WRAPPERS', '-DHAVE_PTHREAD']

bwa_libs = ['z', 'm']
bwa_c_src = ['utils', 'kstring', 'ksw', 'bwt', 'bntseq', 'bwa', 'bwamem', 'bwamem_pair', 'malloc_wrap']
bwa_c_src = [os.path.join('external', 'bwa', i + '.c') for i in bwa_c_src]

ext_modules = cythonize(['vdjalign/*.pyx'],
        aliases={'BWA_SRC': bwa_c_src, 'BWA_FLAGS': compile_args, 'BWA_LIBS': bwa_libs,
                 'BWA_INCLUDE': ['external/bwa']})

setup(
    name='vdjalign',
    packages=find_packages(),
    package_data={'vdjalign': ['imgt/data/*.fasta']},
    ext_modules=ext_modules,
)
