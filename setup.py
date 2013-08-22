from setuptools import setup, find_packages

from Cython.Build import cythonize

ext_modules = cythonize(['vdjalign/*.pyx'])

setup(
    name='vdjalign',
    packages=find_packages(),
    package_data={'vdjalign': ['imgt/data/*.fasta']},
    ext_modules=ext_modules,
    install_requires=['pysam'],
    entry_points={
        'console_scripts': [
            'vdjalign = vdjalign.scripts.vdjalign:main',
    ]},
    test_suite='vdjalign.test',
)
