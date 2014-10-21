from setuptools import setup, find_packages

from Cython.Build import cythonize
import versioneer

versioneer.versionfile_source = 'vdjalign/_version.py'
versioneer.versionfile_build = 'vdjalign/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = 'vdjalign-' # dirname like 'myproject-1.2.0'

ext_modules = cythonize(['vdjalign/*.pyx'])

for i in ext_modules:
    i.define_macros = [('VDJALIGN_VERSION', versioneer.get_version())]

setup(
    name='vdjalign',
    #version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={'vdjalign': ['imgt/data/*.fasta']},
    ext_modules=ext_modules,
    version='0.1.0',
    install_requires=['pysam', 'networkx'],
    entry_points={
        'console_scripts': [
            'vdjalign = vdjalign.scripts.vdjalign:main',
    ]},
    test_suite='vdjalign.test',)
