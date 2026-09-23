#!/usr/bin/env python3

import os
import sys
import platform
import subprocess
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.egg_info import egg_info
from setuptools.dist import Distribution


def check_system_dependencies(required_commands):
    """Check that the given system commands are available."""
    missing = []

    for cmd in required_commands:
        try:
            subprocess.run([cmd, '--version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            missing.append(cmd)

    # Check python3-venv separately since it's a module, not a command
    try:
        subprocess.run([sys.executable, '-m', 'venv', '--help'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        missing.append('python3-venv')

    if missing:
        raise RuntimeError('Missing required system dependencies: %s. Please install all packages listed in the manual (docs/install.md#installation-with-pip), then try the installation again.' % ', '.join(missing))


def build_fasttree_macos(base_dir):
    """Compile FastTree on macOS (it has no prebuilt mac binary; Linux uses the bundled
    one). Backend-independent -- bin/build.sh does the same for the C++ path."""
    if sys.platform != 'darwin':
        return
    ftdir = base_dir / 'packages' / 'FastTree'
    print('Building FastTree (macOS)...')
    subprocess.run(['gcc', '-O3', '-fopenmp-simd', '-funsafe-math-optimizations', '-march=native',
                    '-o', 'FastTree', 'FastTree.c', '-lm'], cwd=str(ftdir), check=True)
    link = base_dir / 'bin' / 'FastTree-macos'
    if link.is_symlink() or link.exists():
        link.unlink()
    os.symlink(ftdir / 'FastTree', link)


def build_compiled_components():
    """Build the backend requested by $PARTIS_BACKEND.

    Default is 'cpp' everywhere except Apple Silicon (macOS arm64), where it's 'zig':
    the C ig-sw uses SSE2/emmintrin and can neither compile nor run on arm64 (issue
    #330), so zig is the only backend that works there.

    'cpp': compile the C++ bcrham + C ig-sw via bin/build.sh (needs the C/C++
      toolchain: scons, gcc, g++).
    'zig': build the in-tree Zig backend via bin/zig-build.sh (fetches its own Zig
      compiler; needs only curl+tar, no C/C++ toolchain / gsl / yaml-cpp). On non-arm64
      platforms this does NOT change the default backend: the C++ backend remains the
      default and zig is selected at run time with --zig. (In a zig-only install the
      C++ backend isn't built, so you must pass --zig -- except on Apple Silicon, where
      bin/partis defaults to zig automatically since the C++ backend can't run.)
    """
    base_dir = Path(__file__).parent.absolute()
    backend = os.environ.get('PARTIS_BACKEND')
    if backend is None:  # Apple Silicon can't build the SSE2 C++ ig-sw (issue #330), so default to zig there
        backend = 'zig' if (sys.platform == 'darwin' and platform.machine() == 'arm64') else 'cpp'
    backend = backend.lower()
    if backend not in ('cpp', 'zig'):
        raise RuntimeError("PARTIS_BACKEND must be 'cpp' or 'zig' (got '%s')" % backend)

    if backend == 'zig':
        print("Building partis zig backend (PARTIS_BACKEND=zig)...")
        if not any((base_dir / 'packages/zig-core').iterdir()):
            raise RuntimeError('Submodule packages/zig-core is empty. Run: git submodule update --init packages/zig-core')
        check_system_dependencies(['curl', 'tar', 'mafft'])  # zig-build.sh fetches its own compiler; no gcc/scons/gsl/yaml-cpp
        result = subprocess.run([str(base_dir / 'bin' / 'zig-build.sh')], cwd=str(base_dir))
        if result.returncode != 0:
            raise Exception('zig-build.sh failed with exit code %d' % result.returncode)
        build_fasttree_macos(base_dir)  # FastTree is backend-independent; no prebuilt mac binary, so compile it (bin/build.sh does this for the cpp path)
        print("✓ Successfully built zig backend (run partis with --zig to use it)")
        return

    print("Building partis compiled components (ig-sw and ham)...")
    for submod in ['packages/ham', 'packages/ig-sw']:
        if not any((base_dir / submod).iterdir()):
            raise RuntimeError('Submodule %s is empty. Run: git submodule update --init %s' % (submod, submod))

    check_system_dependencies(['scons', 'gcc', 'g++', 'python3', 'mafft'])
    build_script = base_dir / "bin" / "build.sh"

    if not build_script.exists():
        raise Exception(f"Build script not found at {build_script}")

    print("Compiling C++ components...")
    result = subprocess.run([str(build_script)], cwd=str(base_dir))

    if result.returncode != 0:
        raise Exception(f"Build script failed with exit code {result.returncode}")

    print("✓ Successfully built ig-sw and ham binaries")


def make_custom_command(base_class):
    """Factory function to create custom command classes that build C++ components."""
    class CustomCommand(base_class):
        def run(self):
            build_compiled_components()
            base_class.run(self)
    return CustomCommand


CustomBuildPy = make_custom_command(build_py)
CustomDevelop = make_custom_command(develop)
CustomInstall = make_custom_command(install)
CustomEggInfo = make_custom_command(egg_info)


class BinaryDistribution(Distribution):
    """Distribution which always forces a binary package with platform name"""
    def has_ext_modules(self):
        return True


# Read the README file
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='partis-bcr',
    use_scm_version=True,
    description='B- and T-cell receptor sequence annotation, simulation, clonal family and germline inference',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/psathyrella/partis',
    author='Duncan Ralph',
    author_email='dkralph@gmail.com',
    license='GPL-3.0-or-later',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
    ],
    keywords='immunology bioinformatics bcr tcr antibody sequence-analysis',
    
    # Packages and package data
    packages=['partis', 'partis.cache', 'partis.scripts'],
    package_dir={'partis': 'partis'},
    python_requires='>=3.7',
    
    # Dependencies
    install_requires=[
        'biopython',
        'colored-traceback',
        'dendropy',
        'fire',
        'matplotlib',
        'numpy',
        'pandas',
        'psutil',
        'pot',
        'pysam',
        'PyYAML',
        'scikit-learn',
        'scipy',
        'seaborn',
        'six',
    ],
    
    # Optional dependencies
    extras_require={
        'plotting': [
            'circlify',
            'ete3',
            'PyQt5',
            'joypy',
            'levenshtein',
        ],
    },
    
    # Entry points
    entry_points={
        'console_scripts': [
            'partis=partis.main:main',
            'cf-germlines.py=partis.scripts.cf_germlines:main',
            'cf-alleles.py=partis.scripts.cf_alleles:main',
            'extract-pairing-info.py=partis.scripts.extract_pairing_info:main',
            'split-loci.py=partis.scripts.split_loci:main',
            'get-naive-probabilities.py=partis.scripts.get_naive_probabilities:main',
            'compare-plotdirs.py=partis.scripts.compare_plotdirs:main',
            'gctree-run.py=partis.scripts.gctree_run:main',
            'make-html=partis.scripts.make_html:main',
            'parse-output.py=partis.scripts.parse_output:main',
            'partis-test.py=partis.scripts.partis_test:main',
        ],
    },

    # Include data files
    include_package_data=True,
    package_data={
        'partis': [
            # Note: Using ../ to include files from outside the package directory
            # These are runtime dependencies that need to be bundled with the package
            '../bin/*',
            '../data/**/*',
            '../test/*.py',
            '../test/*.fa',
            '../test/*.yaml',
            '../test/*.yml',
            '../test/*.sh',
            '../test/*.nwk',
            '../test/paired-data/**/*',
            '../test/ref-results/**/*',
            '../test/paired/ref-results/**/*',
            '../packages/ham/bcrham',
            '../packages/ig-sw/src/ig_align/ig-sw',
            '../packages/bpp/bin/bppseqgen',
            '../packages/bpp/lib/*',
        ],
    },
    
    # Custom build commands
    cmdclass={
        'build_py': CustomBuildPy,
        'develop': CustomDevelop,
        'install': CustomInstall,
        'egg_info': CustomEggInfo,
    },
    
    # Force binary distribution
    distclass=BinaryDistribution,
    
    # URLs
    project_urls={
        'Bug Reports': 'https://github.com/psathyrella/partis/issues',
        'Source': 'https://github.com/psathyrella/partis',
        'Documentation': 'https://github.com/psathyrella/partis/tree/main/docs',
    },
)
