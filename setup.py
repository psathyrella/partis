#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.egg_info import egg_info
from setuptools.dist import Distribution


def check_system_dependencies():
    """Check that required system dependencies are available."""
    required_commands = ['scons', 'gcc', 'g++', 'python3', 'mafft']
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


def build_compiled_components():
    """Build the required C++ components using the build script."""

    print("Building partis compiled components (ig-sw and ham)...")

    # Check system dependencies first
    check_system_dependencies()

    base_dir = Path(__file__).parent.absolute()
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
    packages=['partis', 'partis.cache'],
    package_dir={'partis': 'partis'},
    python_requires='>=3.7',
    
    # Dependencies
    install_requires=[
        'biopython',
        'colored-traceback',
        'dendropy',
        'matplotlib',
        'numpy',
        'pandas',
        'psutil',
        'pysam',
        'PyYAML',
        'scikit-learn',
        'scipy',
        'seaborn',
        'six',
    ],
    
    # Build requirements
    setup_requires=[
        'scons',
        'setuptools-scm',
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
        ],
    },
    
    # Install Python scripts to PATH
    scripts=[
        'bin/cf-germlines.py',
        'bin/cf-alleles.py',
        'bin/extract-pairing-info.py',
        'bin/split-loci.py',
        'bin/get-naive-probabilities.py',
        'bin/compare-plotdirs.py',
        'bin/gctree-run.py',
        'bin/parse-output.py',
        'bin/partis-test.py',
    ],

    # Include data files
    include_package_data=True,
    package_data={
        'partis': [
            # Note: Using ../ to include files from outside the package directory
            # These are runtime dependencies that need to be bundled with the package
            '../bin/*',
            '../data/**/*',
            '../test/**/*',
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
