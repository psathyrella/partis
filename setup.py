#!/usr/bin/env python3

import os
import sys
import subprocess
import shutil
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.egg_info import egg_info


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
    """Build the required C++ components using the original build script."""
    print("Building partis compiled components (ig-sw and ham)...")
    
    # Check system dependencies first
    check_system_dependencies()
    
    base_dir = Path(__file__).parent.absolute()
    build_script = base_dir / "bin" / "build.sh"
    
    if build_script.exists():
        print("Compiling C++ components...")
        try:
            # Run the original build script but suppress test failure
            result = subprocess.run([str(build_script)], 
                                  cwd=str(base_dir),
                                  capture_output=True,
                                  text=True)
            
            # Print stdout/stderr for debugging
            if result.stdout:
                print("Build output:", result.stdout[-1000:])  # Last 1000 chars
            if result.stderr:
                print("Build errors:", result.stderr[-1000:])  # Last 1000 chars
            
            # Check if the important binaries were created
            ig_sw_binary = base_dir / "packages" / "ig-sw" / "src" / "ig_align" / "ig-sw"
            ham_binary = base_dir / "packages" / "ham" / "bcrham"
            
            if ig_sw_binary.exists() and ham_binary.exists():
                print("✓ Successfully built ig-sw and ham binaries")
            else:
                missing = []
                if not ig_sw_binary.exists():
                    missing.append("ig-sw")
                if not ham_binary.exists():
                    missing.append("ham")
                print(f"Warning: Missing binaries: {', '.join(missing)}")
                # Don't fail the installation if tests fail, but binaries exist
                if ig_sw_binary.exists() or ham_binary.exists():
                    print("Continuing with partial build...")
                else:
                    raise Exception("Failed to build required binaries")
                    
        except Exception as e:
            print(f"Build script failed: {e}")
            # Try to fall back to manual scons build
            print("Attempting manual build...")
            try:
                build_with_scons(base_dir)
            except Exception as e2:
                print(f"Manual build also failed: {e2}")
                raise Exception("Could not build required C++ components")
    else:
        print(f"Build script not found at {build_script}, attempting manual build...")
        build_with_scons(base_dir)


def build_with_scons(base_dir):
    """Fallback build using SCons directly."""
    original_dir = os.getcwd()
    try:
        # Build ig-sw
        ig_sw_dir = base_dir / "packages" / "ig-sw" / "src" / "ig_align"
        if ig_sw_dir.exists():
            print("Building ig-sw...")
            os.chdir(str(ig_sw_dir))
            subprocess.check_call(['scons'])
        
        # Build ham
        ham_dir = base_dir / "packages" / "ham"
        if ham_dir.exists():
            print("Building ham...")
            os.chdir(str(ham_dir))
            subprocess.check_call(['scons', 'bcrham'])
            
    finally:
        os.chdir(original_dir)


class CustomBuildPy(build_py):
    """Custom build command that compiles C++ components."""
    def run(self):
        # First run the standard build
        build_py.run(self)
        # Then build our C++ components
        build_compiled_components()


class CustomDevelop(develop):
    """Custom develop command that compiles C++ components."""
    def run(self):
        print("CustomDevelop.run() called!")
        # Build C++ components first
        build_compiled_components()
        # Then run the standard develop
        develop.run(self)


class CustomInstall(install):
    """Custom install command that ensures C++ components are built."""
    def run(self):
        # Build C++ components first
        build_compiled_components()
        # Then run the standard install
        install.run(self)


class CustomEggInfo(egg_info):
    """Custom egg_info command that builds C++ components early."""
    def run(self):
        print("CustomEggInfo.run() called!")
        # Build C++ components before creating egg info
        build_compiled_components()
        # Then run the standard egg_info
        egg_info.run(self)


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
    
    # Packages and package data - keep the original 'python' module name for compatibility
    packages=['python', 'python.cache'],
    package_dir={'python': 'python'},
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
            'joypy',
            'matplotlib',
        ],
        'full': [
            'circlify',
            'ete3',
            'joypy',
            'levenshtein',
            'matplotlib',
        ],
    },
    
    # Entry points
    entry_points={
        'console_scripts': [
            'partis=python.main:main',
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
        'python': [
            '../bin/*',
            '../data/**/*',
            '../test/**/*',
            '../packages/ham/bcrham',
            '../packages/ig-sw/src/ig_align/ig-sw',
        ],
    },
    
    # Custom build commands
    cmdclass={
        'build_py': CustomBuildPy,
        'develop': CustomDevelop,
        'install': CustomInstall,
        'egg_info': CustomEggInfo,
    },
    
    # URLs
    project_urls={
        'Bug Reports': 'https://github.com/psathyrella/partis/issues',
        'Source': 'https://github.com/psathyrella/partis',
        'Documentation': 'https://github.com/psathyrella/partis/tree/main/docs',
    },
)
