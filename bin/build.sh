#!/bin/bash

echo -e "\n--> running $0"
set -eu

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
basedir=`dirname $SCRIPT_DIR`  # go up one level

echo -e "\n--> building ig-sw"
cd $basedir/packages/ig-sw/src/ig_align/ && scons
cd $basedir/

echo -e "\n--> building ham"
cd $basedir/packages/ham/ && scons bcrham
cd $basedir/

if [[ "$OSTYPE" == "darwin"* ]]; then  # FastTree only has linux+windows binaries
    echo -e "\n--> building FastTree"
    cd $basedir/packages/FastTree/ && gcc -O3 -fopenmp-simd -funsafe-math-optimizations -march=native -o FastTree FastTree.c -lm
    cd $basedir/
    ln -sfv $basedir/packages/FastTree/FastTree bin/FastTree-macos
fi

if [ "$*" == "with-simulation" ]; then
    echo -e "\n--> building bpp-newlik (only used for simulation)"
    cd $basedir/packages/bpp-newlik/ && ./install.sh  # the bpp-phyl step is really really incredibly slow
    cd $basedir/
fi
