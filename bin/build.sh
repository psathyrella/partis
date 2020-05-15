#!/bin/bash

echo -e "\n--> running $0"
set -eu

if [ -f /.dockerenv ]; then  # if we're in docker TODO this file may go away in the future, so it'd probably be better to switch to the "grep docker /proc/1/cgroup" method
    basedir=/partis
else
    basedir=$PWD
fi

echo -e "\n--> building ig-sw"
cd $basedir/packages/ig-sw/src/ig_align/ && scons
cd $basedir/

echo -e "\n--> building ham"
cd $basedir/packages/ham/ && scons bcrham
cd $basedir/

if [ "$*" == "with-simulation" ]; then
    echo -e "\n--> building bpp-newlik (only used for simulation)"
    cd $basedir/packages/bpp-newlik/ && ./install.sh
    cd $basedir/
fi

echo -e "\n--> test"
./test/test.py --quick
