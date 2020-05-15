#!/bin/bash

if [ -f /.dockerenv ]; then  # if we're in docker TODO this file may go away in the future, so it'd probably be better to switch to the "grep docker /proc/1/cgroup" method
    basedir=/partis
else
    basedir=$PWD
    echo "using PWD as base partis dir: $PWD"
fi

# install the newlik branch of bpp (another, older version is also pre-installed in partis/packages/bpp/, and unfortunately we do really need both, at least for now), which you need to run simulation with --per-base-mutation set
bppdir=$basedir/packages/bpp-newlik
clone="false"  # if this is false, we assume the source for each subpackage is already here because git submodule was already run; if this is true, we clone each subpackage and then checkout the proper branch (which shouldn't need to happen any more unless i want to update versions or something)

cd $bppdir

packs="eigen3 bpp-core bpp-seq bpp-phyl bpp-popgen bppsuite"
for pack in $packs; do 
    export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$bppdir/$pack/_build
done

for pack in $packs; do
    echo -e "\n$pack"

    if [ $pack == "eigen3" ]; then
    	repo=eigenteam/eigen-git-mirror
    else
    	repo=BioPP/$pack
    fi

    if [ "$clone" == "true" ]; then
	git clone git@github.com:$repo $pack
    fi

    cd $pack
    if ! [ $? -eq 0 ]; then  # don't want to keep going if the dir isn't there
	echo "cd $pack failed"
	exit 1
    fi

    if [ "$clone" == "true" ]; then
	# switch branches
	if [ $pack != "eigen3" ]; then  # crashes on one of the bio repos, which doesn't have a newlik branch
    	    branch=newlik # devel
    	    git branch --track $branch origin/$branch
    	    git checkout $branch
	fi
    fi

    if [ -f ../$pack.patch ]; then  # the version i checked out was missing a couple includes
	git apply --verbose ../$pack.patch
    fi

    # configure/compile
    mkdir -p _build
    cd _build
    echo $PWD
    cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=$bppdir/_build ..  # -DCMAKE_PREFIX_PATH=/home/dralph/work/partis/packages/bpp-src/bpp-core
    make install
    if ! [ $? -eq 0 ]; then  # don't want to keep going if one compile failed
	echo "$pack installation failed"
	exit 1
    fi
    cd ..

    # git status
    # git diff >../$pack.patch

    cd ..
done
