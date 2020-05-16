#!/bin/bash

# install the newlik branch [update: newlik has been merged to master] of bpp (another, older version is also pre-installed in partis/packages/bpp/, and unfortunately we do really need both, at least for now), which you need to run simulation with --per-base-mutation set
action="compile"  # choices: "compile", "clone", "pull" NOTE clone and pull need updating
bppdir=$PWD  # NOTE must be run from within partis/packages/bpp-newlik (well, otherwise partis won't find it afterwards)

packs="eigen3 bpp-core bpp-seq bpp-phyl bpp-popgen bppsuite"
for pack in $packs; do 
    export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$bppdir/_build/$pack
done

for pack in $packs; do
    echo -e "\n$pack"

    if [ $pack == "eigen3" ]; then
    	repo=eigenteam/eigen-git-mirror
    else
    	repo=psathyrella/$pack  # BioPP
    fi

    if [ "$action" == "clone" ]; then
	git clone git@github.com:$repo $pack
    fi

    if ! [ -d $pack ]; then  # don't want to keep going if the dir isn't there
	echo "dir $pack doesn't exist in $PWD (run this script from packages/bpp-newlik/)"
	exit 1
    fi

    # i think i don't need any of this any more, as of now (may 2020) they've merged newlik into master pretty recently, so i think i can just use the master branch
    # i'll need to change this if i want to update versions, of course, but that shouldn't be for a while
    cd $pack
    if [ $pack == "eigen3" ] || [ $pack == "bpp-popgen" ]; then  # crashes on one of the bio repos, which doesn't have a newlik branch
    	branch=master
    else
	branch=newlik
    fi
    if [ "$action" == "clone" ]; then
	echo "need to switch to branch name with tmp- in front for new forks in psathyrella"
	exit 1
    	git branch --track $branch origin/$branch
    	git checkout $branch
    fi
    if [ "$action" == "pull" ]; then
	echo "need to switch to branch name with tmp- in front for new forks in psathyrella"
	exit 1
    	git checkout $branch
	git pull origin $branch
    fi
    cd ..

    # this was necessary for the ~may 2019 version of the newlik branch, which was missing a couple includes
    # if [ -f ../$pack.patch ]; then
    # 	git apply --verbose ../$pack.patch
    # fi

    # cd $pack
    # # git checkout -b tmp-$branch
    # # git push origin tmp-$branch
    # cd ..
    # continue

    if [ "$action" == "compile" ]; then
	mkdir -p _build/$pack
	cd _build/$pack
	cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=$bppdir/_build $bppdir/$pack/
	make install
	if ! [ $? -eq 0 ]; then  # don't want to keep going if one compile failed
	    echo "$pack installation failed"
	    exit 1
	fi
	cd ../..
    fi
done
