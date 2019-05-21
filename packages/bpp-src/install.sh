bd=$PWD #/fh/fast/matsen_e/dralph/code/partis-dev/packages/bpp-src

cd $bd

packs="eigen3 bpp-core bpp-seq bpp-phyl bpp-popgen bppsuite"
for pack in $packs; do 
    export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$bd/$pack/_build
done

for pack in $packs; do
    echo -e "\n$pack"

    # # clone repo
    # if [ $pack == "eigen3" ]; then
    # 	repo=eigenteam/eigen-git-mirror
    # else
    # 	repo=BioPP/$pack
    # fi
    # git clone git@github.com:$repo $pack

    cd $pack

    # # switch branches
    # if [ $pack != "eigen3" ]; then  # crashes on one of the bio repos, which doesn't have a newlik branch
    # 	branch=newlik # devel
    # 	git branch --track $branch origin/$branch
    # 	git checkout $branch
    # fi

    # # configure/compile
    # mkdir -p _build
    # cd _build
    # cmake -Wno-dev -DCMAKE_INSTALL_PREFIX=$bd/_build ..  # -DCMAKE_PREFIX_PATH=/home/dralph/work/partis/packages/bpp-src/bpp-core
    # make install
    # cd ..

    # git status
    git diff >../$pack.patch

    cd ..
done
