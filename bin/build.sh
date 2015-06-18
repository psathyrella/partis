echo "\n--> running $0"
set -eu

if grep -v '/$' /proc/1/cgroup>/dev/null; then
    basedir=/partis  # if we're in docker
else
    basedir=$PWD
fi

echo "\n--> building samtools"
export LEIN_ROOT=1
cd $basedir/packages/samtools/ && make
export PATH=$PWD:$PATH

echo "\n--> building ighutil"
cd $basedir/packages/ighutil/ && make -C clj
pip install --user ./python

echo "\n--> building smctc"
cd $basedir/packages/smctc/ && make

echo "\n--> building ham"
cd $basedir/packages/ham/ && scons bcrham
cd $basedir/

echo "\n--> scons test"
scons test
