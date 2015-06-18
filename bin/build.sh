echo "running $0"
set -eu

if grep -v '/$' /proc/1/cgroup>/dev/null; then
    basedir=/partis  # if we're in docker
else
    basedir=$PWD
fi

echo "building samtools"
export LEIN_ROOT=1
cd $basedir/packages/samtools/ && make
export PATH=$PWD:$PATH

echo "building ighutil"
cd $basedir/packages/ighutil/ && make -C clj
pip install --user ./python

echo "building smctc"
cd $basedir/packages/smctc/ && make

echo "building ham"
cd $basedir/packages/ham/ && scons bcrham
cd $basedir/

echo "testing"
scons test
