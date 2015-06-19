echo -e "\n--> running $0"
set -eu

if grep -v '/$' /proc/1/cgroup>/dev/null; then
    basedir=/partis  # if we're in docker
else
    basedir=$PWD
fi

echo -e "\n--> building samtools"
export LEIN_ROOT=1
cd $basedir/packages/samtools/ && make
export PATH=$PWD:$PATH

echo -e "\n--> building ighutil"
cd $basedir/packages/ighutil/ && make -C clj
pip install --user ./python

echo -e "\n--> building smctc"
cd $basedir/packages/smctc/ && make

echo -e "\n--> building ham"
cd $basedir/packages/ham/ && scons bcrham
cd $basedir/

echo -e "\n--> scons test"
scons test
