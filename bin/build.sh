echo -e "\n--> running $0"
set -eu

if [ -f /.dockerinit ]; then  # if we're in docker
    basedir=/partis
else
    basedir=$PWD
fi

echo -e "\n--> building samtools"
export LEIN_ROOT=1
cd $basedir/packages/samtools/ && make

echo -e "\n--> building ighutil"
cd $basedir/packages/ighutil/ && make -C clj
pip install --user ./python

echo -e "\n--> building ham"
cd $basedir/packages/ham/ && scons bcrham
cd $basedir/

echo -e "\n--> test"
./test/test.py --quick
