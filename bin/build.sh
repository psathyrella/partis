echo -e "\n--> running $0"
set -eu

if [ -f /.dockerinit ]; then  # if we're in docker
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

echo -e "\n--> test"
./test/test.py --quick
