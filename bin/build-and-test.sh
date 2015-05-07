# NOTE this is for use by docker, and as such assumes you're root
echo "running $0"
set -eu

echo "building samtools"
export LEIN_ROOT=1
cd /partis/packages/samtools/ && make
export PATH=$PWD:$PATH

echo "building ighutil"
cd /partis/packages/ighutil/ && make -C clj
pip install --user ./python

echo "building smctc"
cd /partis/packages/smctc/ && make

echo "building ham"
cd /partis/packages/ham/ && scons bcrham
cd /partis/

echo "testing"
scons test
