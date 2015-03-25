set -eu

export LEIN_ROOT=1
echo 'HERE'
cd /partis/packages/samtools/ && make
mkdir -p /partis/_bin
ln -s $PWD/samtools /partis/_bin
cd /partis/packages/ighutil/ && make -C clj
pip install --user ./python
cd /partis/packages/ham/ && scons bcrham
cd /partis/
mkdir -p /plots /v-per-base/plots/ /d-per-base/plots/ /j-per-base/plots/ /partis/_output/test/plots/params/data
export PATH=/partis/_bin:$PATH
scons test
