set -eu

export LEIN_ROOT=1
cd /partis/packages/samtools/ && make
mkdir /partis/_bin
ln -s $PWD/samtools /partis/_bin
cd /partis/packages/ighutil/ && make -C clj
pip install --user ./python
cd /partis/packages/ham/ && scons bcrham
cd /partis/
mkdir -p /true/plots
export PATH=/partis/_bin:$PATH
scons test
