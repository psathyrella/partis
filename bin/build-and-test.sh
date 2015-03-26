# NOTE this is for use by docker, and as such assumes you're root
set -eu

export LEIN_ROOT=1
cd /partis/packages/samtools/ && make
export PATH=$PWD:$PATH
cd /partis/packages/ighutil/ && make -C clj
pip install --user ./python
cd /partis/packages/ham/ && scons bcrham
cd /partis/
scons test
