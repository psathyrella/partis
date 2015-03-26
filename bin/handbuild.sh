cd packages/ham/
scons bcrham
cd ../samtools/
export LEIN_ROOT=1
make
export PATH=$PWD:$PATH
cd ../ighutil/
make -C clj
pip install --user ./python
cd ../../
