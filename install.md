Note that if you only want to run things (not edit the source), it is probably easier to just use the [docker image](https://registry.hub.docker.com/u/psathyrella/partis/).

dependecies
==============
required
--------------
  - scons
  - pip
  - java 7 (for ighutil)
  - ROOT
  - python packages
    - pysam
    - pyyaml
    - cython
    - networkx
    - decorator
    - lxml
    - bs4 (beautiful soup)
    - pandas
    - seaborn

required only for simulator
--------------
  - dendropy
  - TreeSim

included
--------------
  - tclap
  - yaml-cpp
  - ham
  - samtools
  - bppseqgen

installation
==============

```
git clone git@github.com:psathyrella/partis
cd partis/packages/ham/
scons
cd ../samtools/
make
ln -s $PWD/samtools ~/bin/
export PATH=~/bin:$PATH
cd ../ighutil/
make -C clj
pip install --user ./python
cd ..
```
