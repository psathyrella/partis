dependecies
==============
required
--------------
  - bppseqgen (*very* recent version necessary in order to allow per-residue mutation frequency specification)
  - scons
  - pip
  - pysam
  - pyyaml
  - cython
  - networkx
  - decorator
  - TreeSim (R package, needed to generate trees)

optional
--------------
  - ROOT

included in ham/
--------------
  - tclap   
  - yaml-cpp

installation
==============

samtools
--------------

you need a very recent version. Download, say, 1.1: http://sourceforge.net/projects/samtools/files/samtools/1.1/
bunzip2 and untar it, then:
```
cd samtools-1.1
make
cd ..
ln -s $PWD/samtools-1.1/samtools ~/bin/  # get the samtools binary in your path. I do it like this
```

ighutil
--------------
```
git clone git@github.com:cmccoy/ighutil.git
cd ighutil
make -C clj  # NOTE requires java 7
pip install --user ./python  # NOTE without the './' this will do something very different
cd ..
```

ham
--------------
```
git clone git@github.com:psathyrella/ham
cd ham/
scons
cd ..
```

partis
--------------
```
git clone git@github.com:psathyrella/partis
cd partis/
```

run on data to cache parameters and model files in <parameter_dir>
--------------
./runpart.py --cache_parameters --seqfile ./test/data.tsv --is_data --n_bases_skip 9 --v_right_length 56 --parameter_dir ./caches/parameters --ham_dir ../ham --ighutil_dir ~/.local/bin

this does the following:
  - runs smith-waterman on the data in test/data.tsv in order to estimate model parameters, which are written to tmp/parameters/sw_parameters
  - read these parameters and use them to run the viterbi hmm on the same data file. This gives you better estimates of the parameters, which are written to tmp/parameters/hmm_parameters
NOTE
  - ignore all the warnings about 'v right' (this tells you that a lot of the sequences have crappy best-v-matches, but you do not care a.t.m.)
  - should take ten or twenty minutes, maybe. the ham binary should be taking 100 percent cpu for most of that. We cannot really run this part with fewer
      seqs just as a test, because you need lots of seqs to get information on all the positions in all the gene versions
  - ignore all the warnings and errors about bad conserved codons and messed up cysteines and tryptophans. This is just telling you that there is unproductive rearrangements in the data.
      
when it finishes, you can poke around in caches/parameters/hmm_parameters/ and see what is going into the model


run the simulator (using the parameter files in <parameter_dir>) to make six clusters/clones of five sequences each:
--------------
```
./runpart.py --simulate --parameter_dir ./caches/parameters/hmm_parameters --n_max_queries 30 --outdir ./caches/recombinator --debug 0 --ham_dir ../ham
```

Run viterbi on one of the (simulated) seqs in 
--------------
```
./runpart.py --point_estimate --seqfile tmp/recombinator/simu.csv --parameter_dir tmp/parameters/hmm_parameters --n_max_queries 1 --debug 1 --ham_dir ../ham
```

and finally you can run the forward pair algorithm to partition the sequences simulated sequences:
```
./runpart.py --partition --pair --seqfile tmp/recombinator/simu.csv --parameter_dir tmp/parameters/hmm_parameters --ham_dir ../ham
```
The last few rows of the output are the clusters. Since the simulator ran with clones of 5 seqs, and the test file was made with a '$ head -n30', there should
be a five clusters of five, and one of four.

You can also run viterbi pair, to see how pairs of sequences line up (it should be apparent which of these three are from the same event, and which are not)
NOTE I run this on test/simu.csv, which was *not* just made, because otherwise I do not know which query hashes to add to the --queries option
```
./runpart.py --point_estimate --pair --seqfile test/simu.csv --parameter_dir tmp/parameters/hmm_parameters --debug 1 --ham_dir ../ham --queries ' -4938873688098421278:8323195682841878031:40262812797573106'  # NOTE space before minus sign
```
