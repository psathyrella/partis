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

Actually doing stuff, quickly
==============

Run viterbi to get the (by default five-) most likely annotations on one of the (simulated) seqs in test/simu.csv
```
./runpart.py --point_estimate --seqfile ./test/simu.csv --parameter_dir ./caches/parameters/hmm_parameters --n_max_queries 1 --debug 1
```
p.s. did I use 'annotations' correctly? I don\'t know!

And, you can run the forward pair algorithm to partition the simulated sequences:
```
./runpart.py --partition --pair --seqfile ./test/simu.csv --parameter_dir ./caches/parameters/hmm_parameters
```
This runs a hamming preclustering step, then a stripped-down (fast) HMM preclustering step, then a final HMM step.
The last few rows of the output are the clusters. Since the simulator ran with clones of 5 seqs, and the test file was made with a '$ head -n30 ...', there should
be five clusters of five, and one of four.

To see the bayes factors that this is clustering on, you can add '--debug 1', although that\'s a bit verbose, so you can also add '--no_clean', which tells it not to clean up its
temp files in /tmp/$USER/hmms/$pid. The pairwise scores are in files with 'pairscores' in the name.

You can also run viterbi pair, to see how pairs of sequences line up (it should be apparent which of these three are from the same event, and which are not)
```
./runpart.py --point_estimate --pair --seqfile ./test/simu.csv --parameter_dir ./caches/parameters/hmm_parameters --debug 1 --queries ' -4938873688098421278:8323195682841878031:40262812797573106'  # NOTE space before minus sign
```

Doing stuff that takes longer
==============

Now, the previous commands were reading info that I\'ve cached in ./caches, and simulated events in ./test/simu.csv, but you of course want to run from scratch.
I only cache it \'cause to get reliable parameter estimates you have to run over like ten thousand sequences, which takes a bit of time.

So now we go back and create the parameter files that we used in the previous steps, starting from some data in ./test/data.tsv

run on data to cache parameters and model files in parameter_dir
```
./runpart.py --cache_parameters --seqfile ./test/data.tsv --is_data --n_bases_skip 9 --v_right_length 56 --parameter_dir ./caches/new-parameters --ighutil_dir ~/.local/bin
```

this does the following:
  - runs smith-waterman on the data in test/data.tsv in order to estimate model parameters, which are written to parameter_dir/sw_parameters
  - read these parameters and use them to run the viterbi hmm on the same data file. This gives you better estimates of the parameters, which are written to parameter_dir/hmm_parameters
NOTE
  - ignore all the warnings about 'v right' (this tells you that a lot of the sequences have crappy best-v-matches, but you do not care a.t.m.)
  - should take ten or twenty minutes, maybe. the ham binary should be taking 100 percent cpu for most of that. We cannot really run this part with fewer
      seqs just as a test, because you need lots of seqs to get information on all the positions in all the gene versions
  - ignore all the warnings and errors about bad conserved codons and messed up cysteines and tryptophans. This is just telling you that there is unproductive rearrangements in the data.
    Well, and that my germline versions suck.
      
when it finishes, you can poke around in parameter_dir/hmm_parameters/ and see what is going into the model

Now you can run the simulator with these new parameters, too, to make six clusters/clones of five sequences each:
```
./runpart.py --simulate --parameter_dir ./caches/new-parameters/hmm_parameters --n_max_queries 30 --outdir ./caches/recombinator --debug 1
```
