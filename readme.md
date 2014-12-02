dependecies
==============
required
--------------
  - scons
  - pip
  - java 7 (for ighutil)
  - python packages
    - pysam
    - pyyaml
    - cython
    - networkx
    - decorator
    - lxml
    - bs4 (beautiful soup)

required only for simulator
--------------
  - dendropy
  - TreeSim

optional
--------------
  - ROOT

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
cd partis/ham/
scons
cd ../samtools/
make
ln -s $PWD/samtools ~/bin/
cd ../ighutil/
make -C clj
pip install --user ./python
cd ..
```

Actually doing stuff, quickly
==============

Run viterbi to get the (by default five-) most likely annotations on one of the (simulated) seqs in data/test/simu.csv
```
./runpart.py --point_estimate --seqfile data/test/simu.csv --parameter_dir caches/parameters/hmm_parameters --n_max_queries 1 --debug 1
```

And, you can run the forward pair algorithm to partition the simulated sequences (adjust <n_procs> as appropriate):
```
./runpart.py --partition --pair --seqfile data/test/simu.csv --parameter_dir caches/parameters/hmm_parameters --n_procs 10 --n_max_per_region 2
```
This runs a hamming preclustering step, then a stripped-down (fast) HMM preclustering step, then a final HMM step. The last few rows of the output are the clusters.
Since the simulator ran with clones of 5 seqs, there should be six clusters of five sequences each. Note that we reduce <n_max_per_region> from the default 5 to 2
because when we run on pairs of sequences, we take the OR of their best matches.

To see the bayes factors that this is clustering on, you can add '--debug 1', although that\'s a bit verbose, so you can also add '--no_clean', which tells it not to clean up its
temp files in /tmp/$USER/hmms/$pid. The pairwise scores are in files with 'pairscores' in the name.

You can also run viterbi pair, to see how pairs of sequences line up (it should be apparent which of these three are from the same event, and which are not)
```
./runpart.py --point_estimate --pair --seqfile data/test/simu.csv --parameter_dir caches/parameters/hmm_parameters --debug 1 --queries ' -2996237710674838049:1722533587157826946:1516698783101435629'  # NOTE for negative hash, put a space after open single quote but before minus sign
```
Doing stuff that takes longer
==============

Now, the previous commands were reading info that I\'ve cached in ./caches, and simulated events in ./data/test/simu.csv, but you of course want to run from scratch.
I only cache it \'cause to get reliable parameter estimates you have to run over like ten thousand sequences, which takes a bit of time.

So now we go back and create the parameter files that we used in the previous steps, starting from some data in ./data/test/data.tsv

run on data to cache parameters and model files in `parameter_dir` (adjust <n_procs> as appropriate)
```
./runpart.py --cache_parameters --seqfile data/test/data.tsv --is_data --parameter_dir caches/new-parameters --n_procs 10 --skip_unproductive
```

this does the following:
  - runs smith-waterman on the data in data/test/data.tsv in order to estimate model parameters, which are written to `parameter_dir`/sw_parameters
  - read these parameters and use them to run the viterbi hmm on the same data file. This gives you better estimates of the parameters, which are written to `parameter_dir`/hmm_parameters
NOTE
  - ignore all the warnings about 'v right' (this tells you that a lot of the sequences have crappy best-v-matches, but you do not care a.t.m.)
  - should take ten or twenty minutes, maybe. the ham binary should be taking 100 percent cpu for most of that. We cannot really run this part with fewer
      seqs just as a test, because you need lots of seqs to get information on all the positions in all the gene versions
  - ignore all the warnings and errors about bad conserved codons and messed up cysteines and tryptophans. This is just telling you that there is unproductive rearrangements in the data.
    Well, and that my germline versions suck.

when it finishes, you can poke around in `parameter_dir`/hmm_parameters/ and see what is going into the model

Now you can run the simulator with these new parameters, too, to make six clusters/clones of five sequences each.
If you want to run the simulator without running the parameter-generation step above, replace 'new-parameters' with 'parameters'.
```
./runpart.py --simulate --parameter_dir caches/new-parameters/hmm_parameters --n_max_queries 30 --outdir caches/recombinator --debug 1
```
?
