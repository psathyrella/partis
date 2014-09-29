# dependecies:
stochhmm (must be this fork: https://github.com/psathyrella/StochHMM/tree/working-psathyrella)
bppseqgen (*very* recent version necessary in order to allow per-residue mutation frequency specification)
TreeSim (R package)
ROOT

# ----------------------------------------------------------------------------------------
# installation

# stochhmm
git clone git@github.com:psathyrella/StochHMM -b working-psathyrella
./configure
make  # ignore errors/warnings about autoconf being too old
sed -i 's/^CXXFLAGS = \(.*\)/CXXFLAGS = \1 -std=c++0x -Wall/' Makefile src/Makefile  # add back in the cxx flags I *put* in there but stupid make/autoconf keep removing
make  # *then* finish building the stuff that uses c++0x

# partis
git clone git@github.com:psathyrella/StochHMM -b working-psathyrella

# Run this command:
./runpart.py --cache_parameters --seqfile test/data.tsv --is_data --n_bases_skip 9 --v_right_length 56 --parameter_dir tmp/parameters
# this does the following:
#   1) runs smith-waterman on the data in test/data.tsv in order to estimate model parameters, which are written to tmp/parameters/sw_parameters
#   2) read these parameters and use 'em to run the viterbi hmm on the same data file. This gives you better estimates of the parameters, which
#      are written to tmp/parameters/hmm_parameters
# NOTE
#   - ignore all the warnings about 'v right' (this tells you that a lot of the sequences have crappy best-v-matches, but you don't care a.t.m.)
#   - should take ten or twenty minutes, maybe. the stochhmm binary should be taking 100 percent cpu for most of that. We can't really run this part with fewer
#       seqs just as a test, because you need lots of seqs to get information on all the positions in all the gene versions
#   - ignore all the warnings and errors about bad conserved codons and messed up cysteines and tryptophans. This is just telling you that there's unproductive rearrangements in the data.
        
# when it finishes, you can poke around in tmp/parameters/hmm_parameters/ and see what's going into the model

# You can then run viterbi on (one of, unless you increase <n_max_queries>) the simulated seqs in test/simu.csv:
./runpart.py --point_estimate --seqfile test/simu.csv --parameter_dir tmp/parameters/hmm_parameters --n_max_queries 1 --debug 1

# you can also run viterbi pair, to see how pairs of sequences line up
./runpart.py --point_estimate --pair --seqfile test/simu.csv --parameter_dir tmp/parameters/hmm_parameters --debug 1

# and finally you can run the forward pair algorithm to partition the sequences in test/simu.csv
./runpart.py --partition --pair --seqfile test/simu.csv --parameter_dir tmp/parameters/hmm_parameters --debug 1
