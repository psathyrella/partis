# dependecies:
#   bppseqgen (*very* recent version necessary in order to allow per-residue mutation frequency specification)
#   pysam
#   optional:
#     ROOT
#     TreeSim (R package, needed to generate trees)

# ----------------------------------------------------------------------------------------
# installation

# stochhmm
git clone git@github.com:psathyrella/StochHMM -b working-psathyrella
cd StochHMM/
./configure
make  # ignore errors/warnings about autoconf being too old
# and... it'll crash 'cause it doesn't like c++0x
sed -i 's/^CXXFLAGS = \(.*\)/CXXFLAGS = \1 -std=c++0x -Wall/' Makefile src/Makefile  # add back in the cxx flags I *put* in there but stupid make/autoconf keep removing
make  # *then* finish building the stuff that uses c++0x

cd ..

# partis
git clone git@github.com:psathyrella/partis
cd partis/

# Run this command:
./runpart.py --cache_parameters --seqfile test/data.tsv --is_data --n_bases_skip 9 --v_right_length 56 --parameter_dir tmp/parameters --stochhmm_dir ../StochHMM
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

# Then run the simulator (using the parameter files you just made) to make six clusters/clones of five sequences each:
./runpart.py --simulate --parameter_dir tmp/parameters/hmm_parameters --n_max_queries 30 --outdir tmp/recombinator --debug 0 --stochhmm_dir ../StochHMM

# Then run viterbi on one of the simulated seqs you just made:
./runpart.py --point_estimate --seqfile tmp/recombinator/simu.csv --parameter_dir tmp/parameters/hmm_parameters --n_max_queries 1 --debug 1 --stochhmm_dir ../StochHMM

# and finally you can run the forward pair algorithm to partition the sequences simulated sequences:
./runpart.py --partition --pair --seqfile tmp/recombinator/simu.csv --parameter_dir tmp/parameters/hmm_parameters --stochhmm_dir ../StochHMM
# The last few rows of the output are the clusters. Since the simulator ran with clones of 5 seqs, and the test file was made with a '$ head -n30', there should
# be a five clusters of five, and one of four.

# ----------------------------------------------------------------------------------------
# you can also run viterbi pair, to see how pairs of sequences line up (it should be apparent which of these three are from the same event, and which aren't)
# NOTE I run this on test/simu.csv, which was *not* just made, because otherwise I don't know which query hashes to add to the --queries option
./runpart.py --point_estimate --pair --seqfile test/simu.csv --parameter_dir tmp/parameters/hmm_parameters --debug 1 --stochhmm_dir ../StochHMM --queries ' -4938873688098421278:8323195682841878031:40262812797573106'  # NOTE space before minus sign
