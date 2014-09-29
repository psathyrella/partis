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
