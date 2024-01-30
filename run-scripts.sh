#!/bin/bash

odir=$fs/partis/tmp/test-scripts
vsn=v0

for ft in csv fa yaml; do 
    ./bin/parse-output.py test/ref-results/partition-new-simu.yaml $odir/pout.$ft
done
./bin/parse-output.py test/paired/ref-results/partition-new-simu $fs/partis/tmp/potmp --paired
./bin/cf-alleles.py --bases all
./bin/cf-alleles.py --bases 8-51-1
./bin/cf-germlines.py test/ref-results/test/parameters/simu/hmm/germline-sets test/ref-results/test/parameters/data/hmm/germline-sets
./bin/compare-plotdirs.py --outdir $fs/partis/tmp/plots --plotdirs test/ref-results/annotate-new-simu-annotation-performance/hmm/mutation:test/ref-results/annotate-new-simu-annotation-performance/sw/mutation --names hmm:sw
./bin/partis plot-partitions --outfname test/ref-results/partition-new-simu.yaml --plotdir $fs/partis/tmp/tmp-plots --partition-plot-cfg mds
./bin/partis plot-partitions --outfname test/ref-results/partition-new-simu.yaml --plotdir $fs/partis/tmp/tmp-plots --partition-plot-cfg trees
./bin/plot-hmms.py --outdir $fs/partis/tmp/plots --infiles test/ref-results/test/parameters/data/hmm/hmms/IGHD1-20_star_01.yaml 

# TODO copy yaml file before running selection metrics
./bin/partis get-selection-metrics --outfname test/ref-results/partition-new-simu.yaml --tree-inference-method gctree
./bin/read-gctree-output.py --locus igh --gctreedir test/ref-results/partition-new-simu/gctree/iclust-0 --outdir $fs/partis/tmp/pout

rm -r $fs/partis/tmp/bpr-tmp
# single chain bcr-phylo
./bin/bcr-phylo-run.py --base-outdir $fs/partis/tmp/bpr-tmp --overwrite
./bin/smetric-run.py --infname $fs/partis/tmp/bpr-tmp/selection/simu/mutated-simu.yaml --base-plotdir $fs/partis/tmp/plots --metric-method lbi
echo ./bin/partis get-selection-metrics --outfname $fs/partis/tmp/bpr-tmp/selection/simu/mutated-simu.yaml --tree-inference-method iqtree

# paired bcr-phylo-run
./test/cf-paired-loci.py --label test-coar --version $vsn --n-replicates 2 --obs-times-list 15 --n-sim-seqs-per-generation-list 15 --n-sim-events-list 3 --bcr-phylo --perf-metrics naive-hdist --calc-antns --inference-extra-args=--no-indels --plot-metrics iqtree-coar --n-sub-procs 15 --n-max-procs 5 --single-light-locus igk --base-outdir /fh/fast/matsen_e/dralph/partis/paired-loci --actions simu:cache-parameters:partition:write-fake-paired-annotations:iqtree:iqtree-coar
