cmd=./test/cf-germline-inference.py
testopts="--n-tests 1 --n-procs-per-test 2 --no-slurm"

# # varvals
# label=gls-gen-paper-v2
# ntests=3  #50
# for tname in mfreq nsnp multi-nsnp prevalence n-leaves weibull; do  #  alcluster
#     # $cmd $tname --label $label $testopts  # laptop
#     $cmd $tname --n-tests $ntests --label $label --plot &
# done

# # gls-gen
# label=gls-gen-paper-v11  # test-macaque
# glscmd="$cmd gls-gen --label $label"
# for diff in easy; do  # easy hard; do
#     for meth in tigger-default; do #simu partis full tigger-default igdiscover; do  # NOTE can add all methods to --methods arg now, i just keep 'em separate here so the log files are separate
# 	for itest in 0; do #{0..3}; do
# 	    echo $glscmd --methods $meth --n-tests $((itest + 1)) --iteststart $itest --n-procs-per-test 10 --gls-gen-difficulty $diff  #  --species macaque
# 	done
# 	# $glscmd --methods $meth --n-tests 1 --gls-gen-difficulty $diff --plot # --plotcache
#     done
# done

# data
label=gls-gen-paper-v11
# label=v+1  # just for running on davide-gl-valid
# $cmd data --label $label $testopts --n-random-queries 5000  # laptop
# $cmd data --label $label --n-procs-per-test 15 --dry
# $cmd data --label $label --methods partis:tigger-default:igdiscover --plot --method-vs-method #--plotcache --sample-vs-sample
# $cmd data --label $label --write-zenodo-files --methods igdiscover:tigger-default:partis
