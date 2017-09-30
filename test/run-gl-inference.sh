cmd=./test/cf-germline-inference.py
testopts="--n-tests 1 --n-procs-per-test 2 --no-slurm"

# label=XXX
# for tname in alcluster; do #mfreq nsnp multi-nsnp prevalence n-leaves weibull alcluster; do
#     # $cmd $tname --label $label $testopts  # laptop
#     $cmd $tname --n-tests 10 --label $label # --plot
# done

# label=v9
# glscmd="$cmd gls-gen --label $label"
# for diff in easy hard; do
#     for meth in full; do #simu partis full tigger-default igdiscover; do  # NOTE can add all methods to --methods arg now, i just keep 'em separate here so the log files are separate
# 	# for itest in {0..2}; do
# 	#     $glscmd --methods $meth --n-tests $((itest + 1)) --iteststart $itest --n-procs-per-test 12 --gls-gen-difficulty $diff --no-slurm  # --plot
# 	# done
# 	$glscmd --methods $meth --n-tests 3 --n-procs-per-test 20 --gls-gen-difficulty $diff --plot --plotcache
#     done
# done

# label=v10
# $cmd data --label $label $testopts --n-random-queries 5000  # laptop
# $cmd data --label $label --n-procs-per-test 15 --dry

# $cmd data --label $label --methods partis:tigger-default:igdiscover --plot --method-vs-method --plotcache
