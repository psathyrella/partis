cmd=./test/cf-germline-inference.py
label=v9
testopts="--n-tests 1 --n-procs-per-test 2 --no-slurm"

# for tname in alcluster; do #mfreq nsnp multi-nsnp prevalence n-leaves weibull alcluster; do
#     # $cmd $tname --label $label $testopts  # laptop
#     $cmd $tname --n-tests 10 --label $label # --plot
# done


glscmd="$cmd gls-gen --label $label"
# for diff in hard; do
#     for meth in simu tigger-default igdiscover; do #partis full tigger-default tigger-tuned igdiscover; do  # NOTE can add all methods to --methods arg now
# 	$glscmd --methods $meth --n-tests 3 --n-procs-per-test 10 --gls-gen-difficulty $diff --no-slurm  #--dry  # --plot
#     done
# done

for diff in easy hard; do
    for meth in tigger-default igdiscover full; do #simu partis full tigger-default igdiscover; do  # NOTE can add all methods to --methods arg now, i just keep 'em separate here so the log files are separate
	# for itest in {0..2}; do
	#     $glscmd --methods $meth --n-tests $((itest + 1)) --iteststart $itest --n-procs-per-test 12 --gls-gen-difficulty $diff --no-slurm  # --plot
	# done
	echo $glscmd --methods $meth --n-tests 3 --n-procs-per-test 20 --gls-gen-difficulty $diff --plot # --plotcache
    done
done

label=v9
# $cmd data --label $label $testopts --n-random-queries 5000  # laptop
# $cmd data --label $label --n-procs-per-test 15
