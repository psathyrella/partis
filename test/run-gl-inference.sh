cmd=./test/cf-germline-inference.py
label=v8
testopts="--n-tests 1 --n-procs-per-test 2 --no-slurm"

# for tname in alcluster; do #mfreq nsnp multi-nsnp prevalence n-leaves weibull alcluster; do
#     # $cmd $tname --label $label $testopts  # laptop
#     $cmd $tname --n-tests 10 --label $label # --plot
# done


glscmd="$cmd gls-gen --label $label"
for diff in easy hard; do
    for meth in simu tigger-default igdiscover; do #partis full tigger-default tigger-tuned igdiscover; do  # NOTE can add all methods to --methods arg now
	$glscmd --methods $meth --n-tests 3 --n-procs-per-test 10 --gls-gen-difficulty easy --no-slurm  #--dry  # --plot
    done
done
# $cmd data --label $label $testopts --n-random-queries 5000  # laptop
# $cmd data --label $label --n-procs-per-test 15
