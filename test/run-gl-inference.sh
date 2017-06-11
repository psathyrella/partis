cmd=./test/cf-germline-inference.py
label=v3
testopts="--n-tests 1 --n-procs-per-test 2 --no-slurm"

# for tname in alcluster; do #mfreq nsnp multi-nsnp prevalence n-leaves weibull alcluster; do
#     # $cmd $tname --label $label $testopts  # laptop
#     $cmd $tname --n-tests 10 --label $label # --plot
# done

# glscmd="$cmd gls-gen --label $label"
# for meth in partis; do #partis full tigger; do  # NOTE can add all methods to --methods arg now
#     # $glscmd --methods $meth $testopts --gls-gen-events 1000 #--plot  # laptop
#     $glscmd --methods $meth --n-tests 3 --n-procs-per-test 35  #--plot  # quoll
# done

# $cmd data --label $label $testopts --n-random-queries 5000  # laptop
# $cmd data --label $label --n-procs-per-test 15
