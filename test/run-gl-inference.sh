cmd=./test/cf-germline-inference.py
label=v2
testopts="--n-tests 1 --n-procs-per-test 2 --no-slurm"

# for tname in mfreq; do #mfreq nsnp multi-nsnp prevalence n-leaves weibull; do
#     # $cmd $tname --label $label $testopts  # laptop
#     echo $cmd $tname --n-tests 10 --label $label # --plot
# done

# glscmd="$cmd gls-gen --label $label"
# for meth in partis; do #partis full tigger; do  # NOTE can add all methods to --methods arg now
#     # $glscmd --methods $meth $testopts --gls-gen-events 1000 #--plot  # laptop
#     $glscmd --methods $meth --n-tests 3 --n-procs-per-test 35  #--plot  # quoll
# done

# $cmd data $testopts --label $label --n-random-queries 5000
