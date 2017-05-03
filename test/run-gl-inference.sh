cmd=./test/cf-germline-inference.py
label=t2
testopts="--n-tests 1 --n-procs-per-test 2 --no-slurm"

# for tname in mfreq; do #mfreq nsnp multi-nsnp prevalence n-leaves weibull; do
#     # $cmd $tname --label $label $testopts  # laptop
#     $cmd $tname --label $label # --plot
# done

glscmd="$cmd gls-gen --label $label"
for meth in partis; do #partis full tigger; do  # NOTE can add all methods to --methods arg now
    # $glscmd --methods $meth $testopts --gen-gset-events 1000 --plot  # laptop
    $glscmd --methods $meth --n-tests 6 --n-procs-per-test 16 #--plot  # quoll
done
