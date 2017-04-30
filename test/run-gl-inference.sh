cmd=./test/cf-germline-inference.py
label=test

# for tname in mfreq nsnp multi-nsnp prevalence n-leaves weibull; do
#     $cmd $tname --label $label # --plot
# done

glscmd="$cmd gls-gen --label $label"
for meth in partis full; do
    $glscmd --gls-gen-method $meth --n-tests 3 --n-procs-per-test 2 --no-slurm --gen-gset-events 1000 #--plot  # campanellabox
    # $glscmd --gls-gen-method $meth --n-tests 5 --n-procs-per-test 16 #--plot  # quoll
done
