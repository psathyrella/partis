cmd=./test/cf-germline-inference.py
label=test

for tname in gen-gset; do #mfreq nsnp multi-nsnp prevalence n-leaves weibull; do
    $cmd $tname --label $label # --plot
done
