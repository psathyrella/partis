cmd=./test/cf-germline-inference.py
label=test

for tname in weibull; do #mfreq nsnp multi-nsnp prevalence n-leaves; do
    $cmd $tname --label $label # --plot
done
