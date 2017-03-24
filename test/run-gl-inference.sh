cmd=./test/cf-germline-inference.py
label=test

for tname in mfreq nsnp multi-nsnp prevalence; do  # n-leaves
    $cmd $tname --label $label # --plot
done
