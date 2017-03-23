cmd=./test/cf-germline-inference.py
label=test
baseargs="--n-tests 10 --label $label"

$cmd mfreq       --mfreqs 0.1:1.0:2.0     --n-event-list 1000:2000:4000:8000  $baseargs --plot
$cmd nsnp        --nsnp-list 1:2:3:4      --n-event-list 1000:2000:4000:8000  $baseargs --plot
$cmd multi-nsnp  --nsnp-list 1,1:1,3:2,3  --n-event-list 1500:3000:6000:12000 $baseargs --plot
# $cmd prevalence   --prevalence-list 0.1:0.2:0.3  --n-event-list 1000:2000:4000:8000 $baseargs  #--plot
