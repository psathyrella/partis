for i in broken*yaml; do
  echo ">>>>>>>> "$i
  command="../../../hample --seq 666655666613423414513666666666666 --hmmfname $i"
  echo $command
  echo ""
  $command
  echo ""
done

