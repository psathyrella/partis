for i in broken*yaml; do
  echo ">>>>>>>> "$i
  ../../../hample --seq 666655666613423414513666666666666 --hmmfname $i
  echo ""
done

