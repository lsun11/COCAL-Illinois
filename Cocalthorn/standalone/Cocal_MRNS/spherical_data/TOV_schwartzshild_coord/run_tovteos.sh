STARTTIME=$(date +%s)

./tovteos

ENDTIME=$(date +%s)
echo "----------------------------------------------------------------"
echo "It took $(($ENDTIME - $STARTTIME)) seconds to complete this run."

