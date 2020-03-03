for k in {1..184}
do
qsub Rjob_AVbound_k.txt -P dt17 -v var1=$k
done
