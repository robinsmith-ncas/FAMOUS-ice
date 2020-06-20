*ID RRACC_SEG
*/ Change to slurm config for resubmission, SEG May 2019
*DECLARE qsresubmit
*D GSM0U403.2
cp $JOBDIR/qsubmit.qsub /tmp/umuisubmit.$$
*D GLW6U400.22
  sbatch  --output=$OUTFILE --job-name=$JOBNAME \
           /tmp/umuisubmit.$$
*D GLW6U400.23
   echo sbatch  --output=$OUTFILE --job-name=$JOBNAME  \
                 /tmp/umuisubmit.$$
