*ID F90_SEG1
*/ grab the sed script
*DECLARE qsmncompile
*I PXCLLMOD.72
SEG_TMP=$TEMP/seg_tmp_$$
cp $SEDSCRIPT $SEG_TMP
SEG_TMPg=$TEMP/seg_tmpg_$$
echo "$OBJS1 $OBJS2" > $SEG_TMPg
