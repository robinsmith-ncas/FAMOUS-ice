*ID KSH_SEG2A
*/ change actions of echo
*DECLARE mkexecxref
*D mkexecxref.162,mkexecxref.163
  echo -e "\n$section:\talldefs=$alldefs" | expand -t1,16
  echo -e   "$section:\tusedefs=$usedefs" | expand -t1,16
*D PXMKEXEC.8
    chardefs=`echo -e $alldefs | $UM_AWK '{print length}'`
*D PXMKEXEC.14    
    echo -e \
*D mkexecxref.176    
    echo -e \ 
*DECLARE qsconf
*I qsconf.11
gen_sed_strsfe()
{
  line_elements=$2
  prefstr="s^$1^"
  suffstr="^g"

  alltheobjs=$3

  (( counter = line_elements - 1 )) # Initialise counter
  bigstring=$prefstr                # Initialise bigstring

  for thisobj in $alltheobjs        # loop over all string entities
  do
    (( counter = counter + 1 ))
    if [ "$counter" -eq $line_elements ]
    then
     # start a new line
     counter=0
     if [ ! -z "$smallstring" ]
     then
       # smallstring exists, and has the maximum number of elements,
       # add it to bigstring
       bigstring=$bigstring'\\\\\\\\\\\\\\n'$smallstring
     fi
     # start a new line, reinitialise smallstring to current loop
     # element
     smallstring="$thisobj"
    else
     # add current loop element to line being composed
     smallstring="$smallstring $thisobj"
    fi
  done
  if [ ! -z "$alltheobjs" ]
  then
    echo -e $bigstring'\\\\\\\\\\\\\\n'$smallstring$suffstr
  else
    echo -e $bigstring$suffstr
  fi
}
*D GAV3U405.38
autoload fupper flower upper lower fnmycut fnlistopt \
*D qsconf.325,qsconf.326
  echo -e "\nqsconf($SECONDS): ***\tThis is          : \"$0\""
  echo -e "qsconf($SECONDS): ***\tUsing CL options : \"$cline\"\n"
*D qsconf.408  
          echo -e \
*D qsconf.417
          echo -e \
*D qsconf.430
    echo -e "qsconf($SECONDS): ***\tUsing filter     : \"$sectfilt\"\n"
*D qsconf.448,qsconf.453
  echo -e "qsconf($SECONDS): ***\tUM     dir is    : \"$UMDIR\""
  echo -e "qsconf($SECONDS): ***\tUM version is    : \"$VN\""
  echo -e "qsconf($SECONDS): ***\tOutput dir is    : \"$OUTDIR\""
  echo -e "qsconf($SECONDS): ***\tExecXref   is    : \"$EXECXREF\"" 
  echo -e "qsconf($SECONDS): ***\tObjXref    is    : \"$OBJXREF\""
  echo -e "qsconf($SECONDS): ***\tCOMPVARS   is    : \"$COMPVARS\"\n" 
*D qsconf.498
    echo -e "         Use option -help for more information.\n"
*D qsconf.504
    echo -e "         Use option -help for more information.\n"
*D qsconf.507
  echo -e "qsconf($SECONDS): ***\tConfiguring small executables"
*I GSH6U404.730
  ###### TEMP
  # echo -e ">>> loopexec $loopexec"
  # echo -e ">>> UM_GREP $UM_GREP"
  # echo -e ">>> EXECXREF $EXECXREF"
  # echo -e ">>> UM_AWK $UM_AWK"
  # echo -e ">>> sectfilt $sectfilt"
*D qsconf.520
    echo -e "qsconf($SECONDS): ***\t  $exec"
*D qsconf.547
      echo -e "qsconf($SECONDS): ***\t    nupdate (c-code)"
*D qsconf.567
      echo -e "qsconf($SECONDS): ***\t    source comparison"
*D qsconf.592
      echo -e "qsconf($SECONDS): ***\t    nupdate (fortran)"
*D qsconf.612
      echo -e "qsconf($SECONDS): ***\t    source comparison"
*D qsconf.632
    echo -e "qsconf($SECONDS): ***\t    creating Makefile"
*D qsconf.649
            echo -e \
*D qsconf.671
         echo -e \
*D qsconf.680
            echo -e \
*D qsconf.702
          echo -e \
*D GAV3U405.40
    echo -e `gen_sed_strsfe "@objects@" "4" "$objs"`   >> $SEDSCRIPT1
*D qsconf.734,qsconf.744
    echo -e "qsconf($SECONDS): ***\t  creating top Makefile"
    cp $MAKEFILEEXECINMID $MAKETMP
    for exec in $allexecs
    do
      echo "$exec:" >> $MAKETMP
      echo -e "\tcd $exec""_dir ;\c" >> $MAKETMP
      echo -e \
" \$(MAKE) FORT=\$(FORT) LOAD=\$(LOAD) CC=\$(CC)\n" >> $MAKETMP
      echo "install-$exec: $exec"  >> $MAKETMP
      echo -e "\tcd $exec""_dir ;\c" >> $MAKETMP
      echo -e " \$(MAKE) install\n" >> $MAKETMP
*D qsconf.760
  echo -e \
*D qsconf.778
    echo -e "qsconf($SECONDS): ***\tSetting switches"
*D qsconf.792
    echo -e "qsconf($SECONDS): ***\tFinished setting switches\n"
*D qsconf.797
    echo -e "qsconf($SECONDS): ***\tSetting general options"
*D qsconf.815
      echo -e "qsconf($SECONDS): ***\t  Will build COMP_GEN_OPT $switch"
*D qsconf.822
    echo -e "qsconf($SECONDS): ***\tFinished setting general options\n"
*D qsconf.827
    echo -e "qsconf($SECONDS): ***\tLooping over source sections"
*D qsconf.836
    echo -e "qsconf($SECONDS): ***\t  Configuring section $sect"
*D qsconf.869,qsconf.871
          echo -e "\nWARNING: No valid COMP_GEN_OPTS for this section."
          echo -e "WARNING: There may be an error in the ObjXref file."
          echo -e "WARNING: Continuing.\n"
*D qsconf.888
          echo -e "\nERROR: Unknown id $defmsg in $sect in ObjXref">&2
*D qsconf.948
        echo -e "qsconf($SECONDS): ***\t      nupdate (c-code)"
*D qsconf.970
        echo -e "qsconf($SECONDS): ***\t      source comparison"
*D qsconf.998
        echo -e "qsconf($SECONDS): ***\t      nupdate (fortran)"
*D qsconf.1020
        echo -e "qsconf($SECONDS): ***\t      source comparison"
*D qsconf.1054
            echo -e \
*D qsconf.1096
            echo -e \
*D qsconf.1172
    echo -e \
*D qsconf.1192
  echo -e "qsconf($SECONDS): ***\tFinished configuration\n"
*DECLARE qslistobj
*D GSH1U404.34
      echo -e \
*DECLARE qsmain
*D GEX4U402.26
      echo -e \
*D GKR1U404.46
    echo -e "\nqsmain($SECONDS): ***\tRemoving zero length files" >>\
*D GKR1U404.49
    echo -e "qsmain($SECONDS): ***\tFinished removing files\n" >>\
*D GEX4U402.49
    echo -e "\nqsmain($SECONDS): ***\tStarting fortran nupdate"
*D GEX4U402.60
    echo -e "qsmain($SECONDS): ***\tFinished fortran nupdate\n"
*D gex3u405.6
echo -e "\nqsmain($SECONDS): ***\tRenaming nupdate extracted files" >>\
*D gex3u405.8
echo -e "qsmain($SECONDS): ***\t  Renaming $file to $lfile" >>\
*D gex3u405.10
echo -e "qsmain($SECONDS): ***\tFinished renaming files.\n" >>\
*D GEX4U402.75 
      echo -e "qsmain($SECONDS): ***\tStarting c-code nupdate"
*D GEX4U402.87
      echo -e "qsmain($SECONDS): ***\tFinished c-code nupdate\n"
*D gex3u405.12
      echo -e \
*D gex3u405.15
        echo -e "qsmain($SECONDS): ***\t  Renaming $file to $lfile" >>\
*D gex3u405.17
      echo -e "qsmain($SECONDS): ***\tFinished renaming files.\n" >>\
*D gex3u405.19 
    echo -e "\nqsmain($SECONDS): ***\tRunning qsmncompile"
*D GKR3U405.77,GKR3U405.78
      echo -e "qsmain($SECONDS): ***\tCompleted qsmncompile\n"
      echo -e "\nqsmain($SECONDS): ***\tCompleted qsmncompile\n" >>$OUTPUT
*D GKR3U405.80,GKR3U405.81
      echo -e "qsmain($SECONDS): ***\tQsmncompile has failed. Exiting\n"
      echo -e "\nqsmain($SECONDS): ***\tQsmncompile has failed. Exiting\n"\
*D PXCLLMOD.4,PXCLLMOD.5
  echo -e "\nqsmain($SECONDS): ***\tRunning make for compile step"
  echo -e "\nqsmain($SECONDS): ***\tRunning make for compile step\n" \
*D PXCLLMOD.9,PXCLLMOD.10
    echo -e "qsmain($SECONDS): ***\tCompleted make for compile step\n"
    echo -e "\nqsmain($SECONDS): ***\tCompleted make for compile step\n" \
*D PXCLLMOD.12
    echo -e "qsmain($SECONDS): ***\tCompile step make has failed.\
*D PXCLLMOD.14
    echo -e "\nqsmain($SECONDS): ***\tCompile step make has failed.\
*D PXCLLMOD.18,PXCLLMOD.19
  echo -e "\nqsmain($SECONDS): ***\tRunning make for link step"
  echo -e "\nqsmain($SECONDS): ***\tRunning make for link step\n" \
*D PXCLLMOD.26,PXCLLMOD.27
    echo -e "qsmain($SECONDS): ***\tCompleted make for link step\n"
    echo -e "\nqsmain($SECONDS): ***\tCompleted make for link step\n" \
*D PXCLLMOD.30
    echo -e "qsmain($SECONDS): ***\tLink step make has failed.\
*D PXCLLMOD.32
    echo -e "\nqsmain($SECONDS): ***\tLink step make has failed.\
*DECLARE qsmncompile
*I GEX1U402.2
gen_sed_strsfe()
{
  line_elements=$2
  prefstr="s^$1^"
  suffstr="^g"

  alltheobjs=$3

  (( counter = line_elements - 1 )) # Initialise counter
  bigstring=$prefstr                # Initialise bigstring

  for thisobj in $alltheobjs        # loop over all string entities
  do
    (( counter = counter + 1 ))
    if [ "$counter" -eq $line_elements ]
    then
     # start a new line
     counter=0
     if [ ! -z "$smallstring" ]
     then
       # smallstring exists, and has the maximum number of elements,
       # add it to bigstring
       bigstring=$bigstring'\\\\\\\\\\\\\\n'$smallstring
     fi
     # start a new line, reinitialise smallstring to current loop
     # element
     smallstring="$thisobj"
    else
     # add current loop element to line being composed
     smallstring="$smallstring $thisobj"
    fi
  done
  if [ ! -z "$alltheobjs" ]
  then
    echo -e $bigstring'\\\\\\\\\\\\\\n'$smallstring$suffstr
  else
    echo -e $bigstring$suffstr
  fi
}

*D GAV3U405.44
autoload fnnewcom fnmydiff fnmycut fnlistopt
*D GEX1U404.14
echo -e "\nqsmncompile($SECONDS): ***\tCompile vars diff information."
*D GEX1U404.16
echo -e "qsmncompile($SECONDS): ***\t  System differences."
*D GEX1U404.19
  echo -e "qsmncompile($SECONDS): ***\t    $comtag"
*D GEX1U404.24
  echo -e "qsmncompile($SECONDS): ***\t  Last user run differences."
*D GEX1U404.27
    echo -e "qsmncompile($SECONDS): ***\t    $comtag"
*D GEX1U404.31
echo -e "qsmncompile($SECONDS): ***\tEnd of compile vars info."
*D GEX1U402.153
echo -e \
*D GEX1U402.164
echo -e "qsmncompile($SECONDS): ***\tFinished removing symlinks\n"
*D GEX1U402.174
echo -e "\nqsmncompile($SECONDS): ***\tRemoving zero length files"
*D GEX1U402.176
echo -e "qsmncompile($SECONDS): ***\tFinished removing files\n"
*D GEX1U402.192
echo -e "\nqsmncompile($SECONDS): ***\tRemoving previously modified files"
*D GEX1U402.200
    echo -e \
*D GEX1U402.207
        echo -e \
*D GSH1U403.262
        echo -e \
*D GEX1U402.212
        echo -e \
*D GEX1U402.216
        echo -e \
*D GEX1U402.232
  echo -e \
*D GEX1U402.238
echo -e "qsmncompile($SECONDS): ***\tFinished removing modified files\n"
*D GSH1U403.281
  echo -e \
*D GEX1U402.240
echo -e "\nqsmncompile($SECONDS): ***\tSTARTING control source comparison"
*D GSH1U403.288
  echo -e "qsmncompile($SECONDS): ***\t  $name\c"
*D GSH6U404.381
    echo -e "\t(Diff)\c"         # extracted vn. is different from user vn.
*D GSH6U404.383
    echo -e "\t(New) \c"           # no user vn. of the file
*D GSH1U403.300
    echo -e "\t(Same)\c"
*D GSH1U403.303
    echo -e "\t(F)\c"
*D GSH1U403.304
      echo -e "\t(Opts:new)\c"
*D GEX1U402.291
    echo -e \
*D GSH1U403.309
    echo -e "\t(C)\c"
*D GSH1U403.310
      echo -e "\t(Opts:new)\c"
*D GSH1U403.314
    echo -e "$name.o: $file\n\t$thiscc $thisopts -c $file\n" \
*D GEX1U402.324
    echo -e "ERROR: file type \"$type\" for \"$file\" is unknown."
*D GSH1U403.336
      echo -e "\t(O:rm)\c"
*D GEX1U402.348
echo -e "qsmncompile($SECONDS): ***\tFINISHED control source comparison\n"
*D GEX1U402.356
echo -e \
*D GSH1U404.5
  echo -e \
*D GSH1U403.345
    echo -e \
*D GSH1U403.355
        echo -e \
*D GSH1U403.372,GSH1U403.373
          echo -e "\nERROR: Call to script \"$CONFIGURE\" has failed."
          echo -e "       Contact UM team.\n"
*D GSH1U403.376
        echo -e \
*D GSH1U403.378
        echo -e \
*D GSH1U403.383,GSH1U403.384
          echo -e "\nERROR: make(1) has failed."
          echo -e "       Contact UM team.\n"
*D GSH1U403.387
        echo -e \
*D GSH1U403.408
    echo -e "\nqsmncompile($SECONDS): ***\t  Linking $symobjdir/$dir."
*D GSH1U403.414
      echo -e "qsmncompile($SECONDS): ***\t    $name\c"
*D GSH1U403.421
          echo -e "\t(F:ok)\c"
*D GSH1U403.425
            echo -e "\t(F:Sec)\c"
*D GSH1U403.428
            echo -e "\t(F:OSO)\c"
*D GSH1U403.431
            echo -e "\t(F:src)\c"
*D GSH1U403.433
            echo -e "\t(F:OSN)\c"
*D GSH1U403.445
          echo -e "\t(C:ok)\c"
*D GSH1U403.449
            echo -e "\t(C:Sec)\c"
*D GSH1U403.452
            echo -e "\t(C:src)\c"
*D GSH1U403.455
            echo -e "\t(C:OSO)\c"
*D GSH1U403.457
            echo -e "\t(C:OSN)\c"
*D GSH1U403.464
        echo -e "\n\nERROR: Source for $name not found in $symsrcdir/$dir"
*D GSH1U403.511
            echo -e "\t(O:ok)\c"
*D GSH1U403.515
            echo -e "\t(O:rm/ok)\c"
*D GSH1U403.541
              echo -e "\t(O:rm)\c"
*D GSH1U403.546
          echo -e "\t(O:com)\c"
*D GSH1U403.555
              echo -e "\t(Opts:ok)\c"
*D GSH1U403.576
            echo -e \
*D GSH1U403.586
              echo -e "\t(Opts:ok)\c"
*D GSH1U403.607
            echo -e \
*D GSH1U403.611
            echo -e \
*D GEX1U402.419
echo -e "\nqsmncompile($SECONDS): ***\tFinished symlink-ing object code\n"
*D GSH1U403.632
    echo -e "qsmncompile($SECONDS): ***\t  $sd"
*D PXCLLMOD.58
echo -e "\nqsmncompile($SECONDS): ***\
*D PXCLLMOD.62,PXCLLMOD.63
echo -e `gen_sed_strsfe "@objects@" "6" "$OBJS"`       >> $SEDSCRIPT
echo -e `gen_sed_strsfe "@preobjects@" "6" "$PREOBJS"` >> $SEDSCRIPT
*D PXCLLMOD.73
echo -e "qsmncompile($SECONDS): ***\tSed script generated\n"
*D PXCLLMOD.78
echo -e "\nqsmncompile($SECONDS): ***\
*D PXCLLMOD.98
echo -e "\nqsmncompile($SECONDS): ***\
*D PXCLLMOD.107
echo -e "qsmncompile($SECONDS): ***\tSed script generated\n"
*D PXCLLMOD.113
echo -e "\nqsmncompile($SECONDS): ***\
*D GEX1U402.459
echo -e "qsmncompile($SECONDS): ***\tFinished sed\n"
*DECLARE qsprelim
*D GSH4U404.10
       rm -f $RECONVARS
*I GSH4U404.42
#     echo "Trying qsconf"
#     echo "OUTPUT: $OUTPUT"
#     echo "URECONDIR: $URECONDIR"
#     echo "UPDATE_RECON: $UPDATE_RECON"
#     echo "UPDATE_RECONC: $UPDATE_RECONC"
*D GKR3U405.89
              -execs qxrecon_dump
*D GSH4U404.34
           echo -e "qsprelim ($SECONDS): ***\tCannot relink, \
*D GSHAU404.7
echo -e "\nqsprelim($SECONDS): ***\tRemoving zero length files"
*D GSHAU404.9
echo -e "qsprelim($SECONDS): ***\tFinished removing files\n"
*D GKR3U405.92
     echo -e "\nqsprelim: ***\tOutput from reconfiguration make\n" \
*DECLARE qsmaster
*D qsmaster.131
echo -e "\n\n\n"
*D qsmaster.140
echo -e "\n\n"
*D CJC1U404.33
  echo -e "\n\n"
*D CJC1U404.46
    echo -e "\n\n\n"
*D CJC1U404.55
    echo -e "\n\n"
*D qsmaster.164
echo -e "\n\n\n"
*D qsmaster.173
echo -e "\n\n"
*D qsmaster.204
echo -e "\n\n\n"
*D qsmaster.213
echo -e "\n\n"
*DECLARE gen_sed_string
*D gen_sed_string.106
    echo -e $bigstring'\\\\\\\\\\\\\\n'$smallstring$suffstr
*D gen_sed_string.108
    echo -e $bigstring$suffstr
*DECLARE getcomb
*D getcomb.93
echo -e "`ls -l $UPDEFS`\n"
*D getcomb.107
echo -e "\nstring:  $yes\n"
*D getcomb.111
  echo -e "$section: \c"
