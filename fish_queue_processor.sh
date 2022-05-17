#!/usr/bin/bash

. fish_queue_config.sh

echo "NOTE: You have to manually do 'conda activate test_d' before starting this script."
echo "NOTE: Reading extra papermill parameters from fish_papermill_xparams.yaml"

while :
do
  sleepafter=1
  for datadir in $DATADIRS
  do
    timestamp=$(date +%y%m%d)
    for rawdir in $datadir/*
    do
      stdoutfile=$rawdir/logs/${timestamp}$PAPERMILLOUTSUFFIX
      stderrfile=$rawdir/logs/${timestamp}$PAPERMILLERRSUFFIX
      pmllogfile=$rawdir/logs/${timestamp}$PAPERMILLLOGSUFFIX
      ls $rawdir/logs/*$PAPERMILLOUTSUFFIX > /dev/null 2>&1
      stdoutls=$?
      if [[ -f "$rawdir/$TRANSFERDONEFILENAME" && "$stdoutls" -ne 0 ]] ; then
        subdir=${rawdir#/*/}
        rootdir=${rawdir%/$subdir}
        availdisk=$(df $rootdir | tail -1 | tr -s ' ' | cut -d' ' -f 4)
        if [[ "$availdisk" -le "4000000000" ]] ; then
          echo "$(date) WARNING: Probably too little HD space available on $rootdir ($availdisk)."
        fi
        if [[ ! -f $PORTNUMFILE ]] ; then
          echo 5050 > $PORTNUMFILE
        fi
        read portnum < $PORTNUMFILE
        portnum=$((portnum+1))
        echo $portnum > $PORTNUMFILE
        echo "$(date) INFO: Processing of $rawdir starting. Dashboard port: $DASHBOARDPORT"

        sed -i 's/Probe_FASTA/Probes_FASTA/g' /$rawdir/*_config.yaml > /dev/null 2>&1
        [[ -d "$rawdir/notebooks" ]] || mkdir $rawdir/notebooks
        [[ -d "$rawdir/logs" ]] || mkdir $rawdir/logs
        notebookfile=$rawdir/notebooks/${timestamp}-full-run.ipynb

        cmd="papermill -k test_d notebooks/Template_running_pysmFISH_pipeline.ipynb $notebookfile -p experiment_fpath $rawdir -p run_type new -p parsing_type original -p scheduler_port $portnum --start_timeout 6000 -p dashboard_port $DASHBOARDPORT -f fish_papermill_xparams.yaml --log-output --stdout-file $stdoutfile --stderr-file $stderrfile"
        echo "$(date) INFO: Command is $cmd"
        #echo "$(date) INFO: papermill output goes to $pmllogfile"
        $cmd # > $pmllogfile 2&>1
        exitcode=$?
        if [[ "$exitcode" -ne 0 ]] ; then
          echo "$(date) ERROR: papermill quit with exit code $exitcode"
          echo "$(date)        Log files are named $rawdir/logs/${timestamp}_xxx"
          [[ -f "$stdoutfile" ]] && echo "$(date)        You need to delete $stdoutfile to make the pipeline retry."
          echo "$(date) INFO: Now cleaning all started python processes and dask-worker-space:s..."
          killall -u simone python > /dev/null
          ssh monod10 'rm -Rf /tmp/dask-worker-space && killall -u simone python' > /dev/null
          ssh monod11 'rm -Rf /tmp/dask-worker-space && killall -u simone python' > /dev/null
          ssh monod12 'rm -Rf /tmp/dask-worker-space && killall -u simone python' > /dev/null
          ssh monod33 'rm -Rf /tmp/dask-worker-space && killall -u simone python' > /dev/null
          echo "$(date) INFO: Processing of $rawdir failed."
        else
          echo "$(date) INFO: Processing of $rawdir finished successfully."
          rm $rawdir/$TRANSFERDONEFILENAME
          sleepafter=0
        fi
      fi
    done
  done
  if [[ $sleepafter -eq 1 ]] ; then
    sleep $SLEEPMINUTES
  fi
done
