#!/usr/bin/bash

# Dirs to search for new subdirs to process. Put a space between paths.
DATADIRS="/fish/current_folder /datd/sl/fish_rawdata /date/sl/fish_rawdata /datc/sl/fish_rawdata /datb/sl/fish_rawdata"

# Status log file to create in each subdir.
PROCESSTATUSFILENAME="processing.log"

# Name of the file that tells that a subdir is ready to process.
TRANSFERDONEFILENAME="transferdone.txt"

# Wait time between checks for new jobs to process.
SLEEPMINUTES=60

# General log files for papermill
PAPERMILLOUTSUFFIX="_papermill.stdout"
PAPERMILLERRSUFFIX="_papermill.stderr"
PAPERMILLLOGSUFFIX="_papermill.log"

# Name of file that keeps increment for port number of scheduler
PORTNUMFILE="fish_queue_portnumber.txt"

DASHBOARDPORT=25399
