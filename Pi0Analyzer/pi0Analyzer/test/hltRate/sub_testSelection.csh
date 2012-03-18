#!/bin/csh

if ( ${#argv} != 3 ) then
echo "Usage: $0 datasetname evtRange runLocally"
exit
endif


set dataset = $1
set evtRange = $2
set runLocal = $3

set jobfile = jobs/runm_testSelection.$dataset.$evtRange.csh

cp runm_testSelection.csh $jobfile
perl ~/ReplaceString.pl AAA $dataset $jobfile
perl ~/ReplaceString.pl BBB $evtRange $jobfile

#exit


if( $runLocal == 1) then
checkJobs:
set njobcopying = `ps aux | grep root.exe |grep yangyong | wc |gawk '{print $1}'`
if($njobcopying > 11) then
echo "more than 10 root jobs wait 30 seconds."
sleep 30 
goto checkJobs
endif

csh $jobfile &

else

checkbjobs:
set nrun = `bjobs |grep cmscaf1nd | grep RUN |wc |gawk '{print $1}'`
set npend = `bjobs |grep cmscaf1nd | grep PEND |wc |gawk '{print $1}'`
echo "nrun/npend",$nrun,$npend
if( $nrun > 150 || $npend > 100) then
echo "too many bjobs run/pend ",$nrun,$npend
sleep 30 
goto checkbjobs
endif


bsub -q cmscaf1nd -J test <  $jobfile
endif

