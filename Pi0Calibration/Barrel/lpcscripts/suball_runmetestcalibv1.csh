#!/bin/csh

if ( ${#argv} != 6 ) then
echo "input ${#argv} "
echo "Usage: $0 totalEventRange step iter dataflag pizeta runLocally";
exit
endif

set username = `whoami`

set evtrange = $1
set step = $2
set iter = $3
set dataflag = $4
set pizEta = $5

set runLocally = $6


set r = 1
while ($r <= $evtrange)
set filename = calibres/testCalibv1.dflag$dataflag.pe$pizEta.step$step.iter$iter.r$r
if( -e $filename.txt) then
set filesize = `ls $filename.txt  -l |gawk '{print $5}' `
if( $filesize > 10000) then
set filesizeroot = `ls $filename.root  -l |gawk '{print $5}' `
if( $filesizeroot > 5000) then
echo "job already done $filename.txt"
@ r++
continue
endif
endif
endif

set filename = runme_testCalibv1.$r.$step.$dataflag.$pizEta.$iter.csh
set jobfile = jobs/$filename

cp runme_testCalibv1.csh $jobfile

perl ReplaceString.pl AAAA $dataflag $jobfile
perl ReplaceString.pl BBBB $pizEta $jobfile
perl ReplaceString.pl CCCC $step $jobfile
perl ReplaceString.pl DDDD $iter $jobfile
perl ReplaceString.pl EEEE $r $jobfile


##finally submit the jobs
if( $runLocally == 0  ) then
set condjob = jobs/condor_runmtestcalibv1.$r.$step.$dataflag.$iter.$pizEta
cp condor_runmtestcalibv1 $condjob
perl ReplaceString.pl filename $filename $condjob
echo "submitting jobssss"
condor_submit $condjob
else

checkJobs:
set nrootrun = `ps aux | grep root.exe | grep $username | wc |gawk '{print $1}'`
if( $nrootrun > 10) then
echo "toot many root running. wait 30 seconds.."
sleep 30
goto checkJobs
endif
#exit
csh $jobfile &

endif



@ r++
end
