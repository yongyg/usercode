#!/bin/csh

if ( ${#argv} != 7 ) then
echo "input ${#argv} "
echo "Usage: $0 totalEventRange step iterStart iterEnd dataflag pizeta runLocally ";
exit
endif

cmsenv
which root
set checkroot =  `which root |grep not |grep found |wc |gawk '{print $1}'`
if($checkroot == 1) then
echo "no ROOT! quit and cmsenv "
exit
endif


set username = `whoami`

set evtRangeTot =  $1
set step = $2
set iter = $3
set iterTotal = $4


set dataflag = $5 
set pizEta = $6
set runLocally = $7

while ($iter <= $iterTotal )
csh suball_runmetestcalibv1.csh $evtRangeTot $step $iter $dataflag $pizEta $runLocally

echo "checking if all jobs done. step $step iteration $iter"

set nChecked  = 0 
checkJob:
echo "checking job .."
set evtRange = 1
set nfound = 0 
echo `date`
set donejobrange = "r00"
while ($evtRange <= $evtRangeTot)

set checkifthisjobdone = -1 
set filesize = -1 
set filesizeroot = -1
set evtProcessed = 0 

set filename = calibres/testCalibv1.dflag$dataflag.pe$pizEta.step$step.iter$iter.r$evtRange

if( -e $filename.txt) then
set filesize = `ls $filename.txt -l |gawk '{print $5}'`
set checkifthisjobdone = `tail $filename.txt |grep totalEntries |wc |gawk '{print $1}'`
if($checkifthisjobdone == 1) then
set evtProcessed = `tail $filename.txt |grep totalEntries | gawk '{print $4}'`
endif
if(-e $filename.root) then
set filesizeroot = `ls $filename.root -l |gawk '{print $5}'`
endif
endif

##check job done correctly 
if( $filesize > 10000 && $filesizeroot > 5000 && $checkifthisjobdone > 0 && $evtProcessed > 0 ) then
set donejobrange = ${donejobrange}ar${evtRange}a
@ nfound ++ 
endif

@ evtRange ++ 
end

echo "checking job done " `date`

set njobrunning = `condor_q -submitter $username | grep R | grep runme_testCalibv1 | wc |gawk '{print $1}'`
set njobpending = `condor_q -submitter $username | grep I | grep runme_testCalibv1 | wc |gawk '{print $1}'`

if( $nfound != $evtRangeTot ) then
echo "$nfound jobs done.. sleep 30 seconds to check again njobR/I $njobrunning, $njobpending"
sleep 30

@ njobtobedone = $evtRangeTot - $nfound

echo "checked $nChecked times $njobtobedone jobs to be done"

#run locally now if <=2 jobs are still to be done after checking 100 times
if( $njobtobedone > 0 && $njobtobedone <= 2 && $nChecked >= 100 && $nChecked % 50 == 0 ) then
echo $donejobrange

set kk = 1
while ($kk <= $evtRangeTot)
set checkdone = `echo $donejobrange |grep r${kk}a |wc |gawk '{print $1}' `
if($checkdone < 1 ) then

checkJobs:
set nrootrun = `ps aux | grep root.exe | grep $username | wc |gawk '{print $1}'`
if( $nrootrun > 10) then
echo "toot many root running. wait 30 seconds.."
sleep 30
goto checkJobs
endif

echo "run interactively runme_testCalibv1.$evtRange.$step.$dataflag.$pizEta.$iter.csh"
csh jobs/runme_testCalibv1.$kk.$step.$dataflag.$pizEta.$iter.csh &
endif

@ kk++
end  ### end of run locally now

endif



@ nChecked ++

goto checkJob

endif 


sleep 5
echo "$nfound jobs done now deriving correction factor from step $step"
sleep 5
root -b <<EOF
gSystem->Load("deriveCalibConst_C.so")
deriveCalibConst($dataflag,$pizEta,$step,$iter,$evtRangeTot)
.qqqqqq
EOF

echo "correction derived"

@ iter++
end


