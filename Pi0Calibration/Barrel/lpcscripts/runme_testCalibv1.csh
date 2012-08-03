#!/bin/csh


pwd

uname -a


set workdir = MyJobWorkingDirectory

set curdir = `pwd`
echo $curdir
cd $workdir
setenv SCRAM_ARCH slc5_amd64_gcc462
eval `scramv1 runtime -csh`
cd $curdir
pwd


which root

set dataflag = AAAA
set pizEta = BBBB
set step = CCCC
set iter = DDDD
set evtRange = EEEE

root -b <<EOF
gSystem->Load("$workdir/testCalibv1_C.so");
testCalibv1($dataflag,$pizEta,$step,$iter,$evtRange)
.qqqq
EOF


set filename = testCalibv1.dflag$dataflag.pe$pizEta.step$step.iter$iter.r$evtRange
echo "mv output..."
mv $filename.root  $workdir/calibres/
mv $filename.txt  $workdir/calibres/

echo "job done.."

