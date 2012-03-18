#!/bin/csh


set workdir = /afs/cern.ch/cms/cit/yongy/data/CMSSW/v1/CMSSW_4_2_8/src/hltRate


set dataset = AAA
set evtRange = BBB

which root

root -b <<EOF
gSystem->Load("$workdir/testSelection_C.so");
testSelection("$dataset",$evtRange)
.qqq
EOF


set output = testSelection.$dataset.r$evtRange
scp $output.root  yangyong@pccityongnew.cern.ch:/localdata/yangyong/data/hltRate/res
scp $output.txt  yangyong@pccityongnew.cern.ch:/localdata/yangyong/data/hltRate/res
