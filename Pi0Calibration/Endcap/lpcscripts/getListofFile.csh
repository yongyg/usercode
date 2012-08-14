#!/bin/csh

if ( ${#argv} != 4 ) then
echo "Usage: $0 filepath(/pnfs/cms/WAX/resilient/yangyong/data/pizdata) runlst(132440,132444) keywordlist(pi0Analyzer.pi02012Cv1RAW.GR_P_V39 ) outuptfilename(list.txt)";
exit
endif

set dir = $1

set keyword1 = `echo $3 |sed s/","/" "/g |gawk '{print $1}'`

if(-e $4) then
echo "$4 moved to $4.old"
mv $4 $4.old
endif

foreach run (`echo $2 |sed s/","/" "/g`)

echo $run

if(-e $4) then
ls $dir/$run |grep $keyword1 |gawk '{print $1}' >> $4
else
ls $dir/$run |grep $keyword1 |gawk '{print $1}' > $4
endif

end
