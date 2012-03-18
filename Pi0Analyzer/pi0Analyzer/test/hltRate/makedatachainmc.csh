#!/bin/csh

foreach dataset (`cat allrun`)

#echo $dataset
#alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'


if(-e tmplist.txt) then
mv tmplist.txt tmplist.txtold
endif


set datadir = /castor/cern.ch/user/y/yangyong/data/pizdata/datav1

#rfdir $mycastor/data/Run2011A/$dataset |gawk '{print $9}' > tmplist.txt
rfdir $datadir/$dataset |gawk '{print $9}' >  tmplist.txt
set nfile = `rfdir $datadir/$dataset |grep root |wc |gawk '{print $1}'`


if($nfile < 1) then
continue
endif


#@ kk = $nfile / 2 

@ kk = 1000  ## all together


#if($nfile < 10) then
#set kk = $nfile
#endif


perl split_v1.pl  $datadir/$dataset tmplist.txt $kk 1
 

echo "else if(datasetname ==" \"$dataset\""){"
cat temp_tmplist.txt.txt 
echo "}"


end
