#!/bin/csh

set dataset = MinimumBiasRun2011B-v1RAW

set runlocal = 0 

set r = 1 
while ($r <= 165)

csh sub_testSelection.csh $dataset $r $runlocal

@ r++
end
