#!/bin/bash
binL=(5 10 15 20 30 40 60 80 110 140 180 220 270 320)
binH=(10 15 20 30 40 60 80 110 140 180 220 270 320 400)
TTl=8
TTh=9



#chage the time of processing
for (( ihb=0; ihb<${#binL[@]}; ihb++ )); do

   find "/project/projectdirs/alice/krizek/JEWEL_PbPb_norecoil/jewel" -name hjet_*_TT${TTl}_${TTh}.root | grep "jewel/${binL[$ihb]}_${binH[$ihb]}_" > hb${binL[$ihb]}_${binH[$ihb]}_TT${TTl}_${TTh}.txt 
done


