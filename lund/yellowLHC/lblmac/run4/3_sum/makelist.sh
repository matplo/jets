#!/bin/bash
binL=(5 10 15 20 30 40 60 80 110 140 180 220 270 320)
binH=(10 15 20 30 40 60 80 110 140 180 220 270 320 400)
TTl=20
TTh=25



#chage the time of processing
for (( ihb=0; ihb<${#binL[@]}; ihb++ )); do

   find "/project/projectdirs/alice/krizek/PP5/jewel" -name hjet_*_TT${TTl}_${TTh}.root | grep "jewel/${binL[$ihb]}_${binH[$ihb]}_" > pp5_hb${binL[$ihb]}_${binH[$ihb]}_TT${TTl}_${TTh}.txt 
done


