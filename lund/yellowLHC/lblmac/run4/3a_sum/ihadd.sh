# Usage:
# ./ihadd.sh list.txt


inlist=$1  #input txt file with list of filest to be hadd   files in list.txt are with full path!   
root2=`echo ${inlist%%.*}.root`  #output file
nsum=200   #group files together 

###########################################

sum1="qqqq"
sum2="XXXX.sh"
root1="QQQQ"

if [ -e $root2 ]; then
    rm $root2 
fi

totalNr="`cat $inlist | wc -l`"
i=0
k=0
###############################################
#          Distribute to scripts
for file in `cat $inlist`; do
   if [ $((i%nsum)) = "0" ]; then
       k=$((k+1)) 
       printf "hadd $root1$k.root" > $sum1$k.sh
   fi
   
   
   printf "  $file" >> $sum1$k.sh
   
   i=$((i+1)) 
done
##############################################
#             Run hadd scripts  
echo "=================Intermediate hadd =============="
inlist2="dummy2.sh"
ls $sum1*.sh > $inlist2 
for file in `cat $inlist2`; do
chmod +x $file
./$file
done
##############################################
#            hadd the results
echo "=================Final hadd =============="
inlist3="dummy3.sh"
ls $root1* > $inlist3 

printf "hadd $root2" > $sum2
for file in `cat $inlist3`; do
   
   printf "  $file" >> $sum2
   
done

chmod +x $sum2

./$sum2
##############################################
#    cleanup

rm $sum2
rm $inlist2
rm $inlist3
rm $sum1*
rm $root1*


