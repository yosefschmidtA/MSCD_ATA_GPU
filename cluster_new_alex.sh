#!/bin/bash

arquivo=mscdlist.txt;

#Limites

xo=-10.7511;
xf=10.2534;


yo=-10.7651;
yf=10.7651;

zo=-10.3;
zf=10.3;

grep "no katom emiter x y z" $arquivo >> temp.txt;

tipo=(`awk '{print $2}' temp.txt | xargs`);
x=(`awk '{print $4}' temp.txt | xargs`);
y=(`awk '{print $5}' temp.txt | xargs`);
z=(`awk '{print $6}' temp.txt | xargs`);

linhas=`cat temp.txt | wc -l`;
linhas=`echo "$linhas -1" | bc`;

rm -rf temp.txt;

#echo "$linhas +1" | bc > temp.txt;
#echo "" >> temp.txt;

for i in `seq 0 $linhas`;
       do
        tipo[$i]=`echo "${tipo[$i]}" | sed 's/1/Co/'`;
	tipo[$i]=`echo "${tipo[$i]}" | sed 's/2/Ga/'`;
	tipo[$i]=`echo "${tipo[$i]}" | sed 's/3/Ti/'`;
	tipo[$i]=`echo "${tipo[$i]}" | sed 's/4/O/'`;

        varX=`echo "(${x[$i]} >= $xo) && (${x[$i]} <= $xf)" | bc`;
          if [ "1" = "$varX" ]; then
             varY=`echo "(${y[$i]} >= $yo) && (${y[$i]} <= $yf)" | bc`;
               if [ "1" = "$varY" ]; then 
                  varZ=`echo "(${z[$i]} >= $zo) && (${z[$i]} <= $zf)" | bc`;
                     if [ "1" = "$varZ" ]; then
                        echo "${tipo[$i]}  ${x[$i]}  ${y[$i]}  ${z[$i]}" >> temp.txt; 
                        echo "Construindo linha $i de $linhas";
                     fi 
               fi
          fi


       done;

linhas=`cat temp.txt | wc -l`;
echo "$linhas" | bc > cluster.xyz;
echo "" >> cluster.xyz;

awk '{printf "%2s\t%8s\t%8s\t%8s\n",$1,$2,$3,$4}' temp.txt >> cluster.xyz;
rm -rf temp.txt;
echo "Feito!";
