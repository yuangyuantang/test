#!/bin/bash

for file in *.top.csv
    do csvformat -T $file > ${file%%.*}.reformat.top.txt
    awk 'BEGIN{FS="\t"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"|"$6"|"$7"|"$8"|"$9"|"$10"|"$11}' ${file%%.*}.reformat.top.txt > ${file%%.*}.reformat.top.tmp.txt
    sed 's/Kingdom|Phylum|Class|Order|Family|Genus|Species/Species/g' ${file%%.*}.reformat.top.tmp.txt > ${file%%.*}.top.re.txt
    sed 's/unclassified/NA/g' ${file%%.*}.top.re.txt > ${file%%.*}.top.NA.txt
#   cp ${file%%.*}.top.NA.txt converted/${file%%.*}.reformat.top.csv
    sed 's/\t/,/g' ${file%%.*}.top.NA.txt > ${file%%.*}.reformat.top.csv
    rm *top.txt
    rm *tmp.txt
    rm *re.txt
    rm *top.NA.txt
done

