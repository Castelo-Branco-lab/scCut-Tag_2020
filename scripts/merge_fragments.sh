#!/bin/bash

for file in $@; do
  bgzip -cd $file > ${file/.gz/}_temp
done



cat ${@/.gz/}_temp | sort -k1,1 -k2,2n
rm  ${@/.gz}_temp