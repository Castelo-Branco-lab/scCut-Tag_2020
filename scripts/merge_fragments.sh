#!/bin/bash

for file in $@; do
  bgzip -cd $file > ${file/.gz/}
done

cat ${@/.gz} 
rm  ${@/.gz} 