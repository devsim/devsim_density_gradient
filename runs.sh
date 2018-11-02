#!/bin/bash
# use python 2.7 environment
set -e
for i in '1.5' '3'; do \
  for j in 'Electrons' 'Holes'; do \
  echo $i $j ${j}_${i}nm.pdf
  time python test_1d.py $i $j ${j}_${i}nm.pdf
done;
done;
