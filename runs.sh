#!/bin/bash
set -e
export PYTHONHOME=${HOME}/anaconda
export DYLD_INSERT_LIBRARIES=${HOME}/anaconda/lib/libpython2.7.dylib:${HOME}/anaconda/lib/libmkl_rt.dylib
for i in '1.5' '3'; do \
  for j in 'Electrons' 'Holes'; do \
  echo $i $j ${j}_${i}nm.pdf
  time devsim_py test_1d.py $i $j ${j}_${i}nm.pdf
done;
done;
