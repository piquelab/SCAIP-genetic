#!/bin/bash
lda=$1
cat ../datasets.txt | tr -d '\015' |  while read var; do
    echo "--- Sumbitting -$var---";
    sbatch -q primary -N 1-1 -n 2 --mem=50G -t 10000 --job-name=$var --wrap "R --vanilla --args $var $lda < get_dosages.R"
    sleep 0.5
done
