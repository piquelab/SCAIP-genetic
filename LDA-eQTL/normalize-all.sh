#!/bin/bash
lda=$1
cat cell-types.txt | tr -d '\015' |  while read var; do
    echo "--- Sumbitting -$var---";
    sbatch -q primary -n 4 --mem=50G -t 10000 --wrap "R --vanilla --args $var $lda < normalize-all.R"
    sleep 0.5
done
