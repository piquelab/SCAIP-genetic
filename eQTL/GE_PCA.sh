#!/bin/bash

cat datasets.txt | tr -d '\015' |  while read var; do
    echo "--- Sumbitting -$var---";
    sbatch -q express -N1-1 -n 8 --mem=50G -t 1000 --job-name=$var --wrap "R --vanilla --args $var < GE_PCA.R"
    sleep 0.5
done
