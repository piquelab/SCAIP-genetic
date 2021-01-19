#!/bin/bash
cat ../datasets.txt | tr -d '\015' |  while read var; do
for i in $(seq 0 10); do
    sbatch -q primary -N1-1 -n 1 --mem=20G -t 100 --job-name=$i$var --wrap "R --vanilla --args $i $var < process_nominals-auto.R"
sleep 0.5
done;
done;
done
