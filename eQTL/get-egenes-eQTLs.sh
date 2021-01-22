#!/bin/bash
cat ../datasets.txt | tr -d '\015' |  while read var; do
    sbatch -q express -N1-1 -n 1 --mem=20G -t 1000 --job-name=$var --wrap "R --vanilla --args $var 3 < get-egenes-eQTLs.R "
sleep 0.5
done;
done
