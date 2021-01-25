#!/bin/bash
cat ../datasets.txt | tr -d '\015' |  while read var; do
    sbatch -q express -N1-1 -n 3 --mem=40G -t 1000 --job-name=$var --wrap "  R --vanilla --args 3 $var < compare_to_bulk.R "
sleep 5
done;
done
