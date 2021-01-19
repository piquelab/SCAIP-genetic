#!/bin/bash
cat ../datasets_contrasts.txt | tr -d '\015' |  while read var1 var2 var3; do
    sbatch -q primary -N1-1 -n 5 --mem=80G -t 1000 --job-name=$var --wrap "R --vanilla --args $var1 $var2 $var3 < reQTL_lm_permutations.R"
sleep 0.5
done;
done
