#!/bin/bash
cat ../datasets_contrasts.txt | tr -d '\015' |  while read var1 var2 var3; do
    sbatch -q express -N1-1 -n 1 --mem=5G -t 1000 --job-name=$var1 --wrap "R --vanilla --args $var1 $var2 $var3 < reQTL_plots.R "
sleep 0.5
done;
done
