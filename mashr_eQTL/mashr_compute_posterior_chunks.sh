#!/bin/bash
chunks=$1
for i in $(seq 1 $chunks); do
    sbatch -q express -N1-1 -n 1 --mem=10G -t 1000 --job-name=$i --wrap "R --vanilla --args $i $chunks < mashr_compute_posterior_chunks.R"
sleep 0.5
done;
done
