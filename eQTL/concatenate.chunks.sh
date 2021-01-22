cat datasets.txt | tr -d '\015' |  while read var; do
for i in $(seq 0 10); do
for j in $(seq 1 30); do
     zcat eQTL_output/$var.GEPC$i.nominals.chunk$j.txt.gz
done | gzip -c > eQTL_output/$var.GEPC$i.nominals.eQTL.txt.gz;
done;
done
