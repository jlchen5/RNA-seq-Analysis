#!/bin/bash

echo " Writed by Jiale Chen"
echo " trim_galore cut adapters started at $(date)"

cat sample.txt |while read id  ; do

if ([ ! -e ./$id'_1.fastq.gz' ] && [ ! -e ./$id'_2.fastq.gz' ]) ; then
  trim_galore -q 20  --paired $id'_1.fastq.gz' $id'_2.fastq.gz' --gzip
fi

if [ ! -e ./$id.fastq.gz] ; then
  trim_galore -q 20 --phred33 --gzip    
fi

echo "trim_galore cut adapters finished at $(date)"
       
    
done    
    
    
