#!/bin/bash

samtools view -bT $1 $2 > /tmp/tmp.bam
samtools sort /tmp/tmp.bam > $3
rm /tmp/tmp.bam
samtools index $3
