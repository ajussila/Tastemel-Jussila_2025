#!/bin/bash

for name in {MT130_S11_R1,MT131_S12_R1,MT132_S13_R1,MT133_S14_R1}

do

perl ../../parseSeq1_from_fastq.pl <(zcat ../../raw_fastq/2iL_Day10/${name}_R1_001.fastq.gz) ../../raw_fastq/2iL_Day10/${name}_R1_001.fastq.gz

done
