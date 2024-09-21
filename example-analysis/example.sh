#!/bin/bash

# Get complexity for an example region
panct complexity --gbz ../testdata/hprc-v1.1-mc-grch38.gbz \
	--region chr11:119077050-119178859 --out test.tab \
	--metrics sequniq-normwalk,sequniq-normnode

# Make windows across hg38
bedtools makewindows -g hg38.txt -w 100000 > hg38_windows_100kb.bed

# Get complexity for each window
panct complexity --gbz ../testdata/hprc-v1.1-mc-grch38.gbz \
	--region-file hg38_windows_100kb.bed --out hg38_100kb.tab \
	--metrics sequniq-normwalk,sequniq-normnode

# Test out for list of GWAS hits
panct complexity --gbz ../testdata/hprc-v1.1-mc-grch38.gbz \
	--region-file trubetskoy_scz_hg38.bed --out hg38_100kb_trubetskoy.tab \
	--metrics sequniq-normwalk,sequniq-normnode
