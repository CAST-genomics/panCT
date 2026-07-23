#!/bin/bash

# Get complexity for an example region
panct complexity \
	--region chr11:119077050-119178859 --out test.tab \
	--metrics sequniq-normwalk,sequniq-normnode \
	../testdata/hprc-v1.1-mc-grch38.gbz

# Make windows across hg38
bedtools makewindows -g hg38.txt -w 100000 > hg38_windows_100kb.bed

# Get complexity for each window
panct complexity \
	--region hg38_windows_100kb.bed --out hg38_100kb.tab \
	--metrics sequniq-normwalk,sequniq-normnode \
	../testdata/hprc-v1.1-mc-grch38.gbz

# Test out for list of GWAS hits
panct complexity \
	--region trubetskoy_scz_hg38.bed --out hg38_100kb_trubetskoy.tab \
	--metrics sequniq-normwalk,sequniq-normnode \
	../testdata/hprc-v1.1-mc-grch38.gbz
