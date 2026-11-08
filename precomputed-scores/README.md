# Precomputed pangenome complexity scores

Scores are computed for 50kb, 100kb, or 1Mb windows on hg38 in the files:

* `hprc-v1.1-mc-grch38_complexity_50000.tab`
* `hprc-v1.1-mc-grch38_complexity_100000.tab`
* `hprc-v1.1-mc-grch38_complexity_1000000.tab`

Columns give:

* chrom, start, end: genomic location of the window in hg38
* numnodes: number of minigraph-cactus nodes identified in the window
* total_length: sum of the lengths of all nodes in the window (note this will usually be close to but not exactly equal to the length of the hg38 window, since it includes the lengths of non-reference nodes)
* numwalks: number of walks identified through the subgraph. This is based on output of the `query` command from gbz-base and is in some cases more than the total number of assemblies used to build the graph.
* sequniq-normwalk:  sum_n  len(n)*p_n*(1-p_n)/L where L is the average walk length
* sequniq-normnode: sum_n  len(n)*p_n*(1-p_n)/L where L is the average node length

Scores of None indicate no walks were identified through the region.

## Computing scores

Scores were generated using the following commands:

1. Make windows of different sizes across hg38
```
for window in 50000 100000 1000000
do
	bedtools makewindows -g hg38.txt -w ${window} > windows/hg38_windows_${window}.bed
done
```

2. Compute complexity scores for each window based on the HPRC minigraph-cactus v1 graph

Scores are based on hprc-v1.1-mc-grch38 available here: https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38/hprc-v1.1-mc-grch38.gbz
```
for window in 50000 100000 1000000
do
	panct complexity \
		--region windows/hg38_windows_${window}.bed \
		--out hprc-v1.1-mc-grch38_complexity_${window}.tab \
		--metrics sequniq-normwalk,sequniq-normnode \
		../testdata/hprc-v1.1-mc-grch38.gbz
done
```