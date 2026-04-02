"""
Compute population specific uniqueness scores for regions
of a pangenome graph
"""

import time
import logging
from pathlib import Path
from typing import Optional

from collections import Counter
import numpy as np
import pandas as pd
# pandas needs to be added to the environment files

from .logging import getLogger
from . import gbz_utils as gbz
from .data import Region, Regions
from . import graph_utils as gutils

AVAILABLE_METRICS = ['popuniq-normwalk', 'popuniq-normnode']

def main(
    graph_file: Path,
    output_file: Path = Path("/dev/stdout"),
    region_str: str | Path = None,
    metrics: str = "sequniq-normwalk",
    reference: str = "GRCh38",
    exclude_samples: str ="GRCh38,CHM13",
    walk_file : Path= None,
    assemblies_file: Path = None,
    log: logging.Logger = None,
):
    """
    Compute population specific sequence uniqueness 
    scores for regions of a pangenome graph

    If a GFA file is given, compute complexity
    on the entire file.

    If a GBZ file is given, must specify a region
    (or file with list of regions)

    Parameters
    ----------
    graph_file : Path
        Path to GFA or GBZ file
    output_file : str, optional
        Path to output file
    region_str : str|Path, optional
        chrom:start-end of region to process or a BED file of regions
    metrics : str, optional
        Comma-separated list of metrics to compute.
    reference : str, optional
        Sample ID of reference
    walk_file : Path
        Path to associated walk file for assembly.
    assemblies_file : Path
        Path to a .tsv file that contains the columns 'Sample ID', 'Haplotype' and 'Population Abbreviation'.
        Used to assign samples to population groups.
        Assemblies file can be downloaded from the pangenome consortium data explorer.
    log : logging.Logger, optional
        Logger object

    Returns
    -------
    retcode : int
        Return code of the program
    """
    # log file doesn't currently have this.
    if log is None:
        log = getLogger(name="population_uniqueness", level="ERROR")
    start_time = time.time()

    #### Check files and indices #####

    #add the assemblies table
    file_type = None
    if graph_file.suffix == ".gfa":
        # TODO: also handle .gfa.gz
        file_type = "gfa"
    elif graph_file.suffix == ".gbz":
        file_type = "gbz"
        if not gbz.check_gbzbase_installed(log):
            return 1
        if not gbz.check_gbzfile(graph_file, log):
            return 1
    else:
        log.critical("Invalid graph type. Must be .gbz or .gfa")
        return 1

    if (not assemblies_file.exists()):
        log.critical('Assemblies file not found. Assemblies file must be provided' 
        'for sample to population mapping.')
        
    if (file_type == "gbz") and (not walk_file.exists()):
        log.critical('GBZ file provided but walk file not found. Walk file must be provided' 
        'to get accurate node to sample mapping.')

    
    #### Check requested metrics #####
    metrics_list = metrics.split(",")
    for m in metrics_list:
        if m not in AVAILABLE_METRICS:
            log.critical(f"Encountered invalid metric {m}")
            return 1
    #### Import assemblies file #####
    assemblies=pd.read_csv(assemblies_file,sep='\t',usecols=['Sample ID','Haplotype','Population Abbreviation'])
    #make filter an optional list of assemblies to remove.
    exclude_samples=exclude_samples.split(',')
    log.info(f'Filtering out the following samples: {exclude_samples}.'
             "['GRCh38','CHM13', 'HG00272', 'HG03492'] recommended for pangenome v2.0")
    #['GRCh38','CHM13', 'HG00272', 'HG03492'] for v 2
    assemblies=assemblies[~assemblies['Sample ID'].isin(exclude_samples)]
    #dictionary of sample sizes for each population
    asm_count=Counter(assemblies.drop_duplicates()['Population Abbreviation'])
    asm_count['total']=sum(asm_count.values())
    #dictionary of sample ID to population
    asm=assemblies[['Sample ID','Population Abbreviation']].drop_duplicates()
    asm.index=asm['Sample ID']
    asm=asm[["Population Abbreviation"]].to_dict(orient='dict')['Population Abbreviation']
    log.info(f'assemblies file read. {asm_count["total"]} assemblies to be analyzed.')

    ##### Set up output file #####
    outf = open(output_file, "w")
    header = []
    if file_type == "gbz":
        header = ["chrom", "start", "end"]
    #
    header.extend(["numnodes","total_length", "numwalks"] + sorted([f'{metric}_{x}' for metric in metrics_list for x in set(asm.values())]))
    outf.write("\t".join(header) + "\n")



    ##### If GFA, just process the whole graph #####
    if file_type == "gfa":
        if region_str is not None:
            log.warning("Regions are ignored when processing GFA")
        exclude = []
        if reference != "":
            exclude = [reference]
        node_table = gutils.NodeTable(graph_file, exclude,walk_file)
        metric_results = []
        for m in metrics_list:
            metric_results.extend(compute_population_uniqueness(node_table, asm, m))

        items = [
            len(node_table.nodes.keys()),
            node_table.get_total_node_length(),
            node_table.numwalks,
        ] + metric_results
        outf.write("\t".join([str(item) for item in items]) + "\n")
        outf.flush()
        end_time = time.time()
        total_time = end_time - start_time
        log.debug(f"Total time: \t{total_time}\n")
        outf.close()
        return 0

    #### If GBZ: Set up list of regions to process #####
    regions = []
    if region_str is not None:
        if isinstance(region_str, Path):
            regions = Regions.read(region_str, log=log)
        else:
            region = Region.read(region_str)
            regions = Regions((region,), log=log)
    if len(regions) == 0:
        log.critical("Did not detect any regions")
        return 1

    ##### Process each region #####
    for region in regions:
        log.info(
            "Processing region {chrom}:{start}-{end}".format(
                chrom=region.chrom, start=region.start, end=region.end
            )
        )
        # Load node table for the region
        
        node_table = gbz.load_node_table_from_gbz(graph_file, region, reference, exclude_samples, walk_file)
        
        metric_results = []
        log.info('computing population specific sequence uniqueness')
        for m in metrics_list:
            metric_results.extend(compute_population_uniqueness(node_table, asm, asm_count, m,exclude_samples))
        
        items = (
            [region.chrom, region.start, region.end]
            + [
                len(node_table.nodes.keys()),
                node_table.get_total_node_length(),
                node_table.numwalks            ]
            + metric_results
        )

        outf.write("\t".join([str(item) for item in items]) + "\n")
        outf.flush()
        
    ##### Cleanup #####
    end_time = time.time()
    time_per_region = (end_time - start_time) / len(regions)
    log.debug(f"Time per region\t{time_per_region}\n")
    outf.close()
    return 0

#see if this needs to be moved to graph utils
def calc_exp_het(asm_count,anc):
    """
    Helper function used to calculate the expected heterozygosity for a node in each population.
    
    asm_count: dict
        dictionary of haplotype instances from the assemblies table from human pangenome consortium.
    anc : list
        count of how many haplotypes from each population are present for node.
    """
    exp_het={}
    anc['total']=sum(anc.values())
    for k in asm_count.keys():
        p=anc[k]/asm_count[k]
        q=1-p
        exp_het[k]=2*p*q
    return(exp_het) 

def compute_population_uniqueness(node_table: gutils.NodeTable, asm,asm_count, metric:str,exclude_samples=[]):
    """
    Compute population specific uniqueness for a node table. Options:
    popuniq-normwalk
    popuniq-normnode
    
    Parameters
    ----------
    node_table : graph_utils.NodeTable
       Stores info on lengths/walks through each node
    asm: dict
        dictionary that maps sample ID to population. Based on assemblies file.
    asm_count: dict
        dictionary of the sample sizes for each population in the total assembly
    metric : str
       Which metric to compute. See description above
    exclude_samples: list
        List of samples to ignore for analysis (particularly those in assembly file that aren't in the assembly.)
    Returns
    -------
    complexity : float
       Complexity score

    Raises
    ------
    ValueError
       If invalid metric specified
    
    """
    complexity=dict()
    if metric in ('popuniq-normwalk', 'popuniq-normnode'):
        complexity=dict.fromkeys(sorted([f'{metric}_{x}' for x in set(asm.values())]), 0)
        for n in node_table.nodes.keys():
            #get list of samples present for node
            pops = []
            pops.extend(
                asm[s.split('.')[0]]
                for s in node_table.nodes[n].samples
                if not filter or s.split('.')[0] not in exclude_samples)            
            #get dictionary of population instances
            anc_count = Counter(pops)
            anc_count = {key: anc_count.get(key, 0) for key in asm.values()}
            
            #add attributes to class Node
            node_table.nodes[n].anc_count = anc_count
            node_table.nodes[n].exp_het= calc_exp_het(asm_count, node_table.nodes[n].anc_count)
            
            length=node_table.nodes[n].length
            
            for k in pops:
                HT=node_table.nodes[n].exp_het['total']
                HS=node_table.nodes[n].exp_het[k]
                if (HT==0):
                    node_table.nodes[n].Fst[k]=0 
                    # present in all samples therefore completely undifferentiated
                else:
                    FST=(HT-HS)/HT
                    #we have below 0 FST values- apparently known to be an issue from sample sizing problems
                    #Our samples are itty bitty so makes sense
                    if FST<0:
                        FST=0
                    node_table.nodes[n].Fst[k]=FST
                #print(f'FST {n} - {k}: {node_table.nodes[n].Fst[k]}')
                complexity[f'{metric}_{k}']+=length*node_table.nodes[n].Fst[k]
                
        if metric == 'popuniq-normwalk':
                complexity = {key: value / node_table.get_mean_walk_length() for key, value in complexity.items()}
        elif metric == 'popuniq-normnode':
            complexity = {key: value / node_table.get_mean_node_length() for key, value in complexity.items()}
    return list(complexity.values())