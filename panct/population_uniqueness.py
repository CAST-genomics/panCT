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

AVAILABLE_METRICS = ['popuniq-normwalk', 'popuniq-normnode','popuniq-normdegree']

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
        options: popuniq-normwalk, popuniq-normnode, popuniq-normdegree
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
    
    # assumes 'GRCh38','CHM13' if called from command line.
    exclude_samples=exclude_samples.split(',')
    log.info(f'Filtering out the following samples: {exclude_samples}.'
             "['GRCh38','CHM13', 'HG00272', 'HG03492'] recommended for pangenome v2.0")
    #TODO: make recommended exclusion list for v1
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
    header.extend(["numnodes","total_length", "numwalks"] + sorted([f'{metric}_{x}' for metric in metrics_list for x in list(set(asm.values()))+['total']]))
    #add mean degree?
    outf.write("\t".join(header) + "\n")



    ##### If GFA, just process the whole graph #####
    if file_type == "gfa":
        if region_str is not None:
            log.warning("Regions are ignored when processing GFA")
        exclude = []
        if reference != "":
            exclude = [reference]
        node_table = gutils.NodeTable(graph_file, exclude,walk_file)
        #do we need to exclude walks from link file? things to consider...
        if 'popuniq-normdegree' in metrics_list:
            link_table=gutils.LinkTable(graph_file,reference)
            for n in node_table.nodes:
                for l in link_table.links.keys():
                    if ((n==link_table.links[l].node_1) | (n==link_table.links[l].node_2)):
                        node_table.nodes[n].degree+=1   
        else:
            link_table=None

        metric_results = []
        for m in metrics_list:
            metric_results.extend(compute_population_uniqueness(node_table, asm, asm_count, m, exclude_samples))

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
        if 'popuniq-normdegree' in metrics_list:
            if (node_table.gfa_file!=None):
                link_table=gutils.LinkTable(node_table.gfa_file,reference)
            else:
                log.info('Node table does not contain gfa file, using gbz file for link table.')
                link_table = gbz.load_link_table_from_gbz(graph_file, region, reference, exclude_samples, walk_file)
                #put degree into node_table- how would be the best way to do this within the node_table class?
            for n in node_table.nodes:
                for l in link_table.links.keys():
                    if ((n==link_table.links[l].node_1) | (n==link_table.links[l].node_2)):
                        node_table.nodes[n].degree+=1                

        else:
            link_table=None
        metric_results = []
        log.info('computing population specific sequence uniqueness')
        for m in metrics_list:
            metric_results.extend(compute_population_uniqueness(node_table, asm, asm_count, m, exclude_samples))
        
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

def compute_population_uniqueness(node_table: gutils.NodeTable, asm, asm_count, metric:str,exclude_samples=['GRCh38','CHM13']):
    """
    Compute population specific uniqueness for a node table. Options:
    popuniq-normwalk
    popuniq-normnode
    popuniq-normdegree
    
    Parameters
    ----------
    node_table : graph_utils.NodeTable
       Stores info on lengths/walks through each node
    link_table : graph_utils.LinkTable
        Stores info on the links within the pangenome region
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
       List of population uniqueness scores, with 1 score per population, as well as a total score. Total score calculated using mean Hs for all populations.

    Raises
    ------
    ValueError
       If invalid metric specified
    
    """
    complexity=dict()
    populations=sorted(list(asm_count.keys()))+['total']
    #TOTAL HAS TO BE LAST- it is calculated on last loop of populations as an average of the previous values.
    if metric in ('popuniq-normwalk', 'popuniq-normnode','popuniq-normdegree'):
        complexity=dict.fromkeys([f'{metric}_{x}' for x in populations], 0)
        for n in node_table.nodes.keys():
            #get list of samples present for node
            pops = []
            pops.extend(
                asm[s.split('.')[0]]
                for s in node_table.nodes[n].samples
                if s.split('.')[0] not in exclude_samples)
                #list of populations present for node- take the population value from the node to pop dict
            
            #get dictionary of population instances
            anc_count = Counter(pops)
            anc_count = {key: anc_count.get(key, 0) for key in asm.values()}
            #count instances into dictionary
            
            #add attributes to class node
            node_table.nodes[n].anc_count = anc_count
            node_table.nodes[n].exp_het= calc_exp_het(asm_count, node_table.nodes[n].anc_count)
            length=node_table.nodes[n].length
            for k in populations:
                HT=node_table.nodes[n].exp_het['total']
                if k=='total':
                    HS=np.mean(list(node_table.nodes[n].exp_het.values()))
                    #mean is calculated after 0 limited Hs scores
                else:
                    HS=node_table.nodes[n].exp_het[k]
                #print(f'{k}: {HS}')
                if (HT==0):
                    node_table.nodes[n].Fst[k]=0 
                    # present in all samples therefore completely undifferentiated
                else:
                    FST=(HT-HS)/HT
                    #we have below 0 FST values- apparently known to be an issue from sample sizing problems
                    #Standard to set those to 0, so that's what we're doing
                    if FST<0:
                        FST=0
                    node_table.nodes[n].Fst[k]=FST
                
                #calculate degree from link table for degree normalized
                if metric=='popuniq-normdegree':
                    degree=0             
                    complexity[f'{metric}_{k}']+=node_table.nodes[n].degree*node_table.nodes[n].Fst[k]
                    
                else:
                    complexity[f'{metric}_{k}']+=length*node_table.nodes[n].Fst[k]
            ###
        if metric == 'popuniq-normwalk':
                complexity = {key: value / node_table.get_mean_walk_length() for key, value in complexity.items()}
        elif metric == 'popuniq-normnode':
            complexity = {key: value / node_table.get_mean_node_length() for key, value in complexity.items()}
        elif metric == 'sequniq-normdegree':
            complexity={key: value / node_table.get_mean_degree() for key, value in complexity.items()}
    return list(complexity.values())