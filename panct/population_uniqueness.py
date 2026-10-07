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
import os
# pandas needs to be added to the environment files

from .logging import getLogger
from . import gbz_utils as gbz
from .data import Region, Regions
from . import graph_utils as gutils

AVAILABLE_METRICS = ['popuniq-normwalk', 'popuniq-normnode','popuniq-normdegree']


_REFERENCE_ALIASES = {
    'GRCh38': 'GRCh38.0',
    'CHM13': 'CHM13.0',
}


def _expand_reference_aliases(samples):
    """
    Given a list of sample IDs to exclude, add the equivalent GRCh38/GRCh38.0
    and CHM13/CHM13.0 alias for any of those two names that appear, so that
    exclusion works regardless of which form is used in the assemblies file
    or on the command line.
    """
    reverse_aliases = {v: k for k, v in _REFERENCE_ALIASES.items()}
    expanded = set(samples)
    for s in list(expanded):
        if s in _REFERENCE_ALIASES:
            expanded.add(_REFERENCE_ALIASES[s])
        elif s in reverse_aliases:
            expanded.add(reverse_aliases[s])
    return list(expanded)


def main(
    graph_file: Path,
    output_file: Path = Path("/dev/stdout"),
    region_str: str | Path = None,
    metrics: str = "popuniq-normwalk",
    reference: str = "GRCh38",
    exclude_samples: str ="GRCh38,CHM13",
    walk_file : Path= None,
    assemblies_file: Path = None,
    log: logging.Logger = None,
    skip_highnode=False,
    memory_limit_gb=None
):
    """
    Compute population specific sequence uniqueness 
    scores for regions of a pangenome graph

    If a GFA file is given, compute population uniqueness
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
        Path to a .tsv file that contains, at minimum, the columns
        'Sample ID' and 'Population Abbreviation'. Used to assign samples
        to population groups. Two formats are supported:

        "haplotype-removed" format:'Sample ID' values do NOT include a 
        haplotype suffix (e.g. "HG00097"). Both haplotypes of an individual
        are assumed to belong to the same population.
        
        "haplotype-specific" format: no 'Haplotype' column; 'Sample ID'
        values already carry a haplotype suffix (e.g. "HG00097.1"), so
        each haplotype of an individual can be assigned to a different
        population (or excluded independently).

        The format is auto-detected based on whether any 'Sample ID' value
        contains a ".". Rows whose 'Population Abbreviation' is "drop" are
        excluded from analysis, in addition to any samples given via
        `exclude_samples`. "GRCh38"/"GRCh38.0" and "CHM13"/"CHM13.0" are
        treated as equivalent names when applying exclusions.
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

    #inform that skipping certain regions
    if skip_highnode:
        log.info('skip value True passed; therefore, regions with > 5e5 nodes will be skipped.')
    if memory_limit_gb:
        log.info(f"memory limit of {memory_limit_gb} passed; "
        "therefore, regions requiring high memory will be skipped.")

    skipped_bed_file = output_file.with_name(output_file.stem + "_skipped_regions.bed")
    
    #### Import assemblies file #####
    assemblies_full = pd.read_csv(assemblies_file, sep='\t', dtype=str)
    has_haplotype_col = 'Haplotype' in assemblies_full.columns
    usecols = ['Sample ID', 'Population Abbreviation']
    if has_haplotype_col:
        usecols.append('Haplotype')
    assemblies = assemblies_full[usecols].copy()

    #check if the file has assemblies built into Sample ID- recommended 
    haplotype_specific = assemblies['Sample ID'].astype(str).str.contains('.', regex=False).any()
    log.info(
        f"Assemblies file {'has' if haplotype_specific else 'does not have'} "
        "haplotype-specific Sample IDs; haplotypes will be treated "
        f"{'separately' if haplotype_specific else 'together'}."
    )

    if has_haplotype_col or haplotype_specific:
        haplotypes_per_row = 1
    else:
        haplotypes_per_row = 2
        log.info(
            "Assemblies file has neither a 'Haplotype' column nor "
            "haplotype-suffixed Sample IDs; assuming each listed sample is "
            "a diploid individual contributing 2 haplotypes."
        )
    
    # assumes 'GRCh38','CHM13' if called from command line.
    exclude_samples=exclude_samples.split(',')
    exclude_samples = _expand_reference_aliases(exclude_samples)


    drop_mask = assemblies['Population Abbreviation'].astype(str).str.strip().str.lower() == 'drop'
    drop_samples = assemblies.loc[drop_mask, 'Sample ID'].tolist()
    if drop_samples:
        log.info(f"Excluding {len(drop_samples)} sample(s) flagged 'drop' in assemblies file: {drop_samples}")

    exclude_samples = _expand_reference_aliases(list(exclude_samples) + drop_samples)

    log.info(f'Filtering out the following samples: {exclude_samples}.'
            'Note: excluded assembly format should match Sample ID format in assembly file.'
             "['GRCh38','CHM13'] recommended")

    assemblies=assemblies[~assemblies['Sample ID'].isin(exclude_samples)]
    
    #dictionary of sample sizes for each population
    asm_count=Counter(assemblies.drop_duplicates()['Population Abbreviation'])
    if haplotypes_per_row != 1:
        asm_count = Counter({pop: count * haplotypes_per_row for pop, count in asm_count.items()})
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
    outf.write("\t".join(header) + "\n")


    # Create ONE persistent worker pool for the whole run, reused across every
    # region's NodeTable/LinkTable build. Avoids paying interpreter-startup
    # and import cost per region.
    pool = gbz.make_table_pool(memory_limit_gb)

    try:
        ##### If GFA, just process the whole graph #####
        if file_type == "gfa":
            if region_str is not None:
                log.warning("Regions are ignored when processing GFA")
            exclude = []
            if reference != "":
                exclude = [reference]

            try:
                node_table = gbz.build_table_with_limit(
                    gutils.NodeTable,
                    dict(gfa_file=graph_file, exclude_samples=exclude, walk_file=walk_file),
                    memory_limit_gb, log,
                    region_label="whole graph",
                    skipped_bed_file=skipped_bed_file,
                    bed_fields=(graph_file.name, "NA", "NA"),
                    pool=pool,
                )
            except gbz.RegionExtractionFailed:
                log.info("Skipping whole-graph GFA run: exceeded memory cap.")
                outf.close()
                return 0

            if 'popuniq-normdegree' in metrics_list:
                try:
                    link_table = gbz.build_table_with_limit(
                        gutils.LinkTable,
                        dict(gfa_file=graph_file, exclude_samples=exclude),
                        memory_limit_gb, log,
                        region_label="whole graph (link table)",
                        skipped_bed_file=skipped_bed_file,
                        bed_fields=(graph_file.name, "NA", "NA"),
                        pool=pool,
                    )
                except gbz.RegionExtractionFailed:
                    log.info("Skipping whole-graph GFA run: link table exceeded memory cap.")
                    outf.close()
                    return 0
                for n in node_table.nodes:
                    for l in link_table.links.keys():
                        if ((n == link_table.links[l].node_1) | (n == link_table.links[l].node_2)):
                            node_table.nodes[n].degree += 1
            else:
                link_table = None

            metric_results, matched_numwalks = compute_population_uniqueness(
                node_table, asm, asm_count, metrics_list, exclude_samples, haplotype_specific
            )

            items = [
                len(node_table.nodes.keys()),
                node_table.get_total_node_length(),
                matched_numwalks,
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
            log.info(f"Processing region {region.chrom}:{region.start}-{region.end}")
            region_label = f"{region.chrom}:{region.start}-{region.end}"

            gfa_file = gbz.extract_region_from_gbz(
                graph_file, region, reference,
                memory_limit_gb=memory_limit_gb, log=log,
                skipped_bed_file=skipped_bed_file,
            )
            if gfa_file is None:
                # already logged + written to skip bed inside extract_region_from_gbz
                continue

            try:
                if pool is not None:
                    status, summary, metric_results = pool.apply(
                        _pool_popuniq_region_worker,
                        ((gfa_file, exclude_samples, walk_file, reference, metrics_list, asm, asm_count, haplotype_specific),),
                    )
                else:
                    status, summary, metric_results = _pool_popuniq_region_worker(
                        (gfa_file, exclude_samples, walk_file, reference, metrics_list, asm, asm_count, haplotype_specific)
                    )
            except Exception as e:
                log.error(f"Unexpected error processing region {region_label}: {e}")
                with open(skipped_bed_file, "a") as skipf:
                    skipf.write(f"{region.chrom}\t{region.start}\t{region.end}\n")
                continue

            if status != "ok":
                log.warning(
                    f"Region {region_label} failed ({status}); likely exceeded {memory_limit_gb}GB cap."
                    if memory_limit_gb is not None else
                    f"Region {region_label} failed ({status})."
                )
                with open(skipped_bed_file, "a") as skipf:
                    skipf.write(f"{region.chrom}\t{region.start}\t{region.end}\n")
                continue

            log.info('computing population specific sequence uniqueness')
            items = [region.chrom, region.start, region.end] + list(summary.values()) + metric_results
            outf.write("\t".join([str(item) for item in items]) + "\n")
            outf.flush()

        ##### Cleanup #####
        end_time = time.time()
        time_per_region = (end_time - start_time) / len(regions)
        log.debug(f"Time per region\t{time_per_region}\n")
        outf.close()
        return 0

    finally:
        if pool is not None:
            pool.close()
            pool.join()

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


def compute_population_uniqueness(
    node_table: gutils.NodeTable,
    asm,
    asm_count,
    metrics: list[str],
    exclude_samples=['GRCh38','CHM13'],
    haplotype_specific: bool = False,
):
    """
    Compute population specific uniqueness for a node table for one or more
    metrics, accumulating all requested metrics in a single pass over the
    nodes. Options:
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
        Keyed by base Sample ID (no haplotype suffix) when haplotype_specific
        is False, or by full Sample ID (including haplotype suffix) when
        haplotype_specific is True.
    asm_count: dict
        dictionary of the sample sizes for each population in the total assembly
    metrics : list[str]
       Which metrics to compute, in the order they should be returned.
       See description above for valid options.
    exclude_samples: list
        List of samples to ignore for analysis (particularly those in assembly file that aren't in the assembly.)
    haplotype_specific: bool
        If True, `asm` is keyed by full walk sample name (Sample ID already
        includes a haplotype suffix, e.g. "HG00097.1"), so each haplotype is
        looked up and excluded independently rather than being collapsed to
        a base Sample ID.
    Returns
    -------
    popuniq : list[float]
       List of population uniqueness scores. For each metric (in the order
       given in `metrics`), there is 1 score per population, plus a total
       score, concatenated together in metric order. Total score calculated
       using mean Hs for all populations.
    matched_numwalks : int
       Number of distinct walk/sample names seen in `node_table` that are
       both NOT excluded (via `exclude_samples`) and present in `asm` (i.e.
       the intersection of "samples in the node table" and "samples in the
       filtered assemblies file"). Walks belonging to a sample that is
       simply absent from the assemblies file (and not itself excluded) are
       skipped rather than raising a KeyError, and are not counted here.

    Raises
    ------
    ValueError
       If invalid metric specified
    
    """
    for m in metrics:
        if m not in AVAILABLE_METRICS:
            raise ValueError(f"Invalid metric {m}")

    populations = sorted(list(asm_count.keys()))+['total']
    #TOTAL HAS TO BE LAST- it is calculated on last loop of populations as an average of the previous values.

    # Accumulators, one dict per metric, all populated together in a single
    # pass over the nodes (instead of re-looping over the nodes once per metric).
    popuniq = {m: dict.fromkeys([f'{m}_{x}' for x in populations], 0) for m in metrics}

    # Distinct walk/sample names that are both un-excluded and present in the
    # assemblies-derived population mapping. Accumulated here (rather than in
    # a separate pass) since we're already iterating every node's samples.
    matched_samples = set()

    for n in node_table.nodes.keys():
        #get list of samples present for node
        pops = []
        if haplotype_specific:
            for s in node_table.nodes[n].samples:
                if s in exclude_samples:
                    continue
                if s not in asm:
                    # Present in the graph but not in the (filtered)
                    # assemblies file - not part of the intersection.
                    continue
                pops.append(asm[s])
                matched_samples.add(s)
        else:
            for s in node_table.nodes[n].samples:
                base = s.split('.')[0]
                if (base in exclude_samples) or (s in exclude_samples):
                    continue
                if base not in asm:
                    # Present in the graph but not in the (filtered)
                    # assemblies file - not part of the intersection.
                    continue
                pops.append(asm[base])
                matched_samples.add(s)
            #list of populations present for node- take the population value from the node to pop dict

        #get dictionary of population instances
        anc_count = Counter(pops)
        anc_count = {key: anc_count.get(key, 0) for key in asm.values()}
        #count instances into dictionary

        #add attributes to class node
        node_table.nodes[n].anc_count = anc_count
        node_table.nodes[n].exp_het= calc_exp_het(asm_count, node_table.nodes[n].anc_count)
        length=node_table.nodes[n].length
        degree=node_table.nodes[n].degree
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

            fst_val = node_table.nodes[n].Fst[k]
            #calculate degree from link table for degree normalized, length otherwise
            for m in metrics:
                if m=='popuniq-normdegree':
                    popuniq[m][f'{m}_{k}']+=degree*fst_val
                else:
                    popuniq[m][f'{m}_{k}']+=length*fst_val
        ###

    n_nodes = len(node_table.nodes.keys())
    for m in metrics:
        if n_nodes>0:
            if m == 'popuniq-normwalk':
                popuniq[m] = {key: value / node_table.get_mean_walk_length() for key, value in popuniq[m].items()}
            elif m == 'popuniq-normnode':
                popuniq[m] = {key: value / node_table.get_mean_node_length() for key, value in popuniq[m].items()}
            elif m == 'popuniq-normdegree':
                mean_degree = node_table.get_mean_degree()
                popuniq[m] = {
                    key: (value / mean_degree if mean_degree != 0 else np.nan)
                    for key, value in popuniq[m].items()
                }
        else:
            popuniq[m] = {
                key: np.nan for key in popuniq[m]
            }

    results = []
    for m in metrics:
        results.extend(popuniq[m].values())
    return results, len(matched_samples)


def _pool_popuniq_region_worker(args):
    """
    Runs entirely inside the worker (or in-process if no pool): builds
    NodeTable (+ LinkTable if needed), computes population uniqueness
    metrics, cleans up the intermediate GFA file, and returns only small
    summary data — never the NodeTable/LinkTable objects themselves.
    """
    gfa_file, exclude_samples, walk_file, reference, metrics_list, asm, asm_count, haplotype_specific = args
    try:
        node_table = gutils.NodeTable(gfa_file=gfa_file, exclude_samples=exclude_samples, walk_file=walk_file)

        if 'popuniq-normdegree' in metrics_list:
            if node_table.gfa_file is not None:
                link_table = gutils.LinkTable(node_table.gfa_file, exclude_samples)
            else:
                link_table = None
            if link_table is not None:
                for n in node_table.nodes:
                    for l in link_table.links.keys():
                        if ((n == link_table.links[l].node_1) | (n == link_table.links[l].node_2)):
                            node_table.nodes[n].degree += 1

        metric_results, matched_numwalks = compute_population_uniqueness(
            node_table, asm, asm_count, metrics_list, exclude_samples, haplotype_specific
        )

        summary = {
            "numnodes": len(node_table.nodes.keys()),
            "total_length": node_table.get_total_node_length(),
            "numwalks": matched_numwalks,
        }
        return ("ok", summary, metric_results)
    except MemoryError:
        return ("memory_error", None, None)
    except Exception as e:
        return ("error", str(e), None)
    finally:
        try:
            os.remove(gfa_file)
        except OSError:
            pass