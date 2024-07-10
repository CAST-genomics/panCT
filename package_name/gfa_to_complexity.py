# Taken from grant proposal section C.2.1 in Research Strategy
# We will define sequence uniqueness as U = ∑ s in S (|s|*p_s*(1 − p_s))/L 
# S is the set of nodes in a region,
#|s| is the length in bp of node s, 
# p_s is the percent of sequences that go through node s, and 
# L is the average length in bp of all segments traversing the subgraph of interest

# In type 1 we define L as the average length of all segments
# In type 2 we define L as the length of the linear genome used to segment the gfa

# This metric is meant to capture the relative amount of sequence in a region 
# that is shared vs. polymorphic amongst haplotypes in a region.

# TODO: Create easily runnable unit tests for this metric
# TODO: Generate method to calculate average path length from start -> finish.
#    Requires ordering of nodes.

'''
Example Usage:
    python gfa_to_complexity.py [gfa_file] [output_file]
'''

import subprocess
import sys
import os

def main():
    gfa, outp, gfa_input = read_inputs()
    node_map = "/home/wwford/datasets/hprc-v1.1-mc-grch38-node-sample-map-tabix.tsv.gz"

    # Score 1 Function
    complexity = score1(gfa, node_map)

    write_score(outp, gfa_input, complexity)

def read_inputs():
    # Verify inputs
    gfa_input = sys.argv[1]
    out_file = sys.argv[2]

    if not os.path.isfile(gfa_input):
        raise ValueError("Gfa input path is broken!")

    if os.path.isfile(out_file):
        outp = open(out_file, 'w+')
    else:
        outp = open(out_file, 'w')

    gfa = open(gfa_input,'r')

    return gfa, outp, gfa_input

def write_score(outp, gfa_file, complexity):
    coords = "\t".join(gfa_file.split("/")[-1].split(".")[0].rsplit("_", maxsplit=2))
    outp.write(f"{coords}\t{complexity}\n")

def score1(gfa: str, node_map: str, total_haplotypes:int = 94) -> float:
    """
    Score a gfa file, normalizing by the average bp length of each node.

    Args:
        gfa_file (str) : file path
        node_map (str) : node map from corresponding pangenome
        total_haplotypes (int) : default = 90, number of haps in v1.1 of the Human Pangenome

    Returns:
        float : the 'complexity score' of a gfa file
    """

    complexity = 0
    total_length = 0
    number_of_nodes = 0

    for line in gfa.readlines():
        if line.startswith("S"):
            # Get the sequence length and node id
            vars = line.split(sep="\t")
            node = vars[1]
            length = False
            for var in vars[3:]:
                if var.startswith("LN"):
                    length = var
                    break
            
            if not length: 
                print(f"Error! Node {node} has no length")
                return

            length_int = int(length.split(sep=":")[2])

            # Compute Average Segment Length
            total_length += length_int
            number_of_nodes += 1

            # Calculate p_s
            p_s = get_p_s(node, node_map, total_haplotypes)

            # Add complexity
            addition = length_int*p_s*(1-p_s)

            complexity += addition

    average_length = total_length/number_of_nodes
    complexity = complexity/(average_length)

    return complexity

def get_p_s(target_node:str, node_map:str, total_haplotypes:int) -> float:
    """
    Calcualtes the percent of haplotypes that travel through the given node
    Args:
        target_node (str)
        node_map (str) : file path to tabix compressed node map
        total_haplotypes (int) : The total number of haplotypes 
    Returns:
        float: the percentage of total haplotypes that pass through this node
    """
    command = ["tabix", node_map, f":{target_node}-{target_node}"]
    tabix_output = subprocess.run(command, stdout=subprocess.PIPE).stdout.decode('utf-8').strip()
    haplotypes = tabix_output.split('\t')[2:] # list of all haplotypes for node
    return len(haplotypes)/total_haplotypes

def snakemake_call(gfa_input, out_file):

    if not os.path.isfile(gfa_input):
        raise ValueError("Gfa input path is broken!")

    if os.path.isfile(out_file):
        outp = open(out_file, 'w+')
    else:
        outp = open(out_file, 'w')

    gfa = open(gfa_input,'r')

    node_map = "/home/wwford/datasets/hprc-v1.1-mc-grch38-node-sample-map-tabix.tsv.gz"

    complexity = score1(gfa, node_map)

    write_score(outp, gfa_input, complexity)

if __name__ == "__main__":
    try:
        import snakemake
        # Check if the script is executed within a Snakemake rule
        if snakemake:
            # Call my_function with input and output files from Snakemake
            snakemake_call(snakemake.input[0], snakemake.output[0])
            return
    except ImportError:
        # Handle the non snakemake call
        main()