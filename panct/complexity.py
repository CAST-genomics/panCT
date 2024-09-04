# Taken from grant proposal section C.2.1 in Research Strategy
# We will define sequence uniqueness as U = ∑ s in S (|s|*p_s*(1 − p_s))/L
# S is the set of nodes in a region,
# |s| is the length in bp of node s,
# p_s is the percent of sequences that go through node s, and
# L is the average length in bp of all paths traversing the subgraph of interest.

# This metric is meant to capture the relative amount of sequence in a region that is shared vs. polymorphic amongst haplotypes in a region.

import sys
import time
import subprocess


def main():

    # total arguments
    if len(sys.argv) != 3:
        print("Usage: python calculate_complexity.py input.gfa node_map.tsv.gz")
        return

    print(f"Computing complexity of {sys.argv[1]}\t{sys.argv[2]}")

    gfa_file = sys.argv[1]
    node_map = str(sys.argv[2])

    complexity = complexity_score(gfa_file, node_map)

    print(f"Complexity\t{complexity}\n")


def complexity_score(gfa_file: str, node_map: str):
    """
    Compute a complexity score for a given GFA file

    Parameters
    ----------
    gfa_file : str
        The path to the GFA file
    node_map : str
        The path to the node map file

    Returns
    -------
    float
        The complexity score of the GFA file
    """
    start_time = time.time()

    gfa = open(gfa_file, "r")

    complexity = 0
    total_length = 0
    number_of_nodes = 0

    # TODO should this be over the whole gfa or just in the subsetted region?
    # 90 is toatl number of haplotypes in minigraph cactus
    total_haplotypes = 90

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
            addition = length_int * p_s * (1 - p_s)

            # print("addition:", addition)
            complexity += addition

    average_length = total_length / number_of_nodes
    complexity = complexity / (average_length)

    end_time = time.time()
    time_per_node = (end_time - start_time) / number_of_nodes
    print(f"Time per node\t{time_per_node}")

    return complexity


def get_p_s(target_node: str, node_map: str, total_haplotypes: int) -> float:
    """
    Calculate the percent of haplotypes that travel through the given node

    Parameters
    ----------
    target_node : str
        The node to calculate the percent of haplotypes that travel through
    node_map : str
        The path to the node map file
    total_haplotypes : int
        The total number of haplotypes in the dataset

    Returns
    -------
    float
        The percent of haplotypes that travel through the given node
    """
    command = ["tabix", node_map, f":{target_node}-{target_node}"]
    tabix_output = (
        subprocess.run(command, stdout=subprocess.PIPE).stdout.decode("utf-8").strip()
    )
    haplotypes = tabix_output.split("\t")[2:]  # list of all haplotypes for node
    return len(haplotypes) / total_haplotypes
