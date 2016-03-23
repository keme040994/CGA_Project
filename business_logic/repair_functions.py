# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: March, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *
from random import randint
import networkx as nx


# ================
# REPAIR FUNCTIONS
# ================

# FUNCTION: repair_population
def repair_population(population):
    """
    Repairs a DAG when it has cycles to be deleted. This function runs 3 important steps:
        I) Eliminates cycles on the diagonal of the matrix
        II) Eliminates ones that appear more than 3 times on each column of the matrix
        III) Eliminates the edges that appears more in the cycles and breaks simple cycles

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        LIST[Chromosome()]
            A new list of chromosomes with fixed DAG (no cycles)
    """
    repaired_population = []
    for i in range(0, len(population)):
        # I)
        diagonal_loops_breaker(population[i].get_genes())

        # II)
        population[i].set_genes(transpose_matrix(population[i].get_genes()))
        for j in range(0, len(population[i].get_genes())):
            while sum(population[i].get_genes()[j]) > 3:
                list_positions_ones = [n for n, x in enumerate(population[i].get_genes()[j]) if x == 1]
                random_one = randint(0, len(list_positions_ones)-1)
                population[i].get_genes()[j][list_positions_ones[random_one]] = 0
        population[i].set_genes(transpose_matrix(population[i].get_genes()))

        # III)
        di_graph = convert_matrix_to_digraph(population[i].get_genes())
        while not find_most_repeated_cycles(di_graph) == []:
            most_repeated_cycles = find_most_repeated_cycles(di_graph)
            most_repeated_cycles.sort(key=lambda x: x[1], reverse=True)  # A sorting of the list, bigger ones first
            di_graph.remove_edge(most_repeated_cycles[0][0][0], most_repeated_cycles[0][0][1])

        repaired_population.append(Chromosome(convert_digraph_to_matrix(di_graph)))
    return repaired_population


# FUNCTION: diagonal_loops_breaker
def diagonal_loops_breaker(matrix):
    """
    Breaks cycles when the i and j are the same. Example: edge (1,1), edge (2,2), etc.

    Args:
        matrix : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix represents the genes of a chromosome
    """
    for i in range(0, len(matrix)):
        matrix[i][i] = 0


# FUNCTION: find_most_repeated_cycles
def find_most_repeated_cycles(di_graph):
    """
    Returns a list filled with this format for each element: [edge : amount_of_appearances].

    Args:
        di_graph : nx.DiGraph()
            A networkx DiGraph class for representing DAG

    Returns:
        MATRIX[[TUPLE, INT], [TUPLE, INT], [TUPLE, INT], ...]
            If we have at least one edge with one appearance
        MATRIX[]
            If we don't have edges
    """
    list_all_cycles = []
    cycles = list(nx.simple_cycles(di_graph))
    for i in range(0, len(cycles)):
        list_all_cycles.append(find_cycle_edges(cycles[i], di_graph.edges(cycles[i])))
    flatted_edges = sum(list_all_cycles, [])  # This flattens the nested list of edges

    # This list contains a list of edges and their appearances on the list, but only appearances bigger than 0
    checked_edges = []
    while len(flatted_edges) > 0:
        cont = flatted_edges.count(flatted_edges[0])
        if cont > 0:  # Amount of appearances bigger than 1
            checked_edges.append([flatted_edges[0], cont])
        # This remove a value from a list
        flatted_edges[:] = (value for value in flatted_edges if value != flatted_edges[0])
    return checked_edges


# FUNCTION: find_cycle_edges
def find_cycle_edges(nodes, all_edges):
    """
    Returns all the edges between a set of certain nodes.

    Args:
        nodes : LIST[INT, INT, ...]
            The nodes that matters for this filtering edges purpose
        all_edges : LIST[TUPLE, TUPLE, ...]
            All the edges of certain chromosome

    Returns:
        LIST[TUPLE, TUPLE, ...]
            This list contains all the edges that belong to a set of nodes
    """
    edges = []
    for i in range(0, len(all_edges)):
        edge = all_edges[i]

        cont = 0
        for j in range(0, len(nodes)):
            if edge[0] == nodes[j]:
                cont += 1
            if edge[1] == nodes[j]:
                cont += 1

        if cont == 2:
            edges.append(edge)
    return edges
