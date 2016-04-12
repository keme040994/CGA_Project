# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: April, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *
from copy import deepcopy
from random import randint


# ===================
# SELECTION FUNCTIONS
# ===================

# FUNCTION: fitness_calculator
def fitness_calculator(ordered_population):
    """
    Calculates the 'fitness' that will be used for doing a rank based selection.

    Args:
        ordered_population : LIST[Chromosome]
            A list filled with 'Chromosome' objects sorted by the likelihood result
    """
    size = len(ordered_population)
    for i in range(0, size):
        ordered_population[i].fitness = (2*(size+1-(i+1)))/(size*(size+1))


# FUNCTION: selection_function
def selection_function(population, threshold):
    """
    Creates a new generation of chromosomes, based on the CHC Algorithm.

    Args:
        population : LIST[Chromosome]
            A list filled with 'Chromosome' objects
        threshold : INT
            The limit used as a filter in the algorithm

    Returns:
        LIST[Chromosome]
            A list filled with 'Chromosome' objects, containing the new population
    """
    population_visited = [i for i in range(0, len(population))]
    children_population = []
    for i in range(1, len(population)//2):
        random1 = randint(0, len(population_visited)-1)
        random2 = randint(0, len(population_visited)-1)
        while random1 == random2:
            random2 = randint(0, len(population_visited)-1)

        index_a = population_visited[random1]
        index_b = population_visited[random2]

        parent_a = population[index_a]
        parent_b = population[index_b]

        hamming_distance = calculate_hamming_distance(parent_a.get_genes(), parent_b.get_genes())
        if hamming_distance > threshold:
            children = HUX_crossover_function(parent_a, parent_b)
            children_population.append(children[0])
            children_population.append(children[1])
            population_visited.remove(index_a)
            population_visited.remove(index_b)
    return children_population


# FUNCTION: HUX_crossover_function
def HUX_crossover_function(parent_a, parent_b):
    """
    The HUX crossover (Half Uniform Crossover) consists in finding all the bits that differ between two parents and
    swapping the half of this bits, creating two new offsprings.

    Args:
        parent_a : Chromosome
            An object 'Chromosome' that will be used for mating
        parent_b : Chromosome
            An object 'Chromosome' that will be used for mating

    Returns:
        TUPLE(Chromosome, Chromosome)
            A tuple formed by two offsprings, represented by a 'Chromosome' object
    """
    genes_a = deepcopy(parent_a.get_genes())
    genes_b = deepcopy(parent_b.get_genes())

    positions = get_positions(genes_a, genes_b)
    for i in range(0, len(positions)//2):
        index_pos = randint(0, len(positions)-1)
        x = positions[index_pos][0]
        y = positions[index_pos][1]
        genes_a[x][y], genes_b[x][y] = genes_b[x][y], genes_a[x][y]
        positions.remove(positions[index_pos])
    return Chromosome(genes_a), Chromosome(genes_b)
