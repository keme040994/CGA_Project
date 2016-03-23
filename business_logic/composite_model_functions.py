# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: March, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *


# =========================
# COMPOSITE MODEL FUNCTIONS
# =========================

# FUNCTION: composite_model_creator
def composite_model_creator(population, cant_genes):
    """
    Creates the composite model by multiplying every chromosome to its fitness and summing them all in a new matrix.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        cant_genes : INT
            The amount of genes for each chromosome

    Returns:
        MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The composite model
    """
    new_population = fitness_x_genes(population)
    composite_model = [[0 for i in range(cant_genes)] for i in range(cant_genes)]
    for i in range(0, len(new_population)):
        genes = population[i].get_genes()
        for j in range(0, len(genes)):
            for k in range(0, len(genes)):
                composite_model[j][k] += genes[j][k]
    return round_matrix(composite_model, 3)


# FUNCTION: fitness_x_genes
def fitness_x_genes(population):
    """
    This is an auxiliary function that multiplies the fitness by the genes on every chromosome of a given population.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        LIST[Chromosome(), Chromosome(), ...]
            The new population with the genes multiplied by the fitness
    """
    new_population = []
    for i in range(0, len(population)):
        genes = population[i].get_genes()
        fitness = population[i].get_fitness()
        for j in range(0, len(genes)):
            for k in range(0, len(genes)):
                genes[j][k] *= fitness
        new_population.append(Chromosome(genes))
    return new_population
