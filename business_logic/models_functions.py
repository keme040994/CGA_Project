# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: May, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *


# ===============
# MODEL FUNCTIONS
# ===============

# FUNCTION: model_creator
def model_creator(population, num_genes, likelihood_function):
    """
    Creates the composite / amalgamated model by multiplying every chromosome to its fitness and summing them all in
    a new matrix.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        num_genes : INT
            The amount of genes for each chromosome
        likelihood_function: INT
            Indicates the likelihood function that is going to be used
                1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two

    Returns:
        MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The composite / amalgamated model
    """
    new_population = fitness_x_genes(population)
    model = [[0 for i in range(num_genes)] for i in range(num_genes)]
    for i in range(0, len(new_population)):
        genes = population[i].get_genes()
        for j in range(0, len(genes)):
            for k in range(0, len(genes)):
                model[j][k] += genes[j][k]

    # all values in the model are rounded by 3 decimals
    model = round_matrix(model, 3)

    if likelihood_function == 1:
        sum_matrix(model, transpose_matrix(model))

    return model


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
