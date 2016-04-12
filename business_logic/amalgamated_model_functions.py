# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: April, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from business_logic.general_functions import *


# ===========================
# AMALGAMATED MODEL FUNCTIONS
# ===========================

# FUNCTION: amalgamated_model_creator
def amalgamated_model_creator(amalgamated_population, cant_genes):
    """
    Creates the amalgamated model by having an amalgamated population (every last chromosome of every last generation
    of composite model creation), multiplying every chromosome to its relative likelihood result, summing them all in
    a new matrix and dividing every element in the new matrix by the summed relative likelihood result.
    Works for NextStep One and NextStep One-Two paradigms.

    Args:
        amalgamated_population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        cant_genes : INT
            The amount of genes for each chromosome

    Returns:
        MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The amalgamated model
    """
    amalgamated_model = [[0 for i in range(cant_genes)] for i in range(cant_genes)]
    amalgamated_population = relative_likelihood_result_x_genes(amalgamated_population)

    # Calc of the summed relative likelihood results
    summed_rlr = 0.0
    for i in range(0, len(amalgamated_population)):
        summed_rlr += amalgamated_population[i].get_relative_likelihood_result()

    for chromosome in amalgamated_population:
        genes = chromosome.get_genes()
        for i in range(0, cant_genes):
            for j in range(0, cant_genes):
                amalgamated_model[i][j] += genes[i][j]
    return round_matrix(element_divisor(amalgamated_model, summed_rlr), 3)


# FUNCTION: amalgamated_model_creator_cotemporal
def amalgamated_model_creator_cotemporal(amalgamated_population, cant_genes):
    """
    Creates the amalgamated model by having an amalgamated population (every last chromosome of every last generation
    of composite model creation), multiplying every chromosome to its relative likelihood result, summing them all in
    a new matrix and dividing every element in the new matrix by the summed relative likelihood result.
    Works for Cotemporal paradigm.

    Args:
        amalgamated_population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        cant_genes : INT
            The amount of genes for each chromosome

    Returns:
        MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The amalgamated model
    """
    amalgamated_model = [[0 for i in range(cant_genes)] for i in range(cant_genes)]
    amalgamated_population = relative_likelihood_result_x_genes(amalgamated_population)

    # Calc of the summed relative likelihood results
    summed_rlr = 0.0
    for i in range(0, len(amalgamated_population)):
        summed_rlr += amalgamated_population[i].get_relative_likelihood_result()

    for chromosome in amalgamated_population:
        genes = chromosome.get_genes()
        for i in range(0, cant_genes):
            for j in range(0, cant_genes):
                amalgamated_model[i][j] += genes[i][j] + genes[j][i]
    return round_matrix(element_divisor(amalgamated_model, summed_rlr), 3)


# FUNCTION: relative_likelihood_result_x_genes
def relative_likelihood_result_x_genes(population):
    """
    This is an auxiliary function that multiplies the relative likelihood result by the genes on every chromosome of a
    given population.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        LIST[Chromosome(), Chromosome(), ...]
            The new population with the genes multiplied by the relative likelihood result
    """
    new_population = []
    for i in range(0, len(population)):
        genes = population[i].get_genes()
        relative_likelihood_result = population[i].get_relative_likelihood_result()
        for j in range(0, len(genes)):
            for k in range(0, len(genes)):
                genes[j][k] *= relative_likelihood_result
        new_population.append(Chromosome(genes))
        new_population[-1].set_relative_likelihood_result(relative_likelihood_result)
    return new_population


# FUNCTION: element_divisor
def element_divisor(matrix, value):
    """
    Divides every element on a matrix by a given value.

    Args:
        matrix : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given matrix
        value : INT
            A given number

    Returns:
        matrix : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The matrix divided by the given value
    """
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            matrix[i][j] /= value
    return matrix
