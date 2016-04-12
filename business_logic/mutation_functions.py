# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: April, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from copy import deepcopy
from random import randint
import random


# ==================
# MUTATION FUNCTIONS
# ==================

# FUNCTION: cataclysmic_mutation
def cataclysmic_mutation(best_chromosome, pop_size, cant_genes, cant_mutations):
    """
    This function consists on keeping one copy of the best element in the population unchanged, and filling the rest of
    the population by taking that best individual, creating a copy of it, and then flip a random 35% of its bits until
    having N individuals again.

    Args:
        best_chromosome : Chromosome()
            The best chromosome in the population
        pop_size : INT
            The size of the population
        cant_genes : INT
            The amount of genes for each chromosome
        cant_mutations : INT
            The desired amount of mutations per chromosome

    Returns:
         LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects, containing the new population
    """
    new_population = [best_chromosome]
    for i in range(1, pop_size):
        new_chromosome = Chromosome(deepcopy(best_chromosome.get_genes()))
        for j in range(0, cant_mutations):
            random1 = randint(0, cant_genes-1)
            random2 = randint(0, cant_genes-1)
            while random1 == random2:  # for avoiding having two equal random numbers
                random2 = randint(0, cant_genes-1)
            new_chromosome.bit_changer(random1, random2)
        new_population.append(new_chromosome)
    return new_population


# FUNCTION: mutation_function
def mutation_function(population, mutation_prob, cant_mutations):
    """
    Given a probability of mutation, this function applies a mutation (change of a bit in the genes) to every
    chromosome in the population.

    Args:
        population : LIST[Chromosome]
            A list filled with 'Chromosome' objects
        mutation_prob : FLOAT
            The desired probability of mutation
        cant_mutations : INT
            The desired amount of mutations per chromosome
    """
    for i in range(0, len(population)):
        for j in range(0, cant_mutations):
            if random.random() <= mutation_prob:
                random1 = randint(0, len(population[0].get_genes())-1)
                random2 = randint(0, len(population[0].get_genes())-1)
                while random1 == random2:  # for avoiding having two equal random numbers
                    random2 = randint(0, len(population[0].get_genes())-1)
                population[i].bit_changer(random1, random2)
