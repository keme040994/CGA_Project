# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: March, 2016.

# LIBRARIES
from business_logic.chromosome import Chromosome
from copy import deepcopy
from random import randint
import random


# ==================
# MUTATION FUNCTIONS
# ==================

# FUNCTION: cataclysmic mutation
def cataclysmic_mutation(best_chromosome, pop_size, cant_genes, cant_mutations):
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
                    random1 = randint(0, len(population[0].get_genes())-1)
                    random2 = randint(0, len(population[0].get_genes())-1)
                population[i].bit_changer(random1, random2)
