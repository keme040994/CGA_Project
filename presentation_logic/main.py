# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: March, 2016.

# LIBRARIES
from business_logic.composite_model_functions import *
from business_logic.initialization_functions import *
from business_logic.mutation_functions import *
from business_logic.repair_functions import *
from business_logic.selection_functions import *
from copy import deepcopy
import data_logic.data_input as data
from presentation_logic.output_functions import *


# ========
#   MAIN
# ========

# FUNCTION: main
def main():
    """ The main function program, having the structure of the CGA. """

    # Initial Variables
    pop_size = 250  # Population size
    cant_genes = 12  # Amount of genes for each chromosome
    per_ones = 5  # Percentage of ones based on matrix's size
    likelihood_function = 1  # 1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two
    per_mutations = 35  # Percentage of mutations
    per_elitism = 10  # Percentage of elitism
    selection_prop = 0.05  # Probability of selection
    mutation_prop = 0.3  # Probability of mutation
    cant_mutations = 1  # Amount of mutations per chromosome
    per_filter_cm = 0.7  # For filtering values in the composite model
    per_filter_am = 0.7  # For filtering values in the amalgamated model
    cant_matings = 500  # Amount of matings (loops)
    cant_composite_model = 3  # Amount of composite models fro creation

    current_population = seed_population(pop_size, cant_genes, per_ones, likelihood_function, data.rep)
    threshold = ((cant_genes*cant_genes)-cant_genes)//4

    for i in range(0, cant_matings):
        print("Mating #" + str(i+1))

        # Selection
        children_population = repair_population(selection_function(current_population, threshold))

        if len(children_population) == 0:
            threshold -= 1
        else:
            full_population = current_population + children_population
            current_population = select_best_chromosomes(full_population, pop_size, likelihood_function, data.rep)

        if threshold < 0:
            cant_mutations = (cant_genes*cant_genes*per_mutations)//100
            current_population = cataclysmic_mutation(select_the_best_chromosome(current_population,
                                                                                 likelihood_function, data.rep),
                                                      pop_size, cant_genes, cant_mutations)
            current_population = repair_population(current_population)
            threshold = ((cant_genes*cant_genes)-cant_genes)//4

    # Removing from the population the repeated chromosomes
    likelihood_result_calculator(current_population, likelihood_function, data.rep)
    relative_likelihood_result_calculator(current_population)
    relative_likelihood_result_sorting(current_population)
    ranked_selection_calculator(current_population)
    current_population = select_uniques_chromosomes(current_population)

    composite_model = composite_model_creator(current_population, cant_genes)

    if likelihood_function == 1:
        composite_model = sum_matrix(composite_model, transpose_matrix(composite_model))

    view_model(composite_model, "COMPOSITE MODEL (ROUNDED TO 3 DECIMALS)")

    if likelihood_function == 1:
        create_model_image_cotemporal(convert_model_to_digraph(composite_model, per_filter_cm, data.a_protein_names),
                                      data.a_name + " - (Composite Mode)")
    else:
        create_model_image(convert_model_to_digraph(composite_model, per_filter_cm, data.a_protein_names),
                           data.a_name + " - (Composite Mode)")

    write_file(data.a_name + " - (Composite Mode)", str(composite_model))
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

if __name__ == "__main__":
    main()
