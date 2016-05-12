# Created by: Dr. David John & Kenneth Meza.
# Created at: March, 2016.
# Updated at: May, 2016.

# LIBRARIES
import bigfloat as bf
from business_logic.log_pnml_functions import *
import networkx as nx
import numpy as np
import sys


# =================
# GENERAL FUNCTIONS
# =================

# FUNCTION: get_positions
def get_positions(matrix_a, matrix_b):
    """
    Obtains the positions where two matrices differ in a list of tuples.

    Args:
        matrix_a : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

        matrix_b : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

    Returns:
        LIST[TUPLE, TUPLE]
            A list of tuple with the following format: (x, y)
    """
    positions = []
    for i in range(0, len(matrix_a)):
        for j in range(0, len(matrix_a)):
            if matrix_a[i][j] != matrix_b[i][j] and i != j:  # the condition i != j avoids diagonal position
                positions.append((i, j))
    return positions


# FUNCTION: calculate_hamming_distance
def calculate_hamming_distance(matrix_a, matrix_b):
    """
    Calculates the hamming distance between two matrices.

    Args:
        matrix_a : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

        matrix_b : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

    Returns:
        INT
            The hamming distance value
    """
    return len(get_positions(matrix_a, matrix_b))


# FUNCTION: select_best_chromosomes
def select_best_chromosomes(population, pop_size, likelihood_function, rep):
    """
    Given a population, this function returns the best chromosomes.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        pop_size : INT
            The size of the population
        likelihood_function : INT
            Indicates the likelihood function that is going to be used
                1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two
        rep : LIST[rep1, rep2, rep3, ...]
            A repN is a biological data used to calc the likelihood result

    Returns:
        LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects, containing the new population
    """
    likelihood_result_calculator(population, likelihood_function, rep)
    relative_likelihood_result_calculator(population)
    relative_likelihood_result_sorting(population)
    return population[:pop_size]


# FUNCTION: select_the_best_chromosome
def select_the_best_chromosome(population, likelihood_function, rep):
    """
    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        likelihood_function : INT
            Indicates the likelihood function that is going to be used
                1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two
        rep : LIST[rep1, rep2, rep3, ...]
            A repN is a biological data used to calc the likelihood result

    Returns:
        Chromosome()
            The best chromosome in the population
    """
    likelihood_result_calculator(population, likelihood_function, rep)
    relative_likelihood_result_calculator(population)
    relative_likelihood_result_sorting(population)
    return population[0]


# ======================================================================================================================


# FUNCTION: filter_likelihood_selection
def filter_likelihood_selection(num):
    """
    Stops the execution if the likelihood function type is not valid.

    Args:
        num : INT
            A number representing the desired likelihood function to use:
                1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two, n = not valid type
    """
    if num == 1:
        print("* LIKELIHOOD FUNCTION: Cotemporal")
    elif num == 2:
        print("* LIKELIHOOD FUNCTION: Next-Step One")
    elif num == 3:
        print("* LIKELIHOOD FUNCTION: Next-Step One-Two")
    else:
        sys.exit("ERROR: You need to select a valid likelihood function type.")


# FUNCTION: view_population
def view_population(population):
    """
    Shows the genes of every chromosome in a given population.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
    """
    for i in range(0, len(population)):
        print(str(population[i]) + "\n")


# FUNCTION: likelihood_result_calculator
def likelihood_result_calculator(population, likelihood_function, rep):
    """
    Given a population and the replicates, it adds the likelihood result to every chromosome.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
        likelihood_function: INT
            Indicates the likelihood function that is going to be used
                1 = cotemporal, 2 = next_step_one, 3 = next_step_one_two
        rep : LIST[rep1, rep2, rep3, ...]
            A repN is a biological data used to calc the likelihood result
    """
    if likelihood_function == 1:
        for i in range(0, len(population)):
            di_graph = convert_matrix_to_digraph(population[i].get_genes())
            population[i].set_log_likelihood_result(log_pnml_cotemporal(di_graph, rep))
    elif likelihood_function == 2:
        for i in range(0, len(population)):
            di_graph = convert_matrix_to_digraph(population[i].get_genes())
            population[i].set_log_likelihood_result(log_pnml_next_step_one(di_graph, rep))
    elif likelihood_function == 3:
        for i in range(0, len(population)):
            di_graph = convert_matrix_to_digraph(population[i].get_genes())
            population[i].set_log_likelihood_result(log_pnml_next_step_one_two(di_graph, rep))


# FUNCTION: relative_likelihood_result_sorting
def relative_likelihood_result_sorting(population):
    """
    Sorts a population given, in descending order, based on the relative likelihood result of every chromosome.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
    """
    population.sort(key=lambda x: x.relative_likelihood_result, reverse=True)


# FUNCTION: relative_likelihood_result_calculator
def relative_likelihood_result_calculator(population):
    """
    Given a population, this function calculates the relative likelihood result (it is a way to make
    the likelihood result bigger) to every chromosome.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects
    """
    with bf.quadruple_precision:
        total = sum_likelihood_result(population)
        for i in range(0, len(population)):
            log_likelihood_result = bf.exp(bf.BigFloat(str(population[i].get_log_likelihood_result())))
            population[i].set_relative_likelihood_result(float(bf.div(log_likelihood_result, total)))


# FUNCTION: sum_likelihood_result
def sum_likelihood_result(population):
    """
    Sums all the likelihood results values on a population.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        INT
            The sum of all likelihood results on a population
    """
    with bf.quadruple_precision:
        total = bf.BigFloat("0.0")
        for i in range(0, len(population)):
            log_likelihood_result = bf.exp(bf.BigFloat(str(population[i].get_log_likelihood_result())))
            total = bf.add(total, log_likelihood_result)
    return total


# FUNCTION: convert_matrix_to_digraph
def convert_matrix_to_digraph(matrix):
    """
    Converts a python's matrix into a networkx's digraph.

    Args:
        matrix : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

    Returns:
        nx.DiGraph()
            The converted matrix into a digraph
    """
    genes = np.asmatrix(matrix)
    di_graph = nx.to_networkx_graph(genes, create_using=nx.DiGraph())
    return di_graph


# FUNCTION: convert_digraph_to_matrix
def convert_digraph_to_matrix(digraph):
    """
    Converts a networkx's digraph into a python's matrix.

    Args:
        digraph : nx.DiGraph()
            It's a networkx's digraph representing genes

    Returns:
        MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            The converted diagraph into a mtrix
    """
    return [[int(i) for i in j] for j in nx.to_numpy_matrix(digraph).tolist()]


# FUNCTION: transpose_matrix
def transpose_matrix(matrix):
    """
    Given a matrix, this function returns its transpose.

    Args:
        matrix : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

    Returns:
        MATRIX[[INT, INT, ...], [INT, INT, ...], ...]

    """
    t_matrix = []
    for i in range(0, len(matrix)):
        t_matrix.append(extract_column(matrix, i))
    return t_matrix


# FUNCTION: extract_column
def extract_column(matrix, column_number):
    """
    An auxiliary function that extracts a column given a matrix.

    Args:
        matrix : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes
        column_number : INT
            A number between 0 to N-1

    Returns:
        LIST[INT, INT, INT, ...]
            A list that contains the elements of the column
    """
    column = []
    for i in range(0, len(matrix)):
        column.append(matrix[i][column_number])
    return column


# FUNCTION: select_uniques_chromosomes
def select_uniques_chromosomes(population):
    """
    Returns a filtered population with only unique elements. The comparison criteria is the 'relative likelihood
    result'.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        LIST[Chromosome(), Chromosome(), ...]
            A new list filled with 'Chromosome' objects
    """
    unique_population = []
    for i in range(0, len(population)):
        if find_element_on_list(population[i].get_relative_likelihood_result(), unique_population):
            unique_population.append(population[i])
    return unique_population


# FUNCTION: find_element_on_list
def find_element_on_list(relative_likelihood_result, unique_population):
    """
    An auxiliary function that returns a boolean value if the element exists on the 'unique_population' parameter.

    Args:
        relative_likelihood_result : FLOAT
            A number representing the value of certain chromosome's relative likelihood result
        unique_population: LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects where each one is unique

    Returns:
        BOOLEAN
            True it the value is found or False if not
    """
    for i in range(0, len(unique_population)):
        # based on the 'relative likelihood  result'
        if relative_likelihood_result == unique_population[i].get_relative_likelihood_result():
            return False
    return True


# FUNCTION: view_model
def view_model(model, model_title):
    """
    Displays a model (composite or amalgamated) and all its content

    Args:
        model : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given model
        model_title : STRING
            A title to be putted on the displaying
    """
    print(model_title)
    for row in model:
        for val in row:
            print("{:10.3f}".format(val), end="")
        print()


# FUNCTION: convert_model_to_digraph
def convert_model_to_digraph(model, per_max_model, a_protein_names):
    """
    A special conversion function where given a model (composite or amalgamated) and the proteins names it combines all
    in a new digraph that will contains the information necessary for displaying the graph.

    Args:
        model : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given model
        per_max_model : FLOAT
            A number used for filtering purpose
        a_protein_names : LIST[STRING, STRING, STRING, ...]
            The name of all the proteins given in the variable 'a_protein_names' on 'data_logic.data_input.py'

    Returns:
        nx.DiGraph()
            The converted model into the true digraph format for displaying
    """
    cm_digraph = convert_matrix_to_digraph(model)
    di_graph = nx.DiGraph()
    di_graph.add_nodes_from(a_protein_names)
    for i in range(0, len(cm_digraph.edges())):
        edge = cm_digraph.edges()[i]
        if cm_digraph[edge[0]][edge[1]]['weight'] >= per_max_model:
            di_graph.add_edge(a_protein_names[edge[0]], a_protein_names[edge[1]],
                              weight=cm_digraph[edge[0]][edge[1]]['weight'])
    return di_graph


# FUNCTION: fitness_average
def fitness_average(population):
    """
    Given a population, it calculates the average of the fitness value.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        FLOAT
            The average of the fitness value in the population
    """
    fitness_avg = 0
    for item in population:
        fitness_avg += item.get_fitness()
    return fitness_avg/len(population)


# FUNCTION: fitness_variance
def fitness_variance(population):
    """
    Given a population, it calculates the variance of the fitness value.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        FLOAT
            The variance of the fitness value in the population
    """
    fitness_avg = fitness_average(population)
    fitness_var = 0
    for i in range(0, len(population)):
        fitness_var += (population[i].get_fitness()-fitness_avg)**2
    return fitness_var/len(population)


# FUNCTION: rlr_average
def rlr_average(population):
    """
    Given a population, it calculates the average of the relative likelihood result.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        FLOAT
            The average of the relative likelihood result in the population
    """
    rlr_avg = 0
    for item in population:
        rlr_avg += item.get_relative_likelihood_result()
    return rlr_avg/len(population)


# FUNCTION: rlr_variance
def rlr_variance(population):
    """
    Given a population, it calculates the variance of the relative likelihood result.

    Args:
        population : LIST[Chromosome(), Chromosome(), ...]
            A list filled with 'Chromosome' objects

    Returns:
        FLOAT
            The variance of the relative likelihood result in the population
    """
    rlr_avg = rlr_average(population)
    rlr_var = 0
    for i in range(0, len(population)):
        rlr_var += (population[i].get_relative_likelihood_result()-rlr_avg)**2
    return rlr_var/len(population)


# FUNCTION: sum_matrix
def sum_matrix(matrix_a, matrix_b):
    """
    Sums two matrix.

    Args:
        matrix_a : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A matrix to be summed
        matrix_b : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            Another matrix to be summed

    Returns:
        MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            The summed matrix in a new matrix
    """
    for i in range(0, len(matrix_a)):
        for j in range(0, len(matrix_a)):
            matrix_a[i][j] += matrix_b[i][j]
    return matrix_a


# FUNCTION: round_matrix
def round_matrix(matrix, cant_decimals):
    """
    Rounds out every element on a matrix, given an amount of decimals.

    Args:
        matrix : MATRIX[[FLOAT, FLOAT, ...], [FLOAT, FLOAT, ...], ...]
            A given matrix
        cant_decimals : INT
            The amount of decimals desired for the rounding

    Returns:
        A rounded matrix
    """
    for i in range(0, len(matrix)):
        for j in range(0, len(matrix)):
            matrix[i][j] = round(matrix[i][j], cant_decimals)
    return matrix
