# Created by: Dr. David John & Kenneth Meza.
# Created at: January, 2016.
# Updated at: May, 2016.


# CLASS: Chromosome
class Chromosome:
    """
    This class contains the structure that is necessary for every individual in a population. It also has some auxiliary
    methods used by other external functions.

    Args:
        genes : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            This matrix is the representation of the genes

    Attributes:
        genes : MATRIX[[INT, INT, ...], [INT, INT, ...], ...]
            The matrix associated with the class, representing the genes
        likelihood_result : FLOAT
            A value calculated by the likelihood function selected by the user, very important to measure the quality of
            the genes
        relative_likelihood_result : FLOAT
            An important value, calculated based on the likelihood result by external functions
        fitness : FLOAT
            This final value is calculated for selection purpose based on the population sorted by the relative
            likelihood result
    """
    # CONSTRUCTOR
    def __init__(self, genes):
        self.genes = genes
        self.log_likelihood_result = None
        self.relative_likelihood_result = None
        self.fitness = None

    # ACCESSOR METHODS
    def get_genes(self):
        return self.genes

    def set_genes(self, genes):
        self.genes = genes

    def get_log_likelihood_result(self):
        return self.log_likelihood_result

    def set_log_likelihood_result(self, log_likelihood_result):
        self.log_likelihood_result = log_likelihood_result

    def get_relative_likelihood_result(self):
        return self.relative_likelihood_result

    def set_relative_likelihood_result(self, relative_likelihood_result):
        self.relative_likelihood_result = relative_likelihood_result

    def get_fitness(self):
        return self.fitness

    def set_fitness(self, fitness):
        self.fitness = fitness

    # METHODS
    def bit_changer(self, i, j):
        """
        Changes a bit given a position of the matrix of genes (i,j). If the number is 1 it will be changed for 0 and
        vice versa.

        Args:
            i : INT
                Position x of the matrix (row)
            j : INT
                Position y of the matrix (column)
        """
        self.genes[i][j] = (1 - self.genes[i][j])

    def counter_ones(self):
        """
        Returns the number of 1 that it has into the matrix of genes.

        Returns:
            INT
                A number representing the number of 1 that the matrix of genes has
        """
        num_ones = 0
        for i in range(0, len(self.genes)):
            for j in range(0, len(self.genes)):
                num_ones += self.genes[i][j]
        return num_ones

    def __repr__(self):
        msg = ""
        for i in range(0, len(self.genes)):
            msg += str(self.genes[i]) + "\n"
        msg += "Log Likelihood Result: " + str(self.log_likelihood_result) + "\n"
        msg += "Relative Likelihood Result: " + str(self.relative_likelihood_result) + "\n"
        msg += "Fitness: " + str(self.fitness) + "\n"
        msg += "Num. Ones: " + str(self.counter_ones())
        return msg
