# Created by: Dr. David John & Kenneth Meza.
# Created at: February, 2016.
# Updated at: March, 2016.

# LIBRARIES
import pydotplus as pydot


# =======
#   GUI
# =======

# FUNCTION: create_composite_model_image
def create_model_image(di_graph, name):
    """
    Using a networkx's DAG (created based on the composite model or the amalgamated model), this function creates an
    image of the result of the Simple Genetic Algorithm using 'pydotplus' library.

    Args:
        di_graph : nx.DiGraph()
            It's a networkx's digraph representing the composite model or the amalgamated model
        name : STRING
            The name given in the variable 'a_name' on 'data_logic.data_input.py'

    Returns:
        PNG IMAGE
            An image with the result of the problem, saved on the root file
    """
    display_graph = pydot.Dot(graph_type='digraph')

    nodes = di_graph.nodes()
    for node in nodes:
        display_graph.add_node(pydot.Node(node))

    edges = di_graph.edges()
    for edge in edges:
        if di_graph[edge[0]][edge[1]]['weight'] == 1:
            display_graph.add_edge(pydot.Edge(edge[0], edge[1],
                                              label=str(di_graph[edge[0]][edge[1]]['weight']),
                                              color="green"))
        else:
            display_graph.add_edge(pydot.Edge(edge[0], edge[1],
                                              label=str(di_graph[edge[0]][edge[1]]['weight']),
                                              color="red"))
    display_graph.write_png("../" + name + ".png")


# FUNCTION: create_model_image_cotempoal
def create_model_image_cotemporal(di_graph, name):
    """
    Using a networkx's DAG (created based on the composite model or the amalgamated model), this function creates an
    image of the result of the Simple Genetic Algorithm using 'pydotplus' library. Works well for cotemporal likelihood
    function only.

    Args:
        di_graph : nx.DiGraph()
            It's a networkx's digraph representing the composite model or the amalgamated model
        name : STRING
            The name given in the variable 'a_name' on 'data_logic.data_input.py'

    Returns:
        PNG IMAGE
            An image with the result of the problem, saved on the root file
    """
    display_graph = pydot.Dot(graph_type='graph')

    nodes = di_graph.nodes()
    for node in nodes:
        display_graph.add_node(pydot.Node(node))

    edges = di_graph.edges()
    for edge in edges:
        if edge[0] < edge[1]:
            if di_graph[edge[0]][edge[1]]['weight'] == 1:
                display_graph.add_edge(pydot.Edge(edge[0], edge[1],
                                                  label=str(di_graph[edge[0]][edge[1]]['weight']),
                                                  color="green"))
            else:
                display_graph.add_edge(pydot.Edge(edge[0], edge[1],
                                                  label=str(di_graph[edge[0]][edge[1]]['weight']),
                                                  color="red"))
    display_graph.write_png("../" + name + ".png")


# FUNCTION: write_file
def write_file(file_name, text):
    """
    Writes a text on a new or created file given a name.

    Args:
        file_name : STRING
            The name of the file
        text : STRING
            The text that is going to be written on the file
    """
    file = open("../" + file_name + ".txt", "w")
    file.write(text)
    file.close()
