# sbm_generator.py

import networkx as nx
import numpy as np

def generate_sbm(block_sizes, prob_matrix):
    """
    Generate a stochastic block model and return its adjacency matrix.
    
    Parameters:
    - block_sizes: list of ints (number of nodes in each block)
    - prob_matrix: 2D list or numpy array (connection probabilities between blocks)
    - seed: random seed (optional)

    Returns:
    - adj_matrix: adjacency matrix as a numpy array
    """
    G = nx.stochastic_block_model(block_sizes, prob_matrix, directed=True)
    adj_matrix = nx.to_numpy_array(G, dtype=int)
    return adj_matrix
