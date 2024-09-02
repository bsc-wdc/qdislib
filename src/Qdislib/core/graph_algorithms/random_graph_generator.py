#!/usr/bin/python
#
#  Copyright 2002-2024 Barcelona Supercomputing Center (www.bsc.es)
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#

# -*- coding: utf-8 -*-

import networkx as nx
import random
import matplotlib.pyplot as plt

def generate_random_graph(num_nodes, probability=0.2, seed = 10 ,draw=False):
    """
    Generates a random undirected graph with a specified number of nodes.
    
    Parameters:
    - num_nodes: int - The number of nodes in the graph.
    - probability: float - The probability of creating an edge between any pair of nodes.
    
    Returns:
    - G: networkx.Graph - The generated graph with num_nodes nodes.
    """
    random.seed(seed)

    # Create an empty graph
    G = nx.Graph()
    
    # Add nodes to the graph
    G.add_nodes_from(range(num_nodes))
    
    # Add edges randomly based on the given probability
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if random.random() < probability:
                G.add_edge(i, j)
    
    # Ensure the graph is connected (Optional step)
    # This step ensures that the generated graph is a single connected component
    if not nx.is_connected(G) and num_nodes > 1:
        # Keep adding random edges until the graph becomes connected
        while not nx.is_connected(G):
            u, v = random.sample(range(num_nodes), 2)
            G.add_edge(u, v)

    if draw:
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(G)  # Use spring layout for a visually pleasing graph
        nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=500, font_size=10, font_weight='bold')
        plt.show()
    
    return G
    