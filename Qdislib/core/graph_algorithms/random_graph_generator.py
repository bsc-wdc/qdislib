#!/usr/bin/env python3
#
#  Copyright 2002-2025 Barcelona Supercomputing Center (www.bsc.es)
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

"""Graph algorithms."""

import networkx
import random
import matplotlib.pyplot as plt


def generate_random_graph(
    num_nodes: int, probability: float = 0.2, seed: int = 10, draw: bool = False
) -> networkx.Graph:
    """Generate a random undirected graph with a specified number of nodes.

    :param num_nodes: The number of nodes in the graph.
    :param probability: The probability of creating an edge between any pair of nodes.
    :param seed: Random seed.
    :param draw: Plot the generated random graph.
    :return: The generated graph with num_nodes nodes.
    """
    random.seed(seed)

    # Create an empty graph
    random_graph = networkx.Graph()

    # Add nodes to the graph
    random_graph.add_nodes_from(range(num_nodes))

    # Add edges randomly based on the given probability
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            if random.random() < probability:
                random_graph.add_edge(i, j)

    # Ensure the graph is connected (Optional step)
    # This step ensures that the generated graph is a single connected component
    if not networkx.is_connected(random_graph) and num_nodes > 1:
        # Keep adding random edges until the graph becomes connected
        while not networkx.is_connected(random_graph):
            u, v = random.sample(range(num_nodes), 2)
            random_graph.add_edge(u, v)

    if draw:
        plt.figure(figsize=(8, 6))
        pos = networkx.spring_layout(
            random_graph
        )  # Use spring layout for a visually pleasing graph
        networkx.draw(
            random_graph,
            pos,
            with_labels=True,
            node_color="skyblue",
            node_size=500,
            font_size=10,
            font_weight="bold",
        )
        plt.show()

    return random_graph
