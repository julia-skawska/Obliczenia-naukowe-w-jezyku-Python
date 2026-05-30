#Importy
"""

import networkx as nx
import matplotlib.pyplot as plt

"""
## Exe 12.2

G = nx.erdos_renyi_graph(15, 0.3)

pos = nx.spring_layout(G)

degrees = dict(G.degree())
node_colors = [degrees[n] for n in G.nodes]
node_sizes = [100 + 100 * degrees[n] for n in G.nodes]

nx.draw(
    G,
    pos=pos,
    node_color=node_colors,
    node_size=node_sizes,
    edge_color="black",
    cmap=plt.cm.viridis
)

plt.title("random undirected graph G ")
plt.show()
