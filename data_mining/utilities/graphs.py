import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

def show_graph_with_labels(adjacency_matrix, names, lim,title=''):
    adjacency_matrix += lim*np.eye(adjacency_matrix.shape[0])
    rows, cols = np.where(adjacency_matrix < lim)
    edges = zip(rows.tolist(), cols.tolist())
    gr = nx.Graph()#.DiGraph()
    gr.add_nodes_from(range(adjacency_matrix.shape[0]))
    gr.add_edges_from(edges)
    ##
    pos = nx.spring_layout(gr)
    plt.figure()
    #nx.draw(gr, node_size=500, labels=labels, with_labels=True)
    nx.draw(gr, pos, edge_color='black', width=1, linewidths=1,\
            node_size=500, node_color='pink', alpha=0.9,\
            labels={node:names[node] for node in gr.nodes()}, with_labels=True)
    plt.title(title)
    plt.show()
    
    
    
