import networkx as nx
from scipy.spatial.distance import cdist
import numpy as np
import anndata
import random
import pandas as pd
import math
from scipy.spatial import ConvexHull
from matplotlib.path import Path

# from mapped_region import plot_slice_mapping


def create_graph(adata, degree=4):
    """
    Converts spatial coordinates into graph using networkx library.

    param: adata - ST Slice
    param: degree - number of edges per vertex

    return: 1) G - networkx graph
            2) node_dict - dictionary mapping nodes to spots
    """
    D = cdist(adata.obsm['spatial'], adata.obsm['spatial'])
    # Get column indexes of the degree+1 lowest values per row
    idx = np.argsort(D, 1)[:, 0:degree + 1]
    # Remove first column since it results in self loops
    idx = idx[:, 1:]

    G = nx.Graph()
    for r in range(len(idx)):
        for c in idx[r]:
            G.add_edge(r, c)

    node_dict = dict(zip(range(adata.shape[0]), adata.obs.index))
    return G, node_dict


def generate_graph_from_labels(adata, labels_dict):
    """
    Creates and returns the graph and dictionary {node: cluster_label}
    """
    g, node_to_spot = create_graph(adata)
    spot_to_cluster = labels_dict

    # remove any nodes that are not mapped to a cluster
    removed_nodes = []
    for node in node_to_spot.keys():
        if (node_to_spot[node] not in spot_to_cluster.keys()):
            removed_nodes.append(node)

    for node in removed_nodes:
        del node_to_spot[node]
        g.remove_node(node)

    labels = dict(zip(g.nodes(), [spot_to_cluster[node_to_spot[node]] for node in g.nodes()]))
    return g, labels


def spatial_coherence_score(graph, labels):
    g, l = graph, labels
    true_entropy = spatial_entropy(g, l)
    entropies = []
    for i in range(500):
        new_l = list(l.values())
        random.shuffle(new_l)
        labels = dict(zip(l.keys(), new_l))
        entropies.append(spatial_entropy(g, labels))

    return (true_entropy - np.mean(entropies)) / np.std(entropies)


def spatial_entropy(g, labels):
    """
    Calculates spatial entropy of graph
    """
    # construct contiguity matrix C which counts pairs of cluster edges
    cluster_names = np.unique(list(labels.values()))
    C = pd.DataFrame(0,index=cluster_names, columns=cluster_names)

    for e in g.edges():
        C[labels[e[0]]][labels[e[1]]] += 1

    # calculate entropy from C
    C_sum = C.values.sum()
    H = 0
    for i in range(len(cluster_names)):
        for j in range(i, len(cluster_names)):
            if (i == j):
                z = C[cluster_names[i]][cluster_names[j]]
            else:
                z = C[cluster_names[i]][cluster_names[j]] + C[cluster_names[j]][cluster_names[i]]
            if z != 0:
                H += -(z/C_sum)*math.log(z/C_sum)
    return H


def nonconsistent_edge_score(graph, labels):
    g, l = graph, labels
    true_measure = nonconsistent_edge_measure(g, l)
    measures = []
    for i in range(200):
        new_l = list(l.values())
        random.shuffle(new_l)
        labels = dict(zip(l.keys(), new_l))
        measures.append(nonconsistent_edge_measure(g, labels))

    return (true_measure - np.mean(measures)) / np.std(measures)


def nonconsistent_edge_measure(g, labels):
    # construct contiguity matrix C which counts pairs of cluster edges
    cluster_names = np.unique(list(labels.values()))
    C = pd.DataFrame(0, index=cluster_names, columns=cluster_names)

    for e in g.edges():
        C[labels[e[0]]][labels[e[1]]] += 1

    C_sum = C.values.sum()
    diagonal = 0
    for i in range(len(cluster_names)):
        diagonal += C[cluster_names[i]][cluster_names[i]]

    return float(C_sum - diagonal) / C_sum


def calculate_nonconsistent_edge_measure_for_alignment(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    target_split = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split

    g_A, l_A = generate_graph_from_labels(sliceA, sliceA.obs['aligned'])
    measure_A = nonconsistent_edge_measure(g_A, l_A)

    g_B, l_B = generate_graph_from_labels(sliceB, sliceB.obs['aligned'])
    measure_B = nonconsistent_edge_measure(g_B, l_B)
    return measure_A, measure_B


def calculate_convex_hull_noncosistent_edge_measure(sliceA, sliceB, pi):
    sliceA = sliceA.copy()

    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    source_mapped_points = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_mapped_points.append(sliceA.obsm['spatial'][i])
    source_mapped_points = np.array(source_mapped_points)
    source_hull = ConvexHull(source_mapped_points)
    source_hull_path = Path(source_mapped_points[source_hull.vertices])
    source_hull_adata = sliceA[sliceA.obs.index[source_hull_path.contains_points(sliceA.obsm['spatial'])]]

    g_A, l_A = generate_graph_from_labels(source_hull_adata, source_hull_adata.obs['aligned'])
    measure_A = nonconsistent_edge_measure(g_A, l_A)


    sliceB = sliceB.copy()

    target_split = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split

    target_mapped_points = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_mapped_points.append(sliceB.obsm['spatial'][i])
    target_mapped_points = np.array(target_mapped_points)
    target_hull = ConvexHull(target_mapped_points)
    target_hull_path = Path(target_mapped_points[target_hull.vertices])
    target_hull_adata = sliceB[sliceB.obs.index[target_hull_path.contains_points(sliceB.obsm['spatial'])]]

    g_B, l_B = generate_graph_from_labels(target_hull_adata, target_hull_adata.obs['aligned'])
    measure_B = nonconsistent_edge_measure(g_B, l_B)

    return measure_A, measure_B


def calculate_nonconsistent_edge_score_for_alignment(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    target_split = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split

    g_A, l_A = generate_graph_from_labels(sliceA, sliceA.obs['aligned'])
    score_A = np.abs(nonconsistent_edge_score(g_A, l_A))

    g_B, l_B = generate_graph_from_labels(sliceB, sliceB.obs['aligned'])
    score_B = np.abs(nonconsistent_edge_score(g_B, l_B))
    return score_A, score_B


def calculate_spatial_coherence_for_alignment(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    target_split = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split

    g_A, l_A = generate_graph_from_labels(sliceA, sliceA.obs['aligned'])
    score_A = np.abs(spatial_coherence_score(g_A, l_A))

    g_B, l_B = generate_graph_from_labels(sliceB, sliceB.obs['aligned'])
    score_B = np.abs(spatial_coherence_score(g_B, l_B))
    return score_A, score_B

    # g_A, l_A = generate_graph_from_labels(sliceA, sliceA.obs['aligned'])
    # entropy_A = spatial_entropy(g_A, l_A)
    #
    # g_B, l_B = generate_graph_from_labels(sliceB, sliceB.obs['aligned'])
    # entropy_B = spatial_entropy(g_B, l_B)
    # return entropy_A, entropy_B


def calculate_convex_hull_spatial_entropy(sliceA, sliceB, pi):
    sliceA = sliceA.copy()

    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    source_mapped_points = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_mapped_points.append(sliceA.obsm['spatial'][i])
    source_mapped_points = np.array(source_mapped_points)
    source_hull = ConvexHull(source_mapped_points)
    source_hull_path = Path(source_mapped_points[source_hull.vertices])
    source_hull_adata = sliceA[sliceA.obs.index[source_hull_path.contains_points(sliceA.obsm['spatial'])]]

    g_A, l_A = generate_graph_from_labels(source_hull_adata, source_hull_adata.obs['aligned'])
    entropy_A = spatial_entropy(g_A, l_A)


    sliceB = sliceB.copy()

    target_split = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split

    target_mapped_points = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_mapped_points.append(sliceB.obsm['spatial'][i])
    target_mapped_points = np.array(target_mapped_points)
    target_hull = ConvexHull(target_mapped_points)
    target_hull_path = Path(target_mapped_points[target_hull.vertices])
    target_hull_adata = sliceB[sliceB.obs.index[target_hull_path.contains_points(sliceB.obsm['spatial'])]]

    g_B, l_B = generate_graph_from_labels(target_hull_adata, target_hull_adata.obs['aligned'])
    entropy_B = spatial_entropy(g_B, l_B)

    return entropy_A, entropy_B


def calculate_convex_hull_spatial_coherence_score(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    source_mapped_points = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_mapped_points.append(sliceA.obsm['spatial'][i])
    source_mapped_points = np.array(source_mapped_points)
    source_hull = ConvexHull(source_mapped_points)
    source_hull_path = Path(source_mapped_points[source_hull.vertices])
    source_hull_adata = sliceA[sliceA.obs.index[source_hull_path.contains_points(sliceA.obsm['spatial'])]]

    g_A, l_A = generate_graph_from_labels(source_hull_adata, source_hull_adata.obs['aligned'])
    score_A = np.abs(spatial_coherence_score(g_A, l_A))


    sliceB = sliceB.copy()
    target_split = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split

    target_mapped_points = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_mapped_points.append(sliceB.obsm['spatial'][i])
    target_mapped_points = np.array(target_mapped_points)
    target_hull = ConvexHull(target_mapped_points)
    target_hull_path = Path(target_mapped_points[target_hull.vertices])
    target_hull_adata = sliceB[sliceB.obs.index[target_hull_path.contains_points(sliceB.obsm['spatial'])]]

    g_B, l_B = generate_graph_from_labels(target_hull_adata, target_hull_adata.obs['aligned'])
    score_B = np.abs(spatial_coherence_score(g_B, l_B))

    return score_A, score_B


def calculate_convex_hull_unaligned_percentage(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    source_mapped_points = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_mapped_points.append(sliceA.obsm['spatial'][i])
    source_mapped_points = np.array(source_mapped_points)
    source_hull = ConvexHull(source_mapped_points)
    source_hull_path = Path(source_mapped_points[source_hull.vertices])
    source_hull_adata = sliceA[sliceA.obs.index[source_hull_path.contains_points(sliceA.obsm['spatial'])]]

    source_unaligned_cnt = 0
    for i in range(source_hull_adata.obs['aligned'].values.size):
        if source_hull_adata.obs['aligned'].values[i] == "false":
            source_unaligned_cnt += 1
        elif source_hull_adata.obs['aligned'].values[i] != "true":
            print("ERROR")
            exit(1)
    source_unaligned_percentage = float(source_unaligned_cnt) / source_hull_adata.obs['aligned'].values.size

    sliceB = sliceB.copy()
    target_split = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split

    target_mapped_points = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_mapped_points.append(sliceB.obsm['spatial'][i])
    target_mapped_points = np.array(target_mapped_points)
    target_hull = ConvexHull(target_mapped_points)
    target_hull_path = Path(target_mapped_points[target_hull.vertices])
    target_hull_adata = sliceB[sliceB.obs.index[target_hull_path.contains_points(sliceB.obsm['spatial'])]]

    target_unaligned_cnt = 0
    for i in range(target_hull_adata.obs['aligned'].values.size):
        if target_hull_adata.obs['aligned'].values[i] == "false":
            target_unaligned_cnt += 1
        elif target_hull_adata.obs['aligned'].values[i] != "true":
            print("ERROR")
            exit(1)
    target_unaligned_percentage = float(target_unaligned_cnt) / target_hull_adata.obs['aligned'].values.size

    return source_unaligned_percentage, target_unaligned_percentage


def spatial_entropy_distribution(graph, labels):
    g, l = graph, labels
    true_entropy = spatial_entropy(g, l)
    print("True entropy is: " + str(true_entropy))
    entropies = []
    for i in range(500):
        new_l = list(l.values())
        random.shuffle(new_l)
        labels = dict(zip(l.keys(), new_l))
        entropies.append(spatial_entropy(g, labels))

    return entropies


def get_spatial_entropy_distribution(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    source_split = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    g_A, l_A = generate_graph_from_labels(sliceA, sliceA.obs['aligned'])
    return spatial_entropy_distribution(g_A, l_A)
