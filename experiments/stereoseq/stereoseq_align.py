import sys
sys.path.insert(0, '/n/fs/ragr-research/users/xinhao/workspace/code/paste')

import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sklearn
from pamona import Pamona
import tangram as tg
import matplotlib.patches as mpatches
import pandas as pd
import plotly.express as px

from src.paste.fractional_align import partial_pairwise_align
from experiments.helper import plot_slice_pairwise_alignment_mappingcolor, plot_slice_pairwise_alignment_drosophila
from src.paste.visualization import partial_stack_slices_pairwise, stack_slices_pairwise
import src.paste.PASTE as paste
from src.paste.helper import to_dense_array, extract_data_matrix
from src.paste.fractional_model_selection import decide_overlap


cell_to_color_map = {}
cell_to_color_map['CNS'] = sns.color_palette('tab20')[0]
cell_to_color_map['epidermis'] = sns.color_palette('tab20')[1]
cell_to_color_map['carcass'] = sns.color_palette('tab20')[2]
cell_to_color_map['trachea'] = sns.color_palette('tab20')[3]
cell_to_color_map['muscle'] = sns.color_palette('tab20')[4]
cell_to_color_map['midgut'] = sns.color_palette('tab20')[5]
cell_to_color_map['fat body'] = sns.color_palette('tab20')[6]
cell_to_color_map['amnioserosa'] = sns.color_palette('tab20')[7] 
cell_to_color_map['foregut'] = sns.color_palette('tab20')[8]
cell_to_color_map['salivary gland'] = sns.color_palette('tab20')[9]
def plot_slice(slice_,figsize=(8,10),obs_name='annotation',color_map=cell_to_color_map,ax=None, s=200):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    colors = list(slice_.obs[obs_name].astype('str').map(color_map))
    g = sns.scatterplot(x = slice_.obsm['spatial'][:,0],y = slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".",c=colors,ax=ax)
    #plt.legend(handles=[mpatches.Patch(color=cell_to_color_map[slice_.obs[obs_name].cat.categories[i]], label=slice_.obs[obs_name].cat.categories[i]) for i in range(len(slice_.obs[obs_name].cat.categories))],fontsize=10,title='Cell type',title_fontsize=15)
    if not ax: ax=g
    if ax:
        # ax.invert_yaxis()
        # ax.axis('off')
        ax.set(xlim=(-25,25),ylim=(-30,30))

def plot_slices_overlap(slices, cell_to_color_map=cell_to_color_map):
    #plt.figure(figsize=(10,15))
    plt.figure(figsize=(10,10))
    for i in range(len(slices)):
        adata = slices[i]
        colors = list(adata.obs['annotation'].astype('str').map(cell_to_color_map))
        plt.scatter(adata.obsm['spatial'][:,0],adata.obsm['spatial'][:,1],linewidth=0,s=200, marker=".",color=colors)
    #plt.legend(handles=[mpatches.Patch(color=cell_to_color_map[adata.obs['cell_type'].cat.categories[i]], label=adata.obs['cell_type'].cat.categories[i]) for i in range(len(adata.obs['cell_type'].cat.categories))],fontsize=10,title='Cell type',title_fontsize=15,bbox_to_anchor=(1, 1))
    # plt.gca().invert_yaxis()
    plt.axis('off')
    # plt.show()

def compute_alignment_ari_drosophila(sliceA, sliceB, pi):
    mapped_clusters = []
    for j in range(pi.shape[1]):
        mapping = pi[:, j]
        if np.sum(mapping) > 0:
            i = np.argmax(mapping)
            mapped_clusters.append(sliceA.obs['annotation'][i])
        else:
            mapped_clusters.append("NULL")

    assert len(sliceB.obs['annotation']) == len(mapped_clusters)
    true_clusters_mapped_region = []
    mapped_clusters_mapped_region = []
    for i in range(len(sliceB.obs['annotation'])):
        if mapped_clusters[i] != "NULL":
            true_clusters_mapped_region.append(sliceB.obs['annotation'][i])
            mapped_clusters_mapped_region.append(mapped_clusters[i])

    ari = sklearn.metrics.adjusted_rand_score(true_clusters_mapped_region, mapped_clusters_mapped_region)
    return ari


def decide_overlap_drosophila(sliceA_filename, sliceB_filename):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.X = sliceA.layers['raw_counts']
    sliceB.X = sliceB.layers['raw_counts']

    print(decide_overlap(sliceA, sliceB))


def partial_alignment_drosophila(sliceA_filename, sliceB_filename, m):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.X = sliceA.layers['raw_counts']
    sliceB.X = sliceB.layers['raw_counts']
    plot_slice(sliceA)
    plot_slice(sliceB)

    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity='glmpca', norm=True, return_obj=True, verbose=True)
    print(pi.shape)
    print("Total mass transported is: " + str(np.sum(pi)))
    print("ARI is: " + str(compute_alignment_ari_drosophila(sliceA, sliceB, pi)))

    going_out = np.sum(pi, axis=1) > 0
    coming_in = np.sum(pi, axis=0) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    plot_slice(going_out_part)
    plot_slice(coming_in_part)

    # Alignment visualization
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
    plot_slice_pairwise_alignment_mappingcolor(sliceA, sliceB, pi)
    plot_slice_pairwise_alignment_drosophila(sliceA, sliceB, pi)

    # Stacking
    new_slices = partial_stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()
    return pi


def full_alignment_drosophila(sliceA_filename, sliceB_filename):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.X = sliceA.layers['raw_counts']
    sliceB.X = sliceB.layers['raw_counts']
    plot_slice(sliceA)
    plot_slice(sliceB)
    pi, log = paste.pairwise_align(sliceA, sliceB, alpha=0.1, dissimilarity='kl', norm=True, return_obj=True, verbose=True)
    print("ARI is: " + str(compute_alignment_ari_drosophila(sliceA, sliceB, pi)))
    plot_slice_pairwise_alignment_drosophila(sliceA, sliceB, pi)

    """
    stacking
    """
    new_slices = stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()


def pamona_alignment_drosophila(sliceA_filename, sliceB_filename, m):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.X = sliceA.layers['raw_counts']
    sliceB.X = sliceB.layers['raw_counts']
    A_spotnum = sliceA.shape[0]
    B_spotnum = sliceB.shape[0]
    n_shared = int(max(A_spotnum, B_spotnum) * m)
    A_X, B_X = to_dense_array(extract_data_matrix(sliceA, None)), to_dense_array(extract_data_matrix(sliceB, None))

    Pa = Pamona.Pamona(n_shared=[n_shared], output_dim=5)
    integrated_data, T = Pa.run_Pamona([A_X, B_X])
    pi = T[0][:A_spotnum, :B_spotnum]
    # n_shared_largest_value = np.partition(pi.flatten(), -n_shared)[-n_shared]
    # pi[pi < n_shared_largest_value] = 0

    print("ARI is: " + str(compute_alignment_ari_drosophila(sliceA, sliceB, pi)))
    going_out = np.sum(pi, axis=1) > 0
    coming_in = np.sum(pi, axis=0) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    plot_slice(going_out_part)
    plot_slice(coming_in_part)

    # Alignment visualization
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
    plot_slice_pairwise_alignment_mappingcolor(sliceA, sliceB, pi)
    plot_slice_pairwise_alignment_drosophila(sliceA, sliceB, pi)

    # Stacking
    new_slices = partial_stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()


def tangram_alignment_drosophila(sliceA_filename, sliceB_filename):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.X = sliceA.layers['raw_counts']
    sliceB.X = sliceB.layers['raw_counts']

    ad_sc = sliceA.copy()
    ad_sp = sliceB.copy()
    tg.pp_adatas(ad_sc, ad_sp, genes=None)
    ad_map = tg.map_cells_to_space(ad_sc, ad_sp, density_prior='uniform', num_epochs=500)
    pi = ad_map.X/ad_map.X.sum()

    print("ARI is: " + str(compute_alignment_ari_drosophila(sliceA, sliceB, pi)))
    plot_slice_pairwise_alignment_drosophila(sliceA, sliceB, pi)

    """
    stacking
    """
    new_slices = stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()




def partial_alignment_pie_chart(sliceA_filename, sliceB_filename, pi):
    cell_to_color_map = {}
    cell_to_color_map['CNS'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[0])
    cell_to_color_map['epidermis'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[1])
    cell_to_color_map['carcass'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[2])
    cell_to_color_map['trachea'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[3])
    cell_to_color_map['muscle'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[4])
    cell_to_color_map['midgut'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[5])
    cell_to_color_map['fat body'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[6])
    cell_to_color_map['amnioserosa'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[7])
    cell_to_color_map['foregut'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[8])
    cell_to_color_map['salivary gland'] = matplotlib.colors.to_hex(sns.color_palette('tab20')[9])

    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.X = sliceA.layers['raw_counts']
    sliceB.X = sliceB.layers['raw_counts']
    df_A = pd.DataFrame(columns=['Annotation','Count'])
    for i in range(sliceA.shape[0]):
        df_A.loc[len(df_A)] = [sliceA.obs['annotation'][i], 1]
    fig_A = px.pie(df_A, values='Count', names='Annotation', color='Annotation', color_discrete_map=cell_to_color_map)
    fig_A.show()
    df_B = pd.DataFrame(columns=['Annotation','Count'])
    for i in range(sliceB.shape[0]):
        df_B.loc[len(df_B)] = [sliceB.obs['annotation'][i], 1]
    fig_B = px.pie(df_B, values='Count', names='Annotation', color='Annotation', color_discrete_map=cell_to_color_map)
    fig_B.show()

    # pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity='glmpca', norm=True, return_obj=True, verbose=True)
    going_out = np.sum(pi, axis=1) > 0
    coming_in = np.sum(pi, axis=0) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]

    df_A_aligned = pd.DataFrame(columns=['Annotation','Count'])
    for i in range(going_out_part.shape[0]):
        df_A_aligned.loc[len(df_A_aligned)] = [going_out_part.obs['annotation'][i], 1]
    fig_A_aligned = px.pie(df_A_aligned, values='Count', names='Annotation', color='Annotation', color_discrete_map=cell_to_color_map)
    fig_A_aligned.show()
    df_B_aligned = pd.DataFrame(columns=['Annotation','Count'])
    for i in range(coming_in_part.shape[0]):
        df_B_aligned.loc[len(df_B_aligned)] = [coming_in_part.obs['annotation'][i], 1]
    fig_B_aligned = px.pie(df_B_aligned, values='Count', names='Annotation', color='Annotation', color_discrete_map=cell_to_color_map)
    fig_B_aligned.show()




if __name__ == "__main__":
    filename1 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s01.h5ad"
    filename2 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s02.h5ad"
    filename3 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s03.h5ad"
    filename4 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s04.h5ad"
    filename5 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s05.h5ad"
    filename6 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s06.h5ad"
    filename7 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s07.h5ad"
    filename8 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s08.h5ad"
    filename9 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s09.h5ad"
    filename10 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s10.h5ad"
    filename11 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s11.h5ad"
    filename12 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s12.h5ad"
    filename13 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s13.h5ad"
    filename14 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s14.h5ad"
    filename15 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s15.h5ad"
    filename16 = "/n/fs/ragr-data/users/xinhao/stereoseq/E14-16h/s16.h5ad"

    partial_alignment_pie_chart(filename7, filename8, m=0.7)
