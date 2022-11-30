import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sklearn
import matplotlib.patches as mpatches
from pamona import Pamona
import tangram as tg

from src.paste.fractional_align import partial_pairwise_align
from experiments.helper import plot_slice_pairwise_alignment_mappingcolor, plot_slice_pairwise_alignment_mouse
from src.paste.visualization import partial_stack_slices_pairwise, stack_slices_pairwise
import src.paste.PASTE as paste
from src.paste.helper import to_dense_array, extract_data_matrix


cell_to_color_map = {}
cell_to_color_map['neuroplacodal cell'] = sns.color_palette('tab20')[0]
cell_to_color_map['paraxial cell'] = sns.color_palette('tab20')[1]
cell_to_color_map['native cell'] = sns.color_palette('tab20')[2]
cell_to_color_map['premigratory neural crest cell'] = sns.color_palette('tab20')[3]
cell_to_color_map['midbrain dopaminergic neuron'] = sns.color_palette('tab20')[4]
cell_to_color_map['primitive red blood cell'] = sns.color_palette('tab20')[5]
cell_to_color_map['mesodermal cell'] = sns.color_palette('tab20')[6]
cell_to_color_map['neurectodermal cell'] = sns.color_palette('tab20')[7] # E8.5 specific cell type
cell_to_color_map['endodermal cell'] = sns.color_palette('tab20')[8]
cell_to_color_map['spinal cord interneuron'] = sns.color_palette('tab20')[9]
cell_to_color_map['heart valve cell'] = sns.color_palette('tab20')[10]
cell_to_color_map['surface ectodermal cell'] = sns.color_palette('tab20')[11]
cell_to_color_map['hemangioblast'] = sns.color_palette('tab20')[12]
cell_to_color_map['gut endothelial cell'] = sns.color_palette('tab20')[13]
cell_to_color_map['splanchnic mesodermal cell'] = sns.color_palette('tab20')[14]
cell_to_color_map['notochordal cell'] = sns.color_palette('tab20')[15] # E9.5 specific cell type
def plot_slice(slice_,figsize=None,obs_name='cell_type',color_map=cell_to_color_map,ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['X_spatial'].min(axis=0),slice_.obsm['X_spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    colors = list(slice_.obs[obs_name].astype('str').map(color_map))
    g = sns.scatterplot(x = slice_.obsm['X_spatial'][:,0],y = slice_.obsm['X_spatial'][:,1],linewidth=0,s=s, marker=".",c=colors,ax=ax)
    #plt.legend(handles=[mpatches.Patch(color=cell_to_color_map[slice_.obs[obs_name].cat.categories[i]], label=slice_.obs[obs_name].cat.categories[i]) for i in range(len(slice_.obs[obs_name].cat.categories))],fontsize=10,title='Cell type',title_fontsize=15)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')


def plot_slices_overlap(slices, cell_to_color_map=cell_to_color_map):
    #plt.figure(figsize=(10,15))
    plt.figure(figsize=(10,10))
    for i in range(len(slices)):
        adata = slices[i]
        colors = list(adata.obs['cell_type'].astype('str').map(cell_to_color_map))
        plt.scatter(adata.obsm['spatial'][:,0],adata.obsm['spatial'][:,1],linewidth=0,s=100, marker=".",color=colors)
    #plt.legend(handles=[mpatches.Patch(color=cell_to_color_map[adata.obs['cell_type'].cat.categories[i]], label=adata.obs['cell_type'].cat.categories[i]) for i in range(len(adata.obs['cell_type'].cat.categories))],fontsize=10,title='Cell type',title_fontsize=15,bbox_to_anchor=(1, 1))
    plt.gca().invert_yaxis()
    plt.axis('off')
    # plt.show()


def compute_alignment_ari_mouse(sliceA, sliceB, pi):
    mapped_clusters = []
    for j in range(pi.shape[1]):
        mapping = pi[:, j]
        if np.sum(mapping) > 0:
            i = np.argmax(mapping)
            mapped_clusters.append(sliceA.obs['cell_type'][i])
        else:
            mapped_clusters.append("NULL")

    assert len(sliceB.obs['cell_type']) == len(mapped_clusters)
    true_clusters_mapped_region = []
    mapped_clusters_mapped_region = []
    for i in range(len(sliceB.obs['cell_type'])):
        if mapped_clusters[i] != "NULL":
            true_clusters_mapped_region.append(sliceB.obs['cell_type'][i])
            mapped_clusters_mapped_region.append(mapped_clusters[i])

    ari = sklearn.metrics.adjusted_rand_score(true_clusters_mapped_region, mapped_clusters_mapped_region)
    return ari


def partial_alignment_mouse(sliceA_filename, sliceB_filename, m):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.obsm['spatial'] = sliceA.obsm['X_spatial']
    sliceB.obsm['spatial'] = sliceB.obsm['X_spatial']
    plot_slice(sliceA)
    plot_slice(sliceB)

    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity='glmpca', norm=True, return_obj=True, verbose=True)
    print(pi.shape)
    print("Total mass transported is: " + str(np.sum(pi)))
    print("ARI is: " + str(compute_alignment_ari_mouse(sliceA, sliceB, pi)))

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
    plot_slice_pairwise_alignment_mouse(sliceA, sliceB, pi)

    # Stacking
    new_slices = partial_stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()
    return pi



def full_alignment_mouse(sliceA_filename, sliceB_filename):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.obsm['spatial'] = sliceA.obsm['X_spatial']
    sliceB.obsm['spatial'] = sliceB.obsm['X_spatial']
    plot_slice(sliceA)
    plot_slice(sliceB)
    pi, log = paste.pairwise_align(sliceA, sliceB, alpha=0.1, dissimilarity='kl', norm=True, return_obj=True, verbose=True)
    print("ARI is: " + str(compute_alignment_ari_mouse(sliceA, sliceB, pi)))
    plot_slice_pairwise_alignment_mouse(sliceA, sliceB, pi)

    """
    stacking
    """
    new_slices = stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()



def pamona_alignment_mouse(sliceA_filename, sliceB_filename, m):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.obsm['spatial'] = sliceA.obsm['X_spatial']
    sliceB.obsm['spatial'] = sliceB.obsm['X_spatial']
    A_spotnum = sliceA.shape[0]
    B_spotnum = sliceB.shape[0]
    n_shared = int(max(A_spotnum, B_spotnum) * m)
    A_X, B_X = to_dense_array(extract_data_matrix(sliceA, None)), to_dense_array(extract_data_matrix(sliceB, None))

    Pa = Pamona.Pamona(n_shared=[n_shared], output_dim=5)
    integrated_data, T = Pa.run_Pamona([A_X, B_X])
    pi = T[0][:A_spotnum, :B_spotnum]
    # n_shared_largest_value = np.partition(pi.flatten(), -n_shared)[-n_shared]
    # pi[pi < n_shared_largest_value] = 0

    print("ARI is: " + str(compute_alignment_ari_mouse(sliceA, sliceB, pi)))
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
    plot_slice_pairwise_alignment_mouse(sliceA, sliceB, pi)

    # Stacking
    new_slices = partial_stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()


def tangram_alignment_mouse(sliceA_filename, sliceB_filename):
    sliceA = sc.read_h5ad(sliceA_filename)
    sliceB = sc.read_h5ad(sliceB_filename)
    sliceA.obsm['spatial'] = sliceA.obsm['X_spatial']
    sliceB.obsm['spatial'] = sliceB.obsm['X_spatial']

    ad_sc = sliceA.copy()
    ad_sp = sliceB.copy()
    tg.pp_adatas(ad_sc, ad_sp, genes=None)
    ad_map = tg.map_cells_to_space(ad_sc, ad_sp, density_prior='uniform', num_epochs=500)
    pi = ad_map.X/ad_map.X.sum()

    print("ARI is: " + str(compute_alignment_ari_mouse(sliceA, sliceB, pi)))
    plot_slice_pairwise_alignment_mouse(sliceA, sliceB, pi)

    """
    stacking
    """
    new_slices = stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()

