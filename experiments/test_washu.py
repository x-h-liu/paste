import sys
sys.path.insert(0, '/n/fs/ragr-research/users/xinhao/workspace/code/paste')

import numpy as np
import sklearn
import matplotlib.pyplot as plt
import seaborn as sns
from pamona import Pamona
import tangram as tg


from src.paste.visualization import partial_stack_slices_pairwise, stack_slices_pairwise
import src.paste.PASTE as paste
from src.paste.fractional_align import partial_pairwise_align, partial_pairwise_align_expression_and_rgb, partial_pairwise_align_given_cost_matrix
from experiments.helper import plot_slice_pairwise_alignment_mappingcolor, plot_slice_pairwise_alignment_tumor
from src.paste.helper import intersect, to_dense_array, extract_data_matrix, glmpca_distance


def plot_slice(slice_,figsize=None,ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    g = sns.scatterplot(x = slice_.obsm['spatial'][:,0],y = slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".",ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')
    plt.show()


# color_map = {'true': sns.color_palette()[0], 'false': sns.color_palette()[1]}
color_map = {'true': 'orange', 'false': 'blue'}
def plot_slice_mapping(slice_,figsize=None,obs_name='aligned',color_map=color_map,ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    colors = list(slice_.obs[obs_name].astype('str').map(color_map))
    g = sns.scatterplot(x = slice_.obsm['spatial'][:,0],y = slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".",c=colors,ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')
    plt.show()

tumor_color_map = {'Tumor': 'tomato', 'Normal': 'gray'}
def plot_slice_tumor(slice_,figsize=None,obs_name='tumor',color_map=tumor_color_map,ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    colors = list(slice_.obs[obs_name].astype('str').map(color_map))
    g = sns.scatterplot(x = slice_.obsm['spatial'][:,0],y = slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".",c=colors,ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')
    plt.show()


def plot_slices_overlap(slices, cell_to_color_map=tumor_color_map):
    #plt.figure(figsize=(10,15))
    plt.figure(figsize=(10,10))
    for i in range(len(slices)):
        adata = slices[i]
        colors = list(adata.obs['tumor'].astype('str').map(cell_to_color_map))
        plt.scatter(adata.obsm['spatial'][:,0],adata.obsm['spatial'][:,1],linewidth=0,s=100, marker=".",color=colors)
    #plt.legend(handles=[mpatches.Patch(color=cell_to_color_map[adata.obs['cell_type'].cat.categories[i]], label=adata.obs['cell_type'].cat.categories[i]) for i in range(len(adata.obs['cell_type'].cat.categories))],fontsize=10,title='Cell type',title_fontsize=15,bbox_to_anchor=(1, 1))
    plt.gca().invert_yaxis()
    plt.axis('off')
    # plt.show()


def compute_alignment_ari_ht(sliceA, sliceB, pi):
    mapped_clusters = []
    for j in range(pi.shape[1]):
        mapping = pi[:, j]
        if np.sum(mapping) > 0:
            i = np.argmax(mapping)
            mapped_clusters.append(sliceA.obs['tumor'][i])
        else:
            mapped_clusters.append("NULL")

    assert len(sliceB.obs['tumor']) == len(mapped_clusters)
    true_clusters_mapped_region = []
    mapped_clusters_mapped_region = []
    for i in range(len(sliceB.obs['tumor'])):
        if mapped_clusters[i] != "NULL":
            true_clusters_mapped_region.append(sliceB.obs['tumor'][i])
            mapped_clusters_mapped_region.append(mapped_clusters[i])

    ari = sklearn.metrics.adjusted_rand_score(true_clusters_mapped_region, mapped_clusters_mapped_region)
    return ari

    
def ht_partial_align_visualization(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    # pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)

    source_split = []
    source_mass = np.sum(pi, axis=1)
    A_max_y = np.max(np.array(sliceA[source_mass > 0].obsm['spatial'][:,1])) * 0.9
    for i in range(len(source_mass)):
        if sliceA.obsm['spatial'][i][1] < A_max_y:
            source_split.append("true")
        else:
            source_split.append("false")
    # for i in range(len(source_mass)):
    #     if source_mass[i] > 0:
    #         source_split.append("true")
    #     else:
    #         source_split.append("false")
    sliceA.obs["aligned"] = source_split

    target_split = []
    target_mass = np.sum(pi, axis=0)
    B_min_y = np.min(np.array(sliceB[target_mass > 0].obsm['spatial'][:,1])) * 1.3
    for i in range(len(target_mass)):
        if sliceB.obsm['spatial'][i][1] > B_min_y:
            target_split.append("true")
        else:
            target_split.append("false")
    # for i in range(len(target_mass)):
    #     if target_mass[i] > 0:
    #         target_split.append("true")
    #     else:
    #         target_split.append("false")
    sliceB.obs["aligned"] = target_split

    # going_out = np.sum(pi, axis=1) > 0
    # coming_in = np.sum(pi, axis=0) > 0
    # going_out_part = sliceA[sliceA.obs.index[going_out]]
    # coming_in_part = sliceB[sliceB.obs.index[coming_in]]

    plot_slice_mapping(sliceA)
    plot_slice_mapping(sliceB)
    # plot_slice_tumor(going_out_part)
    # plot_slice_tumor(coming_in_part)
    plot_slice_pairwise_alignment_mappingcolor(sliceA, sliceB, pi)
    plot_slice_pairwise_alignment_tumor(sliceA, sliceB, pi)


def fig1_visualization_leftdown(sliceA, sliceB, pi):
    source_split = []
    source_mass = np.sum(pi, axis=1)
    A_max_y = np.max(np.array(sliceA[source_mass > 0].obsm['spatial'][:,1])) * 0.9
    for i in range(len(source_mass)):
        if sliceA.obsm['spatial'][i][1] < A_max_y:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    target_split = []
    target_mass = np.sum(pi, axis=0)
    B_min_y = np.min(np.array(sliceB[target_mass > 0].obsm['spatial'][:,1])) * 1.3
    for i in range(len(target_mass)):
        if sliceB.obsm['spatial'][i][1] > B_min_y:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split


def ht_partial_align_visualization_leftup(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    source_split = []
    source_mass = np.sum(pi, axis=1)
    A_min_y = np.min(np.array(sliceA[source_mass > 0].obsm['spatial'][:,1])) * 1.3
    for i in range(len(source_mass)):
        if sliceA.obsm['spatial'][i][1] > A_min_y:
            source_split.append("true")
        else:
            source_split.append("false")
    # for i in range(len(source_mass)):
    #     if source_mass[i] > 0:
    #         source_split.append("true")
    #     else:
    #         source_split.append("false")
    sliceA.obs["aligned"] = source_split

    target_split = []
    target_mass = np.sum(pi, axis=0)
    B_max_y = np.max(np.array(sliceB[target_mass > 0].obsm['spatial'][:,1])) * 0.9
    for i in range(len(target_mass)):
        if sliceB.obsm['spatial'][i][1] < B_max_y:
            target_split.append("true")
        else:
            target_split.append("false")
    # for i in range(len(target_mass)):
    #     if target_mass[i] > 0:
    #         target_split.append("true")
    #     else:
    #         target_split.append("false")
    sliceB.obs["aligned"] = target_split

    # going_out = np.sum(pi, axis=1) > 0
    # coming_in = np.sum(pi, axis=0) > 0
    # going_out_part = sliceA[sliceA.obs.index[going_out]]
    # coming_in_part = sliceB[sliceB.obs.index[coming_in]]

    plot_slice_mapping(sliceA)
    plot_slice_mapping(sliceB)
    # plot_slice_tumor(going_out_part)
    # plot_slice_tumor(coming_in_part)
    plot_slice_pairwise_alignment_mappingcolor(sliceA, sliceB, pi)
    plot_slice_pairwise_alignment_tumor(sliceA, sliceB, pi)



def fig1_visualization_leftup(sliceA, sliceB, pi):
    source_split = []
    source_mass = np.sum(pi, axis=1)
    A_min_y = np.min(np.array(sliceA[source_mass > 0].obsm['spatial'][:,1])) * 1.3
    for i in range(len(source_mass)):
        if sliceA.obsm['spatial'][i][1] > A_min_y:
            source_split.append("true")
        else:
            source_split.append("false")
    sliceA.obs["aligned"] = source_split

    target_split = []
    target_mass = np.sum(pi, axis=0)
    B_max_y = np.max(np.array(sliceB[target_mass > 0].obsm['spatial'][:,1])) * 0.9
    for i in range(len(target_mass)):
        if sliceB.obsm['spatial'][i][1] < B_max_y:
            target_split.append("true")
        else:
            target_split.append("false")
    sliceB.obs["aligned"] = target_split



def ht_partial_align(sliceA, sliceB, m):
    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)
    return pi


def partial_alignment_exprgb_ht(sliceA, sliceB, m):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    pi, log = partial_pairwise_align_expression_and_rgb(sliceA, sliceB, alpha=0.1, m=m, armijo=False, norm=True, return_obj=True, verbose=True)
    print(pi.shape)
    print("Total mass transported is: " + str(np.sum(pi)))
    print("ARI is: " + str(compute_alignment_ari_ht(sliceA, sliceB, pi)))

    going_out = np.sum(pi, axis=1) > 0
    coming_in = np.sum(pi, axis=0) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    plot_slice_tumor(going_out_part)
    plot_slice_tumor(coming_in_part)

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
    plot_slice_pairwise_alignment_tumor(sliceA, sliceB, pi)

    # Stacking
    new_slices = partial_stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()
    return pi


def get_transport_cost_matrix(sliceA, sliceB):
    # subset for common genes
    common_genes = intersect(sliceA.var.index, sliceB.var.index)
    sliceA = sliceA[:, common_genes]
    sliceB = sliceB[:, common_genes]
    # Get transport cost matrix
    A_X, B_X = to_dense_array(extract_data_matrix(sliceA, None)), to_dense_array(extract_data_matrix(sliceB, None))
    M = glmpca_distance(A_X, B_X, latent_dim=50, filter=True, verbose=True)
    return M


def partial_alignment_exp_ht(sliceA, sliceB, m, M):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()
    plot_slice(sliceA)
    plot_slice(sliceB)

    # pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)
    pi, log = partial_pairwise_align_given_cost_matrix(sliceA, sliceB, M=M, alpha=0.1, m=m, armijo=False, norm=True, return_obj=True, verbose=True)
    print(pi.shape)
    print("Total mass transported is: " + str(np.sum(pi)))

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

    plt.show()
    return pi






def full_alignment_ht(sliceA, sliceB):
    pi, log = paste.pairwise_align(sliceA, sliceB, alpha=0.1, dissimilarity='kl', norm=True, return_obj=True, verbose=True)
    print("ARI is: " + str(compute_alignment_ari_ht(sliceA, sliceB, pi)))
    plot_slice_pairwise_alignment_tumor(sliceA, sliceB, pi)

    """
    stacking
    """
    new_slices = stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()



def pamona_alignment_ht(sliceA, sliceB, m):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    A_spotnum = sliceA.shape[0]
    B_spotnum = sliceB.shape[0]
    n_shared = int(max(A_spotnum, B_spotnum) * m)
    A_X, B_X = to_dense_array(extract_data_matrix(sliceA, None)), to_dense_array(extract_data_matrix(sliceB, None))

    Pa = Pamona.Pamona(n_shared=[n_shared], output_dim=5)
    integrated_data, T = Pa.run_Pamona([A_X, B_X])
    pi = T[0][:A_spotnum, :B_spotnum]
    # n_shared_largest_value = np.partition(pi.flatten(), -n_shared)[-n_shared]
    # pi[pi < n_shared_largest_value] = 0

    print("ARI is: " + str(compute_alignment_ari_ht(sliceA, sliceB, pi)))
    going_out = np.sum(pi, axis=1) > 0
    coming_in = np.sum(pi, axis=0) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    plot_slice_tumor(going_out_part)
    plot_slice_tumor(coming_in_part)

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
    plot_slice_pairwise_alignment_tumor(sliceA, sliceB, pi)

    # Stacking
    new_slices = partial_stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()



def tangram_alignment_ht(sliceA, sliceB):
    ad_sc = sliceA.copy()
    ad_sp = sliceB.copy()
    tg.pp_adatas(ad_sc, ad_sp, genes=None)
    ad_map = tg.map_cells_to_space(ad_sc, ad_sp, density_prior='uniform', num_epochs=500)
    pi = ad_map.X/ad_map.X.sum()

    print("ARI is: " + str(compute_alignment_ari_ht(sliceA, sliceB, pi)))
    plot_slice_pairwise_alignment_tumor(sliceA, sliceB, pi)

    """
    stacking
    """
    new_slices = stack_slices_pairwise([sliceA, sliceB], [pi])
    plot_slices_overlap(new_slices)

    plt.show()


# if __name__ == "__main__":
#     slice1_ge_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/filtered_feature_bc_matrix.h5"
#     slice1_loc_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/spatial/tissue_positions_list.csv"
#     slice2_ge_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice2/filtered_feature_bc_matrix.h5"
#     slice2_loc_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice2/spatial/tissue_positions_list.csv"
#     slice3_ge_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice3/filtered_feature_bc_matrix.h5"
#     slice3_loc_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice3/spatial/tissue_positions_list.csv"
#     slice1 = read_HT(slice1_ge_file, slice1_loc_file)
#     slice2 = read_HT(slice2_ge_file, slice2_loc_file)
#     slice3 = read_HT(slice3_ge_file, slice3_loc_file)

#     sliceA = slice2
#     sliceB = slice3
#     m = 0.9
#     #pi, log = paste.pairwise_align(slice1, slice2, alpha=0.1, dissimilarity="kl", norm=True, return_obj=True, verbose=True)
#     pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)
#     going_out = np.sum(pi, axis=1) > 0
#     coming_in = np.sum(pi, axis=0) > 0
#     going_out_part = sliceA[sliceA.obs.index[going_out]]
#     coming_in_part = sliceB[sliceB.obs.index[coming_in]]
#     ax1 = plt.figure(1)
#     plot_slice(going_out_part, ax=ax1)
#     ax2 = plt.figure(2)
#     plot_slice(coming_in_part, ax=ax2)

#     # plot_slice_pairwise_alignment_nolabel(slice1, slice2, pi)

