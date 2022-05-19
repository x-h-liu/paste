import sys
sys.path.insert(0, '/n/fs/ragr-research/users/xinhao/workspace/code/paste')
import numpy as np
from preprocessing.data_preprocessing_washu import read_HT, plot_slice
import src.paste.PASTE as paste
from src.paste.fractional_align import partial_pairwise_align
from helper import plot_slice_pairwise_alignment_mappingcolor, plot_slice_pairwise_alignment_tumor
import matplotlib.pyplot as plt
import seaborn as sns



# def ht_partial_align(sliceA, sliceB, m):
#     pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)
#     going_out = np.sum(pi, axis=1) > 0
#     coming_in = np.sum(pi, axis=0) > 0
#     going_out_part = sliceA[sliceA.obs.index[going_out]]
#     coming_in_part = sliceB[sliceB.obs.index[coming_in]]
#     plot_slice(sliceA)
#     plot_slice(going_out_part)
#     plot_slice(sliceB)
#     plot_slice(coming_in_part)


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

    
def ht_partial_align_visualization(sliceA, sliceB, pi):
    sliceA = sliceA.copy()
    sliceB = sliceB.copy()

    # pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)

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

    going_out = np.sum(pi, axis=1) > 0
    coming_in = np.sum(pi, axis=0) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]

    plot_slice_mapping(sliceA)
    plot_slice_mapping(sliceB)
    plot_slice_tumor(going_out_part)
    plot_slice_tumor(coming_in_part)
    plot_slice_pairwise_alignment_mappingcolor(sliceA, sliceB, pi)
    plot_slice_pairwise_alignment_tumor(sliceA, sliceB, pi)


def ht_partial_align(sliceA, sliceB, m):
    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)
    return pi


if __name__ == "__main__":
    slice1_ge_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/filtered_feature_bc_matrix.h5"
    slice1_loc_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice1/spatial/tissue_positions_list.csv"
    slice2_ge_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice2/filtered_feature_bc_matrix.h5"
    slice2_loc_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice2/spatial/tissue_positions_list.csv"
    slice3_ge_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice3/filtered_feature_bc_matrix.h5"
    slice3_loc_file = "/n/fs/ragr-data/datasets/DingLab/HT225/slice3/spatial/tissue_positions_list.csv"
    slice1 = read_HT(slice1_ge_file, slice1_loc_file)
    slice2 = read_HT(slice2_ge_file, slice2_loc_file)
    slice3 = read_HT(slice3_ge_file, slice3_loc_file)

    sliceA = slice2
    sliceB = slice3
    m = 0.9
    #pi, log = paste.pairwise_align(slice1, slice2, alpha=0.1, dissimilarity="kl", norm=True, return_obj=True, verbose=True)
    pi, log = partial_pairwise_align(sliceA, sliceB, alpha=0.1, m=m, armijo=False, dissimilarity="glmpca", norm=True, return_obj=True, verbose=True)
    going_out = np.sum(pi, axis=1) > 0
    coming_in = np.sum(pi, axis=0) > 0
    going_out_part = sliceA[sliceA.obs.index[going_out]]
    coming_in_part = sliceB[sliceB.obs.index[coming_in]]
    ax1 = plt.figure(1)
    plot_slice(going_out_part, ax=ax1)
    ax2 = plt.figure(2)
    plot_slice(coming_in_part, ax=ax2)

    # plot_slice_pairwise_alignment_nolabel(slice1, slice2, pi)

