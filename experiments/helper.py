import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
from matplotlib import colors
import sklearn



layer_to_color_map = {'Layer{0}'.format(i+1):sns.color_palette()[i] for i in range(6)}
layer_to_color_map['WM'] = sns.color_palette()[6]
def plot_slice(slice_,figsize=None,obs_name='layer_guess_reordered',color_map=layer_to_color_map,ax=None, s=100):
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


def plot_slice_umi(slice_, figsize=None, ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    min_umi = np.min(slice_.obs["sum_umi"])
    max_umi = np.max(slice_.obs["sum_umi"])
    norm = colors.Normalize(vmin=min_umi, vmax=max_umi, clip=True)
    # cmap = sns.color_palette("rocket_r", as_cmap=True)
    cmap = sns.color_palette("YlOrBr", as_cmap=True)
    c = [cmap(norm(i)) for i in slice_.obs["sum_umi"]]
    g = sns.scatterplot(x = slice_.obsm['spatial'][:,0],y = slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".", c=c, ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')


def plot_slices_overlap(slices, layer_to_color_map=layer_to_color_map):
    #plt.figure(figsize=(10,15))
    plt.figure(figsize=(10,10))
    for i in range(len(slices)):
        adata = slices[i]
        colors = list(adata.obs['layer_guess_reordered'].astype('str').map(layer_to_color_map))
        plt.scatter(adata.obsm['spatial'][:,0],adata.obsm['spatial'][:,1],linewidth=0,s=100, marker=".",color=colors)
    plt.legend(handles=[mpatches.Patch(color=layer_to_color_map[adata.obs['layer_guess_reordered'].cat.categories[i]], label=adata.obs['layer_guess_reordered'].cat.categories[i]) for i in range(len(adata.obs['layer_guess_reordered'].cat.categories))],fontsize=10,title='Cortex layer',title_fontsize=15,bbox_to_anchor=(1, 1))
    plt.gca().invert_yaxis()
    plt.axis('off')
    plt.show()


def compute_alignment_ari(sliceA, sliceB, pi):
    mapped_clusters = []
    for j in range(pi.shape[1]):
        mapping = pi[:, j]
        if np.sum(mapping) > 0:
            i = np.argmax(mapping)
            mapped_clusters.append(sliceA.obs['layer_guess_reordered'][i])
        else:
            mapped_clusters.append("NULL")

    assert len(sliceB.obs['layer_guess_reordered']) == len(mapped_clusters)
    true_clusters_mapped_region = []
    mapped_clusters_mapped_region = []
    for i in range(len(sliceB.obs['layer_guess_reordered'])):
        if mapped_clusters[i] != "NULL":
            true_clusters_mapped_region.append(sliceB.obs['layer_guess_reordered'][i])
            mapped_clusters_mapped_region.append(mapped_clusters[i])

    ari = sklearn.metrics.adjusted_rand_score(true_clusters_mapped_region, mapped_clusters_mapped_region)
    return ari



"""
Code for visualizing alignment
"""
def largest_indices(ary, n):
    """Returns the n largest indices from a numpy array."""
    flat = ary.flatten()
    indices = np.argpartition(flat, -n)[-n:]
    indices = indices[np.argsort(-flat[indices])]
    return np.unravel_index(indices, ary.shape)


def plot2D_samples_mat(xs, xt, G, thr=1e-8, alpha=0.2, top=1000, weight_alpha=False, **kwargs):
    if ('color' not in kwargs) and ('c' not in kwargs):
        kwargs['color'] = 'k'
    mx = G.max()
    #     idx = np.where(G/mx>=thr)
    idx = largest_indices(G, top)
    print(len(idx[0]))
    for l in range(len(idx[0])):
        plt.plot([xs[idx[0][l], 0], xt[idx[1][l], 0]], [xs[idx[0][l], 1], xt[idx[1][l], 1]],
                 alpha=alpha * (1 - weight_alpha) + (weight_alpha * G[idx[0][l], idx[1][l]] / mx), c='k')


def plot_slice_pairwise_alignment(slice1, slice2, pi, thr=1 - 1e-8, alpha=0.05, top=1000, name='', save=False,
                                  weight_alpha=False):
    coordinates1, coordinates2 = slice1.obsm['spatial'], slice2.obsm['spatial']
    offset = (coordinates1[:, 0].max() - coordinates2[:, 0].min()) * 1.1
    temp = np.zeros(coordinates2.shape)
    temp[:, 0] = offset
    plt.figure(figsize=(20, 10))
    plot2D_samples_mat(coordinates1, coordinates2 + temp, pi, thr=thr, c='k', alpha=alpha, top=top,
                       weight_alpha=weight_alpha)
    plt.scatter(coordinates1[:, 0], coordinates1[:, 1], linewidth=0, s=100, marker=".", color=list(
        slice1.obs['layer_guess_reordered'].map(
            dict(zip(slice1.obs['layer_guess_reordered'].cat.categories, slice1.uns['layer_guess_reordered_colors'])))))
    plt.scatter(coordinates2[:, 0] + offset, coordinates2[:, 1], linewidth=0, s=100, marker=".", color=list(
        slice2.obs['layer_guess_reordered'].map(
            dict(zip(slice2.obs['layer_guess_reordered'].cat.categories, slice2.uns['layer_guess_reordered_colors'])))))
    plt.gca().invert_yaxis()
    plt.axis('off')
    if save: plt.savefig("/Users/ronzeira/Documents/GitHub/paste/paste_output/DLPFC/figs/{}.pdf".format(name),
                         bbox_inches='tight')
    plt.show()


def plot_slice_pairwise_alignment_nolabel(slice1, slice2, pi, thr=1 - 1e-8, alpha=0.05, top=1000, name='', save=False,
                                  weight_alpha=False):
    coordinates1, coordinates2 = slice1.obsm['spatial'], slice2.obsm['spatial']
    offset = (coordinates1[:, 0].max() - coordinates2[:, 0].min()) * 1.1
    temp = np.zeros(coordinates2.shape)
    temp[:, 0] = offset
    plt.figure(figsize=(20, 10))
    plot2D_samples_mat(coordinates1, coordinates2 + temp, pi, thr=thr, c='k', alpha=alpha, top=top,
                       weight_alpha=weight_alpha)
    plt.scatter(coordinates1[:, 0], coordinates1[:, 1], linewidth=0, s=100, marker=".")
    plt.scatter(coordinates2[:, 0] + offset, coordinates2[:, 1], linewidth=0, s=100, marker=".")
    plt.gca().invert_yaxis()
    plt.axis('off')
    if save: plt.savefig("/Users/ronzeira/Documents/GitHub/paste/paste_output/DLPFC/figs/{}.pdf".format(name),
                         bbox_inches='tight')
    plt.show()


# align_region_color_map = {'true': sns.color_palette()[0], 'false': sns.color_palette()[1]}
align_region_color_map = {'true': 'orange', 'false': 'blue'}
def plot_slice_pairwise_alignment_mappingcolor(slice1, slice2, pi, thr=1 - 1e-8, alpha=0.05, top=1000, name='', save=False,
                                  weight_alpha=False):
    coordinates1, coordinates2 = slice1.obsm['spatial'], slice2.obsm['spatial']
    offset = (coordinates1[:, 0].max() - coordinates2[:, 0].min()) * 1.1
    temp = np.zeros(coordinates2.shape)
    temp[:, 0] = offset
    plt.figure(figsize=(20, 10))
    plot2D_samples_mat(coordinates1, coordinates2 + temp, pi, thr=thr, c='k', alpha=alpha, top=top,
                       weight_alpha=weight_alpha)
    plt.scatter(coordinates1[:, 0], coordinates1[:, 1], linewidth=0, s=100, marker=".", color=list(slice1.obs['aligned'].astype('str').map(align_region_color_map)))
    plt.scatter(coordinates2[:, 0] + offset, coordinates2[:, 1], linewidth=0, s=100, marker=".", color=list(slice2.obs['aligned'].astype('str').map(align_region_color_map)))
    plt.gca().invert_yaxis()
    plt.axis('off')
    if save: plt.savefig("/Users/ronzeira/Documents/GitHub/paste/paste_output/DLPFC/figs/{}.pdf".format(name),
                         bbox_inches='tight')
    plt.show()


tumor_color_map = {'Tumor': 'tomato', 'Normal': 'gray'}
def plot_slice_pairwise_alignment_tumor(slice1, slice2, pi, thr=1 - 1e-8, alpha=0.05, top=1000, name='', save=False,
                                  weight_alpha=False):
    coordinates1, coordinates2 = slice1.obsm['spatial'], slice2.obsm['spatial']
    offset = (coordinates1[:, 0].max() - coordinates2[:, 0].min()) * 1.1
    temp = np.zeros(coordinates2.shape)
    temp[:, 0] = offset
    plt.figure(figsize=(20, 10))
    plot2D_samples_mat(coordinates1, coordinates2 + temp, pi, thr=thr, c='k', alpha=alpha, top=top,
                       weight_alpha=weight_alpha)
    plt.scatter(coordinates1[:, 0], coordinates1[:, 1], linewidth=0, s=100, marker=".", color=list(slice1.obs['tumor'].astype('str').map(tumor_color_map)))
    plt.scatter(coordinates2[:, 0] + offset, coordinates2[:, 1], linewidth=0, s=100, marker=".", color=list(slice2.obs['tumor'].astype('str').map(tumor_color_map)))
    plt.gca().invert_yaxis()
    plt.axis('off')
    if save: plt.savefig("/Users/ronzeira/Documents/GitHub/paste/paste_output/DLPFC/figs/{}.pdf".format(name),
                         bbox_inches='tight')
    plt.show()
"""
Code for visualizing alignment ends
"""
