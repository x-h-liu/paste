import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


color_map = {'true': 'orange', 'false': 'blue'}
def plot_slice_mapping(slice_,figsize=None,obs_name='aligned',color_map=color_map,ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    colors = list(slice_.obs[obs_name].astype('str').map(color_map))
    g = sns.scatterplot(x=slice_.obsm['spatial'][:, 0], y=slice_.obsm['spatial'][:, 1], linewidth=0, s=s, marker=".", c=colors, ax=ax)
    if not ax: ax=g
    if ax:
        ax.invert_yaxis()
        ax.axis('off')
    plt.show()


def mapped_region_visualization(sliceA, sliceB, pi):
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

    print("slice A")
    plot_slice_mapping(sliceA)
    print("slice B")
    plot_slice_mapping(sliceB)

