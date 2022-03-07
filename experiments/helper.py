import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib import colors



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
