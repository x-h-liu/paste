import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import seaborn as sns


def plot_slice(slice_,figsize=None,ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0),slice_.obsm['spatial'].max(axis=0)
    len_x,len_y=max_x-min_x,max_y-min_y
    if not figsize: figsize=(10*(len_x/max(len_x,len_y)),10*(len_y/max(len_x,len_y)))
    if not ax: plt.figure(figsize=figsize)
    g = sns.scatterplot(x=slice_.obsm['spatial'][:,0],y=slice_.obsm['spatial'][:,1],linewidth=0,s=s, marker=".",ax=ax)
    if not ax: ax=g
    if ax:
        #ax.invert_yaxis()
        ax.axis('off')
    plt.show()


def plot_slice_umi(slice_, figsize=None, ax=None, s=100):
    (min_x,min_y),(max_x,max_y) = slice_.obsm['spatial'].min(axis=0), slice_.obsm['spatial'].max(axis=0)
    len_x, len_y = max_x-min_x, max_y-min_y
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
        #ax.invert_yaxis()
        ax.axis('off')


# filename = "/Users/xinhaoliu/Desktop/Research/Data/PASTE/mouse_embryo/Puck_190926_03.digital_expression.txt"
# bead_location_filename = "/Users/xinhaoliu/Desktop/Research/Data/PASTE/mouse_embryo/Puck_190926_03_bead_locations.csv"
filename = "/Users/xinhaoliu/Desktop/Research/Data/PASTE/mouse_embryo/Puck_190926_02.digital_expression.txt"
bead_location_filename = "/Users/xinhaoliu/Desktop/Research/Data/PASTE/mouse_embryo/Puck_190926_02_bead_locations.csv"
adata = sc.read_text(filename, first_column_names=True).transpose()  # TODO maybe add spot names
adata.obsm['spatial'] = np.loadtxt(bead_location_filename, delimiter=',', skiprows=1, usecols=(1, 2))
sc.pp.calculate_qc_metrics(adata, inplace=True)

# compute UMI for each spot
gene_expression_matrix = adata.X
spot_umi_counts = np.sum(gene_expression_matrix, axis=1)
adata.obs['sum_umi'] = spot_umi_counts

# # remove outlier
# to_keep = (spot_umi_counts < np.partition(spot_umi_counts, -10)[-10])
# adata = adata[adata.obs.index[to_keep]]


plot_slice_umi(adata)
plt.figure()
counts, edges, bars = plt.hist(adata.obs['sum_umi'])
plt.bar_label(bars)
plt.xlabel("UMI")
plt.ylabel("Number of spots")


plt.show()







