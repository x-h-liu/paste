import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from matplotlib.path import Path

from mapped_region import plot_slice_mapping

def convex_hull_area(sliceA, sliceB, pi):
    source_mapped_points = []
    source_mass = np.sum(pi, axis=1)
    for i in range(len(source_mass)):
        if source_mass[i] > 0:
            source_mapped_points.append(sliceA.obsm['spatial'][i])
    source_mapped_points = np.array(source_mapped_points)
    source_mapped_region_hull = ConvexHull(source_mapped_points)
    source_mapped_region_hull_area = source_mapped_region_hull.volume

    target_mapped_points = []
    target_mass = np.sum(pi, axis=0)
    for i in range(len(target_mass)):
        if target_mass[i] > 0:
            target_mapped_points.append(sliceB.obsm['spatial'][i])
    target_mapped_points = np.array(target_mapped_points)
    target_mapped_region_hull = ConvexHull(target_mapped_points)
    target_mapped_region_hull_area = target_mapped_region_hull.volume

    return source_mapped_region_hull_area, target_mapped_region_hull_area


def visualize_convex_hull(sliceA, sliceB, pi):
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

    plot_slice_mapping(source_hull_adata)
    plot_slice_mapping(target_hull_adata)



