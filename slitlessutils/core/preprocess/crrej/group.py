import math
import numpy as np
from scipy.spatial import distance
from sklearn.cluster import AgglomerativeClustering

from slitlessutils.core.wfss import WFSSCollection
from ....logger import LOGGER


def angular_distance(angle1, angle2, degrees=True):
    '''
    Return angular distance between two angles
    '''
    if degrees:
        circumference = 360
    else:
        circumference = 2 * math.pi
    arc_length = abs(angle2 - angle1)
    distance = min(arc_length, circumference - arc_length)
    return distance


def group_by_position_angle(files, degrees=True, max_pa_diff=0.2, **kwargs):
    '''
    Group input files by position angle using Agglomerative Clustering, for input
    to cosmic ray rejection routine.

    Parameters
    ----------
    files : list of str
        List of input filenames to be grouped for cosmic ray rejection.
    degrees : bool, optional
        Set to False if PA is in radians instead of degrees.
    max_pa_diff : float, optional
        The maximium difference between PAs for two files to be considered
        part of the same group.
    '''
    data_collection = WFSSCollection.from_list(files)
    position_angles = data_collection.get_pas()
    # Input needs to be 2D
    position_angles = np.reshape(position_angles, (len(position_angles), 1))

    # Precompute distance matrix for input to clustering
    distance_matrix = distance.pdist(position_angles, angular_distance)
    # We need an NxN matrix instead of the condensed distance matrix
    distance_matrix = distance.squareform(distance_matrix)

    # Fit clustering model. Resulting clusters are integer values in model.labels_.
    model = AgglomerativeClustering(n_clusters=None, linkage="single", metric="precomputed",
                                    distance_threshold=max_pa_diff)
    model.fit(distance_matrix)

    # Return list of grouped filenames
    grouped_files = []
    for i in range(0, np.max(model.labels_)):
        members = files[np.where(model.labels_ == i)]
        if len(members) == 1:
            LOGGER.warning(f"The file: {members[0]} was not grouped with any others.")
        grouped_files.append()

    return grouped_files
