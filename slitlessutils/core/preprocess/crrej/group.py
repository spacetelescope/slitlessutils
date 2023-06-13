import numpy as np

from slitlessutils.core.wfss import WFSSCollection
from ....logger import LOGGER


def group_by_visit(files, return_unique_visits=False, **kwargs):
    """
    Group a list of files based on their unique visits.

    Parameters
    ----------
    files : list of str
        List of filenames

    return_unique_visits : bool, optional
        Indicating whether to return unique visits along with grouped files.
        If `True`, a tuple containing grouped files and unique visits is returned.
        If `False`, only grouped files are returned. Default is `False`.

    **kwargs
        Additional keyword arguments

    Returns
    ----------
    grouped_files : list or tuple
        Grouped files based on visits. If `return_unique_visits` is `True`, a tuple
        containing grouped files and unique visits is returned. If `return_unique_visits`
        is `False`, only grouped files are returned.
    """

    data_collection = WFSSCollection.from_list(files)
    visits = data_collection.get_visits()

    visits = np.asarray(visits)
    files = np.asarray(files)

    unique_visits, indices, counts = np.unique(
        visits, return_inverse=True, return_counts=True
    )
    rev = np.split(np.argsort(indices), np.cumsum(counts[:-1]))

    grouped_files = []
    for r in rev:
        grouped_files.append(list(files[r]))
        if len(r) == 1:
            LOGGER.warning("It's a singleton. The file: {}".format(files[r]))

    if return_unique_visits:
        return grouped_files, unique_visits
    else:
        return grouped_files
