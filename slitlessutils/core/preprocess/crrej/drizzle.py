import numpy as np
from pathlib import Path

from astropy.io import fits
from drizzlepac import astrodrizzle

from slitlessutils.core.wfss import WFSSCollection
from slitlessutils import LOGGER


def _get_instrument_defaults(file):
    # Args taken from INS HSTaXe FullFrame Cookbooks
    instrument_args = {"ir": {}, "uvis": {}, "acs": {}}

    with fits.open(file, mode="readonly") as hdul:
        h = hdul[0].header
        telescope = h["TELESCOP"]
        if telescope == "HST":
            instrument = h["INSTRUME"]
            if instrument == "WFC3":
                csmid = h["CSMID"]
                if csmid == "IR":
                    return instrument_args["ir"]
                elif csmid == "UVIS":
                    return instrument_args["uvis"]
                else:
                    raise ValueError(f"Invalid CSMID: {csmid}")
            elif instrument == "ACS":
                if h["DETECTOR"] == "WFC":
                    return instrument_args["acs"]
                elif h["DETECTOR"] == "SBC":
                    raise ValueError("SBC does not need drizzling")
            else:
                raise ValueError(f"'{instrument}' is not supported")
        else:
            raise ValueError(f"'{telescope}' is not supported")


def drizzle(files, outdir=Path().absolute(), **kwargs):
    """
    Runs AstroDrizzle as part of Cosmic Ray handling. Passes any additional arguments
    not listed here directly to AstroDrizzle. Any user-specified arguments will override


    Parameters
    ----------
    files : str or list of str
        The list of files to be processed by AstroDrizzle

    outdir : str or `pathlib.Path`, optional
        A directory to write the final rectified mosaics to.
        By default, the current working directory
    """
    # Start with known defaults
    drizzle_kwargs = {
        "output": "su_drizzle",
        "overwrite": False,
        "preserve": False,
        "clean": True,
    }
    # Apply instrument-specific defaults
    if isinstance(files, list):
        file_to_check = files[0]
    elif isinstance(files, str):
        with open(files, "r") as f:
            file_to_check = f.readline()
    drizzle_kwargs.update(_get_instrument_defaults(file_to_check))
    # Finally override any args with the ones the user supplied
    drizzle_kwargs.update(kwargs)
    # Prepend outdir to output
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    drizzle_kwargs["output"] = str(outdir / drizzle_kwargs["output"])

    return astrodrizzle.AstroDrizzle(files, **drizzle_kwargs)


def group_by_visit(files, return_unique_visits=False, **kwargs):
    """
    Group a list of files based on their unique visits.

    Parameters
    ----------
    files : list of str
        List of filenames

    return_unique_visits : bool, optional
        Indicates whether to return unique visits along with grouped files.
        If `True`, a list containing lists of grouped files and an array of unique
        visits are returned.
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
            LOGGER.warning(f"The file {files[r][0]} was not grouped with any others.")

    if return_unique_visits:
        return grouped_files, unique_visits
    else:
        return grouped_files
