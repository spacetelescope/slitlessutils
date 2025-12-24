from importlib.metadata import metadata


def get_metadata(package=None):
    '''
    Function to fetch the metadata of a package.

    Returns
    -------
    meta : `dict`
       A dictionary containing the package metadata
    '''

    if package is None:
        package = __package__.split('.')
        package = package[0]

    meta = metadata(package)
    if meta:
        return meta
    else:
        raise ModuleNotFoundError(f'Invalid package: {package}')
