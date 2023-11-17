import numpy as np

"""
A set of methods to work with indices, which are `np.ndarray`s of
integer-like data.

"""


def compress(indices):
    """
    Method to compress an integer index to the lowest possible set of
    integers.

    Parameters
    ----------
    indices : `np.ndarray`
       The array of indices to compress

    Returns
    -------
    compind : `np.ndarray`
       The lowest set of integers that maintain the uniqueness of the input
       array.  This will always have dtype=int. See examples below

    uniqind : `np.ndarray`
       The unique possible values of the input array.  This will have the
       same dtype as the input `indices`.  See examples below

    Examples
    --------
    >>> import numpy as np
    >>> from slitlessutils.core.utilities import indices
    >>> x=np.array([4,6,4,8,9,2,10,100,5], dtype=int)
    >>> ycomp, yuniq = indices.compress(x)
    >>> ycomp  # doctest: +IGNORE_OUTPUT
    array([1, 3, 1, 4, 5, 0, 6, 7, 2])
    >>> yuniq  # doctest: +IGNORE_OUTPUT
    array([  2,   4,   5,   6,   8,   9,  10, 100])

    """

    unique_indices = np.unique(indices)
    compressed_indices = np.digitize(indices, unique_indices) - 1
    return compressed_indices, unique_indices


def uniq(indices):
    """
    Method to return only the unique indices

    Parameters
    ----------
    indices : `np.ndarray`
        indices to find unique values of

    Returns
    -------
    uniqind : `np.ndarray`
        the unique indices

    Notes
    -----
    The order of the indices is preserved


    Examples
    --------
    >>> import numpy as np
    >>> from slitlessutils.core.utilities import indices
    >>> x=np.array([4,6,4,8,9,2,10,100,5], dtype=int)
    >>> u=indices.uniq(x)
    >>> u
    array([  4,   6,   8,   9,   2,  10, 100,   5])

    """

    idx = np.unique(indices, return_index=True)[1]
    unique_indices = indices[np.sort(idx)]
    return unique_indices


def reverse(ints, ignore=()):
    """
    Method to obtain the reverse indices

    Parameters
    ----------
    ints : `np.ndarray`
       the integers to find the reverse of

    ignore : tuple
       A set of indices to ignore when computing reverse

    Returns
    -------
    revind : dict
       The reversed indices:
       The keywords will be the unique indices
       The values will be the indices in the input array of that value

    Examples
    --------
    >>> import numpy as np
    >>> from slitlessutils.core.utilities import indices
    >>> x=np.array([4,6,4,8,9,2,10,100,5], dtype=int)
    >>> ri=indices.reverse(x)
    >>> for i,j in ri.items():
    ...     print(i,j)  # doctest: +IGNORE_OUTPUT
    2 (array([5]),)
    4 (array([0, 2]),)
    5 (array([8]),)
    6 (array([1]),)
    8 (array([3]),)
    9 (array([4]),)
    10 (array([6]),)
    100 (array([7]),)

    Notice here, that the first column are the unique values of `x`,
    whereas the second column are the indices of `x` that have the
    first column as a value.  These are tuples of arrays, as `x` may be a
    multi-dimensional array.  If one computed the length of the second
    column, then this would be effectively be a histogram (or the same
    as what is returned by `np.bincount()`)

    """

    uniq, ind, cnt = np.unique(ints, return_inverse=True, return_counts=True)
    rev = np.split(np.argsort(ind), np.cumsum(cnt[:-1]))

    # changed to this so ints can be a list
    ri = {u: np.unravel_index(r, np.shape(ints)) for u, r in zip(uniq, rev) if u not in ignore}

    return ri


def decimate(val, *indices, dims=None, unravel=True, return_factor=False):
    """
    Method to decimate a vector over the repeated indices, where
    'decimation' is defined as summing the values if they have the
    same indices.

    Parameters
    ----------
    val : `np.ndarray`
        The values to decimate.  May be of any datatype.

    indices : tuple
        A tuple of `np.ndarray`, where each `np.ndarray` must be of an
        integer dtype.

    dims : tuple, list, `np.ndarray` or None, optional
        The maximum possible value of all the input indices.  If this is
        `None`, then the *1+np.amax()* value of each of the `indices` tuple
        is computed.  If this is an iterable type (tuple, list, `np.ndarray`)
        then it should have the same length as `*indices`.  Default is `None`.

    unravel : bool, optional
        Flag to unravel the output indices such that they match the
        dimensionality of the input indices.  Internally, the multi-
        dimensional indices are grouped into a single one-dimensional
        index, which is used for all calculations.  For this reason, the
        dims must be used.  If unravel=False, then the internally-used
        one-dimensional index is returned, otherwise, the decimated
        multi-dimensional indices are returned.  Default is True

    return_factor : bool, optional
        Flag to return the compression factor achieved by decimating.
        Default is False

    Returns
    -------
    decval : `np.ndarray`
        The decimated array.  Will be same dtype as `val`.

    outind : tuple
        The decimated indices. Will be the same length as `*indices` and
        each entry will be of same dtype as `*indices`.

    factor : float
        The compression factor: the ratio of the sizes of the decimated
        arrays to the input arrays.  This is only computed and returned if
        `return_factor = True`

    Examples
    --------
    We will create a dataset of (x,y,l) triplets, with some duplicates, and
    an array of values v, which will be decimated.

    >>> import numpy as np
    >>> from slitlessutils.core.utilities import indices
    >>> x=np.array([1,1,2,2,3,3,1,1,3,3,4],dtype=np.uint16)
    >>> y=np.array([1,1,2,2,2,2,1,1,3,4,5],dtype=np.uint16)
    >>> l=np.array([1,2,2,2,3,2,1,1,3,3,6],dtype=np.uint16)
    >>> v=np.array([2,4,6,8,7,5,3,1,8,6,4],dtype=np.float64)

    Notice here that the triplets (x,y,l)=(2,2,2) and (1,1,1) appear
    two and three times, respectively.  Therefore, the decimated
    values should have those triplets appearing only once, but have their
    values of (6,8) and (2,3,1) summed (ie. 14 and 6), respectively.
    Whereas all triplets that appear a single time will be just their
    given value:

    >>> vv,xx,yy,ll = indices.decimate(v,x,y,l)
    >>> vv,xx,yy,ll  # doctest: +IGNORE_OUTPUT
    array([ 6.,  4., 14.,  5.,  7.,  8.,  6.,  4.])
    array([1, 1, 2, 3, 3, 3, 3, 4], dtype=uint16)
    array([1, 1, 2, 2, 2, 3, 4, 5], dtype=uint16)
    array([1, 2, 2, 2, 3, 3, 3, 6], dtype=uint16)

    Furthermore, it is important to notice, we did not specify the dimensions
    of any of the integer axes, and so the `decimate()` computed them
    as the *one plus max* of each dimension: so dim=(5,6,7).  But we could
    set them if the values are known (for example, if x,y represent
    pixel coordinates in a WFSS image, then they'd be shape of the image).
    Importantly, if one dimension is set, then *ALL* dimensions must be set

    >>> vv,xx,yy,ll = indices.decimate(v,x,y,l,dims=(5,6,7))

    Finally, the number of dimensions is arbitrary:

    >>> vv,xx,yy = indices.decimate(v,x,y,dims=(5,6))
    >>> vv,xx,yy  # doctest: +IGNORE_OUTPUT
    array([10., 14., 12.,  8.,  6.,  4.])
    array([1, 2, 3, 3, 3, 4], dtype=uint16)
    array([1, 2, 2, 3, 4, 5], dtype=uint16)

    But here notice that the coordinate pair (x,y)=(1,1) appears *FOUR*
    times and so the values and lengths of the decimated arrays are different


    Notes
    -----
    This is a fundamental piece to the creation of matrices for linear
    extraction techniques.  Please see for more details on the need
    for this function.

    """

    # we need to convert a tuple into a one-dimensional.
    # could be done by hand (see snippet below) or with ravelling
    # If we don't pass dimensions, then grab that from the max value
    # of the dimension.  Passing dimensions will be faster.
    if dims is None:
        dims = [np.amax(index) + 1 for index in indices]
    idx = np.ravel_multi_index(indices, dims, order='F')

    # find the unique indices and unique inverses
    out, uind, cinv = np.unique(idx, return_index=True, return_inverse=True)

    # sum the values over the compressed index
    vv = np.bincount(cinv, weights=val)

    # get the unique indices as another tuple of arrays
    if unravel:
        out = tuple(index[uind] for index in indices)
        ret = (vv, *out)
    else:
        ret = (vv, out, dims)

    if return_factor:
        factor = float(len(val)) / float(len(vv))
        ret += (factor,)
    return ret


def span(val, *indices, dims=None, unravel=True):
    if dims is None:
        dims = [np.amax(index) + 1 for index in indices]
    idx = np.ravel_multi_index(indices, dims, order='F')

    out, uind, cnt = np.unique(idx, return_counts=True, return_inverse=True)
    rev = np.split(np.argsort(uind), np.cumsum(cnt[:-1]))

    v = np.zeros_like(indices[0], dtype=type(val))
    for i, r in enumerate(rev):
        v[i] = np.amax(val[r]) - np.amin(val[r])

    if unravel:
        out = tuple(index[uind] for index in indices)
        ret = (v, *out)
    else:
        ret = (v, out, dims)

    return ret


if __name__ == '__main__':

    x = np.array([1, 2, 2, 2, 2, 2, 3], dtype=np.uint16)
    y = np.array([1, 1, 2, 3, 2, 2, 3], dtype=np.uint16)
    l = np.arange(len(x), dtype=int)

    # l = np.array([1, 1, 2, 2, 2, 2, 3], dtype=np.uint16)
    # v = np.array([1, 1, 2, 2, 2, 2, 3], dtype=np.float64)

    x = np.array([0, 0, 1, 1, 1, 2, 3, 3])
    y = np.array([1, 1, 0, 3, 2, 2, 3, 4])
    l = np.array([0, 0, 1, 1, 6, 6, 3, 4])
    ri = reverse(l)

    for ll, g in ri.items():
        print(x[g].shape, x[g[0]].shape)

    i = np.array([[1, 1, 1, 1, 2],
                  [2, 2, 4, 4, 9],
                  [9, 8, 2, 3, 1]], dtype=int)
    ri = reverse(i)

    print(ri)

    x = np.array([1, 1, 2, 2, 2, 2, 3], dtype=np.uint16)
    y = np.array([1, 1, 2, 2, 2, 2, 3], dtype=np.uint16)
    l = np.array([1, 1, 2, 2, 2, 2, 3], dtype=np.uint16)
    v = np.array([1, 1, 2, 2, 2, 2, 3], dtype=np.float64)

    # x=np.array([],dtype=np.uint16)
    # y=np.array([],dtype=np.uint16)
    # l=np.array([],dtype=np.uint16)
    # v=np.array([],dtype=np.uint16)
    dims = (10, 10, 10)
    vv, xx, yy, ll = decimate(v, x, y, l, dims=dims)
    # print(vv,xx,yy,ll)

    vv, xx = decimate(v, x)
    print(vv, xx)
    print(v, x)

    segids = [[1, 1, 1, 1, 1],
              [2, 2, 0, 0, 5],
              [2, 2, 0, 0, 4]]
    ri = reverse(segids)
    print(ri)
