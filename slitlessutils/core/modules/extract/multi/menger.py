import numpy as np


def menger(xyj, xyk, xyl):
    """
    Function to compute the Menger curvature from three (x,y) pairs

    Parameters
    ----------
    xyj : 2-tuple
       The lower pair of points (x,y)

    xyk : 2-tuple
       The central pair of points (x,y)

    xyl : 2-tuple
       The upper pair of points  (x,y)

    Returns
    -------
    curv : float
       The curvature at the central pair of points

    Notes
    -----
    1) The Menger curvature is a measure of the local curvature that is equal
       to the reciprocal of the radius of a circle that contains the
       current point and adjacent points.

    2) see https://en.wikipedia.org/wiki/Menger_curvature


    """
    # xyj=np.array(xyj)
    # xyk=np.array(xyk)
    # xyl=np.array(xyl)

    if np.allclose(xyj, xyk) or np.allclose(xyk, xyl) or np.allclose(xyj, xyl):
        # if any of the points are the same, then the curvature should be zero
        curv = 0.
    else:

        # could maybe consider subtracting off xyk from xyj and xyl to
        # improve the precision?

        num = 2. * np.abs(xyj[0] * (xyk[1] - xyl[1])
                          + xyk[0] * (xyl[1] - xyj[1])
                          + xyl[0] * (xyj[1] - xyk[1]))

        djk = np.hypot(xyk[0] - xyj[0], xyk[1] - xyj[1])
        dkl = np.hypot(xyl[0] - xyk[0], xyl[1] - xyk[1])
        dlj = np.hypot(xyj[0] - xyl[0], xyj[1] - xyl[1])
        den = djk * dkl * dlj

        curv = num / den

    return curv


if __name__ == '__main__':

    j = (200., 300.)
    k = (400., 500.)
    l = (600., 500.)

    print(menger(j, k, l))

    j = (j[0] - k[0], j[1] - k[1])
    l = (l[0] - k[0], l[1] - k[1])
    k = (0, 0)

    print(menger(j, k, l))
