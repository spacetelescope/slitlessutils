import numpy as np


def as_iterable(d):
    """
    Method to take a variable and recast it as an interable.

    Issue here is that sometimes a user may pass a variable like
    orders = '+1', which is technically an iterable.  But we need to be
    careful and distinguish between orders=('+1','+2','+3') or so on.
    This will recast the first case as: orders=('+1',), so when it
    is iterated upon like:

    for order in orders:
       do_something()

    then things work as expected.

    Parameters
    ----------
    d : all types
       an iterable version

    """

    if d is None:
        return ()
    else:
        if isinstance(d, (str, int, float)):
            return (d,)
        elif isinstance(d, (list, np.ndarray)):
            return tuple(d)
        elif isinstance(d, tuple):
            return d
        else:
            raise TypeError('unsupported data type')


if __name__ == '__main__':
    print(mkiter('+1'))
    print(mkiter(('+1', '+2')))
    print(mkiter(['+1', '+2']))
    print(mkiter(np.array(('+1', '+2'))))
    print(mkiter(1))
