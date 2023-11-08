import sys

from .attributes import Attributes
from .region import Region


class DS9Regions(list):
    ''' Primary class for handling a collection of ds9-style regions. '''

    HEADER = '# Region file format: DS9 version 4.1'

    def __init__(self, filename='ds9.reg', fk5=False, **kwargs):

        # filename
        self.filename = filename

        # units of file
        self.fk5 = fk5

        # set the att
        self.attrs = Attributes(**kwargs)

        # just to keep track
        self.nwritten = 0

        # set the default printing to screen
        self.fp = sys.stdout

    @property
    def header(self):
        pars = str(self.attrs)
        # pars=' '.join([format_value(k,v) for k,v in self.attrs.items()])
        head = f'{self.HEADER}\nglobal {pars}'
        if self.fk5:
            head += "\nfk5"

        return head

    def append(self, reg):
        if isinstance(reg, Region):
            super().append(reg)

    def write(self, filename=None):
        if isinstance(filename, str):
            self.filename = filename
        else:
            filename = self.filename

        # with self as fp:
        #    print(fp)
        #    for reg in self:
        #        self.write_region(reg)

        with open(filename, 'w') as fp:
            self.fp = fp
            print(self.header, file=self.fp)
            for reg in self:
                self.write_region(reg)
        self.fp = sys.stdout

    def write_region(self, reg):
        print(reg.ds9format(attrs=self.attrs, fk5=self.fk5), file=self.fp)
        self.nwritten += 1

    def __str__(self):
        s = f'ds9 region file: {self.filename}'
        return s

    def __enter__(self):
        self.fp = open(self.filename, 'w')
        print(self.header, file=self.fp)
        return self

    def __exit__(self, etype, eval, etb):
        self.fp.close()
        self.fp = sys.stdout

    @classmethod
    def from_arrays(cls, form, *args):
        obj = cls()
        for pars in zip(*args):
            obj.append(ds9region(form, *pars))
        return obj


def ds9region(form, *args, **kwargs):
    ''' A factory function to create a region object '''

    form = form.lower()

    if form == 'point':
        from .onedregion import Point
        return Point(*args, **kwargs)
    elif form == 'text':
        from .onedregion import Text
        return Text(*args, **kwargs)

    elif form == 'circle':
        from .enclosedregion import Circle
        return Circle(*args, **kwargs)
    elif form == 'ellipse':
        from .enclosedregion import Ellipse
        return Ellipse(*args, **kwargs)
    elif form == 'box':
        from .enclosedregion import Box
        return Box(*args, **kwargs)
    elif form == 'polygon':
        from .enclosedregion import Polygon
        return Polygon(*args, **kwargs)

    elif form == 'line':
        from .linearregion import Line
        return Line(*args, **kwargs)
    elif form == 'vector':
        from .linearrregion import Vector
        return Vector(*args, **kwargs)
    elif form == 'projection':
        from .linearrregion import Projection
        return Projection(*args, **kwargs)
    elif form == 'segment':
        from .linearrregion import Segment
        return Segment(*args, **kwargs)

    elif form == 'compass':
        from .astrorregion import Compass
        return Compass(*args, **kwargs)
    elif form == 'ruler':
        from .astrorregion import Ruler
        return Ruler(*args, **kwargs)

    elif form == 'circularannulus':
        from .annularrregion import CircularAnnulus
        return CircularAnnulus(*args, **kwargs)
    elif form == 'ellipticalannulus':
        from .annularrregion import EllipticalAnnulus
        return EllipticalAnnulus(*args, **kwargs)
    elif form == 'boxannulus':
        from .annularrregion import BoxAnnulus
        return BoxAnnulus(*args, **kwargs)

    else:
        raise NotImplementedError(f"The region form {form} is not supported.")


if __name__ == '__main__':

    x = DS9Regions()
    x.append(ds9region('circle', 1, 2, 3))

    print(x)
