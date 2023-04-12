from .region import Region


class AstroRegion(Region):
    def __init__(self, **kwargs):
        Region.__init__(self, **kwargs)


class Ruler(AstroRegion):
    # ALLOWED=('image image','fk5 degrees','fk5 arcmin','fk5 arcsec')
    DELIM = ''

    def __init__(self, x0, y0, x1, y1, ruler='image image', **kwargs):

        self.x0, self, x1 = x0, x1
        self.y0, self.y1 = y0, y1
        self.ruler = ruler
        AstroRegion.__init__(self, **kwargs)

    @property
    def region(self):
        return f'# ruler({self.x0},{self.y0},{self.x1},{self.y1}) ruler={self.ruler}'


class Compass(AstroRegion):
    DELIM = ''

    def __init__(self, x0, y0, compass='fk5', north=True, east=True, **kwargs):
        self.x0 = x0
        self.y0 = y0

        self.compass = compass
        self.N = north
        self.E = east

        AstroRegion.__init__(self, **kwargs)

    @property
    def region(self):
        compass = f'compass={self.compass} {{{N}}} {{{E}}} {int(self.N)} {int(self.E)}'

        return f'# compass({self.x0},{self.y0},{self.size}{self.unit}) {compass}'
