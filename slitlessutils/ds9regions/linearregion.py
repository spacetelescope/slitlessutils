from .region import Region


class LinearRegion(Region):
    def __init__(self, **kwargs):
        Region.__init__(self, **kwargs)


class Line(LinearRegion):
    DELIM = ''

    def __init__(self, x0, y0, x1, y1, left=False, right=False, **kwargs):
        self.x0, self, x1 = x0, x1
        self.y0, self.y1 = y0, y1
        self.arrow = {'left': left, 'right': right}
        LinearRegion.__init__(self, **kwargs)

    @property
    def region(self):
        l, r = int(self.arrow['left']), int(self.arrow['right'])
        return f'line({self.x0},{self.y0},{self.x1},{self.y1}) # line={l} {r}'


class Vector(LinearRegion):
    def __init__(self, x0, y0, length, angle, **kwargs):
        self.x0 = x0
        self.y0 = y0
        self.length = length
        self.angle = angle
        LinearRegion.__init__(self, **kwargs)

    @property
    def region(self):
        return f'# vector({self.x0},{self.y0},{self.length}{self.unit},{self.angle}) vector=1'


class Projection(LinearRegion):
    def __init__(self, x0, y0, length, angle, **kwargs):
        self.x0 = x0
        self.y0 = y0
        self.length = length
        self.angle = angle
        LinearRegion.__init__(self, **kwargs)

    @property
    def region(self):
        return f'# projection({self.x0},{self.y0},{self.length},{self.angle})'


class Segment(LinearRegion):
    def __init__(self, x0, y0, length, angle, **kwargs):
        self.x0 = x0
        self.y0 = y0
        self.length = length
        self.angle = angle
        LinearRegion.__init__(self, **kwargs)

    @property
    def region(self):
        xy = ','.join([f'{x},{y}' for x, y in zip(self.x, self.y)])
        return f'# segment({xy})'
