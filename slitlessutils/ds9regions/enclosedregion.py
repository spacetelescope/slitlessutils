from .region import Region


class EnclosedRegion(Region):
    def __init__(self, **kwargs):
        Region.__init__(self, **kwargs)


class Polygon(EnclosedRegion):
    def __init__(self, x, y, **kwargs):
        self.x = x
        self.y = y
        Region.__init__(self, **kwargs)

    @property
    def region(self):
        xy = ','.join([f'{x},{y}' for x, y in zip(self.x, self.y)])
        return f'polygon({xy})'


class Circle(EnclosedRegion):
    def __init__(self, x, y, r, **kwargs):
        self.x = x
        self.y = y
        self.r = r

        EnclosedRegion.__init__(self, **kwargs)

    @property
    def region(self):
        return f'circle({self.x},{self.y},{self.r}{self.unit})'


class Ellipse(EnclosedRegion):
    def __init__(self, x, y, a, b, t, **kwargs):
        self.x = x
        self.y = y
        self.a = a
        self.b = b
        self.t = t

        EnclosedRegion.__init__(self, **kwargs)

    @property
    def region(self):
        return f'ellipse({self.x},{self.y},{self.a}{self.unit},{self.b}{self.unit},{self.t})'


class Box(EnclosedRegion):
    def __init__(self, x, y, dx, dy, t, **kwargs):
        self.x = x
        self.y = y
        self.dx = dx
        self.dy = dy
        self.t = t
        EnclosedRegion.__init__(self, **kwargs)

    @property
    def region(self):
        return f'box({self.x},{self.y},{self.dx}{self.unit},{self.dy}{self.fk5},{self.t})'
