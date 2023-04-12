from .region import Region


class AnnularRegion(Region):
    def __init__(self, **kwargs):
        Region.__init__(self, **kwargs)


class CircularAnnulus(AnnularRegion):
    def __init__(self, x0, y0, radii, **kwargs):
        self.x0 = x0
        self.y0 = y0

        self.radii = radii
        AnnularRegion.__init__(self, **kwargs)

    @property
    def region(self):
        unit = self.unit
        radii = [f'{radius}{unit}' for radius in self.radii]
        radii = ','.join(radii)
        return f'annulus({self.x0},{self.y0},{radii})'


class EllipticalAnnulus(AnnularRegion):
    def __init__(self, x0, y0, majors, minors, angle, **kwargs):
        self.x0 = x0
        self.y0 = y0
        self.majors = majors
        self.minors = minors
        self.angle = angle
        AnnularRegion.__init__(self, **kwargs)

    @property
    def region(self):

        unit = self.unit
        tmp = [f'{a}{unit},{b}{unit}' for a, b in zip(self.majors, self.minors)]
        tmp = ','.join(tmp)

        return f'ellipse({self.x0},{self.y0},{tmp},{self.angle})'


class BoxAnnulus(AnnularRegion):
    def __init__(self, x0, y0, dx, dy, angle, **kwargs):
        self.x0 = x0
        self.y0 = y0
        self.dx = dx
        self.dy = dy
        self.angle = angle

        AnnularRegion.__init__(self, **kwargs)

    @property
    def region(self):

        unit = self.unit
        tmp = [f'{dx}{unit},{dy}{unit}' for dx, dy in zip(self.dx, self.dy)]
        tmp = ','.join(tmp)
        return f'box({self.x0},{self.y0},{tmp},{self.angle})'
