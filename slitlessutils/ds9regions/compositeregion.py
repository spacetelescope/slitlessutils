from .region import Region


class Composite(Region, list):
    def __init__(self, **kwargs):
        Region.__init__(self, **kwargs)

    @property
    def region(self):

        x = 0.0
        y = 0.0
        r = []

        for reg in self:
            x += reg.x
            y += reg.y
            # print(reg.region)
            # r.append(reg.region)
