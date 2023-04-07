from .region import Region

class PandaRegion(Region):
    def __init__(self,**kwargs):
        raise NotImplementedError("Panda regions are not implemented.")


class CircularPanda(PandaRegion):
    def __init__(self,*args,**kwargs):
        PandaRegion.__init__(self,**kwargs)

        
class EllipticalPanda(PandaRegion):
    def __init__(self,*args,**kwargs):
        PandaRegion.__init__(self,**kwargs)


class BoxPanda(PandaRegion):
    def __init__(self,*args,**kwargs):
        PandaRegion.__init__(self,**kwargs)
