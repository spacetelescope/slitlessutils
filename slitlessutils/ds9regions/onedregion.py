from .region import Region

class OneDRegion(Region):
    def __init__(self,**kwargs):
        Region.__init__(self,**kwargs)



class Point(OneDRegion):
    ALLOWED=('circle','box','diamond','cross','x','arrow','boxcircle')
    DELIM=''
    def __init__(self,x,y,point=None,**kwargs):
        self.x=x
        self.y=y
        if point not in self.ALLOWED:
            point='x'
        self.point=point
            
        OneDRegion.__init__(self,**kwargs)

    @property
    def region(self):
        return f'point({self.x},{self.y}) # point={self.point}'



class Text(OneDRegion):
    DELIM=''
    def __init__(self,x,y,text,**kwargs):
        self.x=x
        self.y=y
        kwargs['text']=text
        OneDRegion.__init__(self,**kwargs)
        
    @property
    def region(self):
        return f'# {self.include}text({self.x},{self.y})'
