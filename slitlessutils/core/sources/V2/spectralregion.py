



class SpectralRegion:
    def __init__(self,x,y,w):
        
        self.x=x
        self.y=y
        self.w=w
        self.sed=SED()

    def __len__(self):
        return len(x)

    def __iter__(self):
        yield from zip(self.x,self.y,self.w)

    
