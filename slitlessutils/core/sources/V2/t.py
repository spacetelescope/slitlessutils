class Region:
    def __init__(self,x,y,w):
        self.x=x
        self.y=y
        self.w=w
    def __iter__(self):
        yield from zip(self.x,self.y,self.w)
        
class Source:
    def __init__(self):
        self.regions=[]

        self.regions.append(Region([1,2,3],[2,3,4],[.1,.1,.2]))
        self.regions.append(Region([11,31],[22,42],[.2,.4]))


    def __iter__(self):
        for region in self.regions:
            yield from region

if __name__=='__main__':
    s=Source()

    for x,y,w in s:
        print(x,y,w)
        

        
