import numpy as np

from .spatialpolynomial import SpatialPolynomial


class ParametricPolynomial(list):
    ''' Class to implement the parametric orders '''

    def __init__(self,name):
        self.name=name
        self.order=-1
        self.invert = lambda xy,f: None    # provide a fallback

    def append(self,coefs):
        ''' add another order to the list '''

        poly=SpatialPolynomial(coefs)
        if poly:
            super().append(poly)
            self.order+=1
        
            if self.order==1:
                self.invert=self._first
            elif self.order==2:
                self.invert=self._second
            else:
                self.invert=self._nth

    def _first(self,xy,f):
        ''' analytically invert a 1st order polynomial '''
        coefs=self.coefs(xy)
        return (f-coefs[0])/coefs[1]

    def _second(self,xy,f):
        ''' analytically invert a 2nd order polynomial '''
        coefs=self.coefs(xy)


        disc=coefs[1]*coefs[1]-4*(coefs[0]-f)*coefs[2]
        sqrt=np.sqrt(disc)
        nump=(-coefs[1]+sqrt)/(2*coefs[2])
        #numm=(-coefs[1]-sqrt)/(2*coefs[2])

        return nump


    def _nth(self,xy,f):
        ''' numerically invert an arbitrary order polynomial '''
        raise NotImplementedError('arbitrary order not implemented.')

    
    def coefs(self,xy):#,order=None):
        ''' compute the polynomial coeffs at some position (x,y) '''
        #if order is None:
        coefs=[poly.evaluate(xy) for poly in self]
        #else:
        #    coefs=self[order].evaluate(xy)
        return coefs

    def evaluate(self,xy,t):
        ''' evaluate the polynomial at some position (x,y) and parameter '''
        return sum(p.evaluate(xy)*t**i for i,p in enumerate(self))

    def deriv(self,xy,t):
        ''' compute the derivative of the polynomial '''
        return sum(p.evaluate(xy)*i*t**(i-1) for i,p in enumerate(self[1:],start=1))
