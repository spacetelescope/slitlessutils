import numpy as np


class SpatialPolynomial(dict):
    ''' Class to implement 2d spatial polynomial '''
    def __init__(self,values):
        self.order=None
        if np.isscalar(values):
            self.order=0
            self[(0,0)]=values
        else:
            # the coefs are stored as cantor pairs.
            # https://en.wikipedia.org/wiki/Pairing_function#Inverting_the_Cantor_pairing_function]

            n=(np.sqrt(1+8*len(values))-1)/2
            if n.is_integer():
                n=int(n)
                self.order=n-1

                # old way of decoding cantor pairing
                i=0
                for j in range(n):
                    for k in range(j+1):
                        self[(j-k,k)]=values[i]
                        i+=1
            else:
                print("ERROR!  input must be an array whose length is a triangular number")


            
    def __str__(self):
        return f'Spatial polynomial of order = {self.order}'

    def evaluate(self,xy):
        ''' evaluate the 2d polynomial at a scalar position (x,y) '''
        return sum(coef*xy[0]**i*xy[1]**j for (i,j),coef in self.items())



if __name__=='__main__':
    values=[0,1,2,3,4,5]
    p=SpatialPolynomial(values)
    if p:
        print('okay')
    print(p)
