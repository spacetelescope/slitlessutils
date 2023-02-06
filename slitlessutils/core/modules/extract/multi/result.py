import numpy as np

from ....utilities import headers


class Result:
    """
    Class to hold the results of a matrix inversion

    Parameters
    ----------
    method : str
       The function in `scipy.sparse.linalg` for inversion
    
    x : `np.ndarray`
       The array of spectral values, conceptually given by x = A^-1.b

    istop : int
       An integer indicating the stopping criterion, see 
       `scipy.sparse.linalg.lsqr()` or `scipy.sparse.linalg.lsmr()`

    itn : int
       The number of iterations

    r1norm : float
       The norm of A^-1.b-x (ie. excluding the damping)

    r2norm : float
       The norm including damping, which is given by:
       sqrt( norm(r)^2  +  damp^2 * norm(x - x0)^2 )
    
    anorm : float
       Estimated Frobenius norm of Abar = [[A]; [damp*I]]

    acond : float
       Estimate of the condition of Abar

    arnorm : float
        Estimate of ``norm(AT.r - damp^2*(x - x0))``.

    xnorm : float
        The norm of the results

    unc : `np.ndarray`
        The estimated uncertainties.

    damp : float
        The damping value       

    Notes
    -----
    This is unlikely to be directly initialized, but often used.    
    """    
    
    
    def __init__(self,method,x,istop,itn,r1norm,r2norm,anorm,acond,
                 arnorm,xnorm,unc,damp):
        self.method=method
        self.x=x
        self.istop=istop
        self.itn=itn
        self.r1norm=r1norm
        self.r2norm=r2norm
        self.anorm=anorm
        self.acond=acond
        self.arnorm=arnorm
        self.xnorm=xnorm
        self.lo=unc.copy()
        self.hi=unc.copy()
        self.damp=damp

    

    @property
    def logdamp(self):
        if self.damp==0:
            return np.NINF
        else:
            return np.log10(self.damp)

    @property
    def xy(self):
        return np.log10(self.r1norm),np.log10(self.xnorm)

    @property
    def chi2(self):
        return self.r1norm


    
    def update_header(self,hdr):
        """
        Method to update a fits header
        
        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update

        """
        

        hdr['METHOD']=(self.method,'method for inversion')
        hdr['ACOND']=(self.acond,'condition of the regularized matrix')
        hdr['ANORM']=(self.anorm,'norm of the regularized matrix')
        hdr['ARNORM']=(self.arnorm,'norm of the regularized matrix')
        hdr['R1NORM']=(self.r1norm,'standard chi2')
        hdr['R2NORM']=(self.r2norm,'regularized chi2')
        hdr['XNORM']=(self.xnorm,'norm of the result vector')
        hdr['ISTOP']=(self.istop,'stopping condition')
        hdr['ITER']=(self.itn,'number of iterations')
        hdr['DAMPING']=(self.damp,'regularization parameter')
        headers.add_stanza(hdr,'Inversion Results',before='METHOD')


    def __len__(self):
        return len(self.x)

    def __str__(self):
        return f'result for log(\u2113)={self.damp}'
        
    def __add__(self,dx):
        """
        Method to implement the addition to the optimal values
        
        Parameters
        ----------
        dx : `np.ndarray`
            The value to add
        
        Returns
        -------
        self an updated self object        
        """

        self.x+=dx
        return self

    def __radd__(self,dx):
        return self.__add__(dx)


    def __iadd__(self,dx):
        return self.__add__(dx)

        
