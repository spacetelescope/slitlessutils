import numpy as np

from . import cpolyclip


# NEVER CHANGE THESE
INT=np.int32
FLT=np.float32




def multi(x,y,nx,ny,axis=1):
    
    # compute bounding box
    l=np.clip(np.floor(x.min(axis=axis)-1),0,nx).astype(INT)
    r=np.clip(np.floor(x.max(axis=axis)+1),0,nx).astype(INT)
    b=np.clip(np.floor(y.min(axis=axis)-1),0,ny).astype(INT)
    t=np.clip(np.floor(y.max(axis=axis)+1),0,ny).astype(INT)
    npix=sum((r-l+1)*(t-b+1))
    
    # compute polygon indices.  Here assumign that x is a 2d array,
    # all the pixels must have same number of vertices.  Therefore,
    # this is quickly given as:
    polyinds=np.linspace(0,x.shape[0]*x.shape[1],x.shape[0]+1,dtype=INT)


    # compute polygon indices --- technically, this does't work because
    # at this point x is a numpy array, so they have the same
    # number of polygon edges.
    npoly=len(x)
    polyinds2=np.zeros(npoly+1,dtype=INT)
    for i,xx in enumerate(x):
        polyinds2[i+1]=polyinds[i]+len(xx)
        
    if not np.allclose(polyinds,polyinds2):
        raise RuntimeError
        
        
    # the number of output pixels must be an array (this is a C-gotcha)
    nclip=np.array([0],INT)
    
    # output arrays 
    areas=np.zeros(npix,dtype=FLT)
    xx=np.zeros(npix,dtype=INT)
    yy=np.zeros(npix,dtype=INT)
    
    
    # call the compiled C-code
    cpolyclip.multi(l,r,b,t,x.astype(FLT).ravel(),y.astype(FLT).ravel(),\
                    npoly,polyinds,xx,yy,nclip,areas)
    
    
    # trim the results
    nclip=nclip[0]
    areas=areas[:nclip]
    xx=xx[:nclip]
    yy=yy[:nclip]
    return xx,yy,areas,polyinds




def single(x,y,nx,ny,**kwargs):

    # compute bounding box for the pixel
    l=np.array(np.clip(np.floor(np.min(x)),0,nx),dtype=INT)
    r=np.array(np.clip(np.ceil(np.max(x))+1,0,nx),dtype=INT)
    b=np.array(np.clip(np.floor(np.min(y)),0,ny),dtype=INT)
    t=np.array(np.clip(np.ceil(np.max(y))+1,0,ny),dtype=INT)

    # get number of vertices for the polygon
    nverts=np.array([len(x)],dtype=INT)

    # number of pixels that might be affected
    npix=(r-l+1)*(t-b+1)
    
    # recast some things for C
    nclip=np.array([0],dtype=INT)
    ri=np.zeros(npix+1,dtype=INT)
    
    # output polygon indices
    px_out=np.zeros((nverts[0]+24)*npix,dtype=FLT)
    py_out=np.zeros((nverts[0]+24)*npix,dtype=FLT)
    
    # main outputs (area, pixel coords and reverse indices)
    areas=np.zeros(npix,dtype=FLT)
    inds=np.zeros((npix,2),dtype=INT)
    ri_out=np.zeros(npix+1,dtype=INT)
    
    # call the pologyon clipper
    cpolyclip.single(l,r,b,t,np.array(x,dtype=FLT),np.array(y,dtype=FLT),\
                     nverts,px_out,py_out,inds,nclip,areas,ri_out)
    
    # extract data
    nclip=nclip[0]
    px_out=px_out[:nclip]
    py_out=py_out[:nclip]
    ri_out=ri_out[:nclip]
    
    # main outputs 
    xx=inds[:nclip,0]
    yy=inds[:nclip,1]
    areas=areas[:nclip]
    
    return xx,yy,areas




def polyclip(x,y,nx,ny,**kwargs):
    print('inside polyclip')

    
