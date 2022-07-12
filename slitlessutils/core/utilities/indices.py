import numpy as np

def compress(indices):
    ''' compress an integer index to lowst possible integers

    They will be defined such that:

    indices = unique_indices[compressed_indices]

    '''


    unique_indices=np.unique(indices)
    compressed_indices=np.digitize(indices,unique_indices)-1
    return compressed_indices,unique_indices



def uniq(indices):
    ''' return the integers from an input array of integers '''

    idx = np.unique(indices, return_index=True)[1]
    unique_indices=indices[np.sort(idx)]
    return unique_indices



def reverse(ints,ignore=()):
    ''' find the reverse indices as a dictionary

    1) first find all unique values
    2) return dict whose keys are the unique values and elements are 
       the indices from input list
    '''

    uniq,ind,cnt=np.unique(ints,return_inverse=True,return_counts=True)
    rev=np.split(np.argsort(ind),np.cumsum(cnt[:-1]))

    # changed to this so ints can be a list
    ri={u:np.unravel_index(p,np.shape(ints)) for u,p in zip(uniq,rev) if u not in ignore}

    return ri

def decimate(val,*indices,dims=None,unravel=True,return_factor=False):
    ''' sum a vector `v` over repeated tuples.

    val is an array
    indices is a list of integer arrays of equal length

    therefore the tuples that are preserved are the paired sets:

    tuples = list(zip(*indices))

    returns the summed vector and unique indices

    '''
    # number of dimensions
    ndim=len(indices)

    # we need to convert a tuple into a one-dimensional.
    # could be done by hand (see snippet below) or with ravelling
    # If we don't pass dimensions, then grab that from the max value
    # of the dimension.  Passing dimensions will be faster.
    if dims is None:
        dims=[np.amax(index)+1 for index in indices]
    idx=np.ravel_multi_index(indices,dims,order='F')

    # find the unique indices and unique inverses
    out,uind,cinv=np.unique(idx,return_index=True,return_inverse=True)

    # sum the values over the compressed index
    vv=np.bincount(cinv,weights=val)

  
    
    # get the unique indices as another tuple of arrays
    if unravel:
        out=tuple(index[uind] for index in indices)
        ret=(vv,*out)
    else:
        ret=(vv,out,dims)


    if return_factor:
        factor=float(len(val))/float(len(vv))
        ret+=(factor,)
    return ret

        




if __name__=='__main__':

    x=np.array([0,0,1,1,1,2,3,3])
    y=np.array([1,1,0,3,2,2,3,4])
    l=np.array([0,0,1,1,6,6,3,4])
    ri=reverse(l)

    for ll,g in ri.items():
        print(x[g].shape,x[g[0]].shape)
    lkdjf



    i=np.array([[1,1,1,1,2],
                [2,2,4,4,9],
                [9,8,2,3,1]],dtype=int)
    ri=reverse(i)


    print(ri)

    lkdjf

    x=np.array([1,1,2,2,2,2,3],dtype=np.uint16)
    y=np.array([1,1,2,2,2,2,3],dtype=np.uint16)
    l=np.array([1,1,2,2,2,2,3],dtype=np.uint16)
    v=np.array([1,1,2,2,2,2,3],dtype=np.float64)


    #x=np.array([],dtype=np.uint16)
    #y=np.array([],dtype=np.uint16)
    #l=np.array([],dtype=np.uint16)
    #v=np.array([],dtype=np.uint16)
    dims=(10,10,10)
    vv,xx,yy,ll=decimate(v,x,y,l,dims=dims)
    #print(vv,xx,yy,ll)

    vv,xx=decimate(v,x)
    print(vv,xx)
    print(v,x)


    segids=[[1,1,1,1,1],
            [2,2,0,0,5],
            [2,2,0,0,4]]
    ri=reverse(segids)
    print(ri)


    ldj


















#if __name__=='__main__':
#    print('hi')
#
###
#
#    x=np.array([[0,0,1,1,0],
#                [1,0,2,2,2],
#                [0,5,5,10,0]],dtype=int)
#
#
#
#    #x=np.ravel(x)
#    r=reverse(x)
#    print(x.shape)
#    #print(r)
#    for ind,pix in r.items():
#        #pix=np.unravel_index(p,x.shape)
#        print(pix)
#
#        print(x[pix])
#        #print(ind,x[pix],pix)
