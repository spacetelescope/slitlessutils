import tqdm
import multiprocessing as mp
import psutil as ps
from functools import partial

from ...logger import LOGGER



class Pool:
    def __init__(self,func,ncpu=None,desc=None,quiet=False):
        ''' instantiate the pool '''

        # get some settings for the processing
        ncpus=ps.cpu_count(logical=False)
        ncores=ps.cpu_count(logical=True)
        nthreads=ncores//ncpus
        nmax=ncpus-nthreads
        
        
        # set a default to the max
        if ncpu is None or ncpu<=0:
            self.ncpu=nmax
        else:
            self.ncpu=min(max(ncpu,1),nmax)    # force this to be in range
        self.desc=desc
        self.func=func
        self.quiet=quiet

        

    def __zip__(self,itrs,*args):
        ''' internal generator to zip iterables to scalars '''
        for itr in itrs:
            yield (itr,*args)
        
    def __worker__(self,args):
        ''' internal method to unpack arguments '''
        return self.func(*args)


    #def __enter__(self):
    #    return self


    #def __exit__(self,etype,eval,etb):
    #    pass

    

    def __str__(self):
        lines=['Pool object with:',
               'NCPU = {}'.format(self.ncpu),
               'FUNC = {}'.format(self.func)]
        return '\n'.join(lines)

    #@suppress_warnings
    def __call__(self,itrs,*args,total=None,**kwargs):
        ''' call the pool '''

        # number of iterations to do
        if total is None:
            total=len(itrs)    

        # get the number of CPUs to use
        ncpu=min(total,self.ncpu)
                    
        # start multiprocessing as necessary
        if ncpu==1:
            LOGGER.info('Serial processing')
            results=[self.func(i,*args,**kwargs) for i in
                     tqdm.tqdm(itrs,total=total,desc=self.desc)]
            
        else:
            LOGGER.info(f'Parallel processing: {total} jobs with {ncpu} processes')
            if kwargs is not None:
                func=self.func    # save it for later
                self.func=partial(self.func,**kwargs)


            
            # actually start the processing pool:                
            with mp.Pool(processes=ncpu) as p:
                imap=p.imap(self.__worker__,self.__zip__(itrs,*args))

                results=list(tqdm.tqdm(imap,total=total,desc=self.desc))

            #p=mp.Pool(processes=self.ncpu)
            #imap=p.imap(self.__worker__,self.__zip__(itrs,*args))
            #results=list(tqdm.tqdm(imap,total=total,desc=self.desc))
                        
            if kwargs is not None:
                self.func=func
            
        return results
