import subprocess
import os

from ...logger import LOGGER
        
def gzip(filename):
    r=subprocess.call(['gzip','-f',filename])
    if r != 0:
        LOGGER.warning(f'Cannot gzip {filename} with error: {r}')
    return filename+'.gz'

def gunzip(filename):
    r=subprocess.call(['gunzip',filename])
    if r != 0:
        LOGGER.warning(f'Cannot gunzip {filename} with error: {r}')


    newfilename,exten=os.path.splitext(filename)

    return newfilename



if __name__=='__main__':
    orig='/Users/rryan/t.t'
    os.system(f'touch {orig}')
    tmp=gzip(orig)
    new=gunzip(tmp)
    print(tmp,new)
    print(new==orig)
    os.system(f'rm {orig}')
