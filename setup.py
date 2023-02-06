from glob import glob
import numpy as np
import os
import requests
from setuptools import setup,Extension



if __name__ == '__main__':

    #filename='latest.tar.gz'
    #url=''
    #path=''
    #r=requsts.get(url,allow_redirects=True)
    #with open(os.path.join(path,filename),'wb') as fp:
    #    fp.write(r.content)

    
    tokens=('slitlessutils','core','wfss','config','polyclip')
    ext=Extension('.'.join(tokens)+'.cpolyclip',
                  glob(os.path.join(*tokens,'src','*.c')),
                  include_dirs=[os.path.join(*tokens,'include'),
                                np.get_include()])
    
    setup(ext_modules=[ext])
