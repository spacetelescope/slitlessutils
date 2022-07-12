import numpy as np
from setuptools import setup,Extension
import os
from glob import glob

if __name__ == '__main__':

    tokens=('slitlessutils','core','instrument','XML','polyclip')
    ext=Extension('.'.join(tokens)+'.cpolyclip',
                  glob(os.path.join(*tokens,'src','*.c')),
                  include_dirs=[os.path.join(*tokens,'include'),
                                np.get_include()])
    
    setup(ext_modules=[ext])
