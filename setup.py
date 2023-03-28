from glob import glob
import numpy as np
import os
from setuptools import setup, Extension


if __name__ == '__main__':

    tokens=('slitlessutils','core','wfss','config','polyclip')
    ext=Extension('.'.join(tokens)+'.cpolyclip',
                  glob(os.path.join(*tokens,'src','*.c')),
                  include_dirs=[os.path.join(*tokens,'include'),
                                np.get_include()])

    setup(ext_modules=[ext])
