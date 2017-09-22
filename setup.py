import os
import numpy
from distutils.core import setup
from distutils.extension import Extension

setup(
      name="distcorr",
      ext_modules=[Extension("distcorr",["distcorr.c"])],
      include_dirs=[numpy.get_include(),
                    os.path.join(numpy.get_include(), 'numpy')]
)
