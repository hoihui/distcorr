import os
import numpy
from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

cmdclass = { }
ext_modules = [ ]

if use_cython:
    ext_modules += [
        Extension("distcorr", [ "distcorr.pyx" ]),
    ]
    cmdclass.update({ 'build_ext': build_ext })
else:
    ext_modules += [
        Extension("distcorr", [ "distcorr.c" ]),
    ]

setup(
    name="distcorr",
    version='1',
    description='Simple random number generators',
    author='Hoi Hui',
    author_email='hoiyinhui@gmail.com',
    url='https://github.com/hoihui/distcorr',
    cmdclass = cmdclass,
    ext_modules=ext_modules,    
    include_dirs=[numpy.get_include(),
                  os.path.join(numpy.get_include(), 'numpy')]
)
