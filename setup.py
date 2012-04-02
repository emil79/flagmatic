import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

if not os.environ.has_key('SAGE_ROOT'):
    print "    ERROR: The environment variable SAGE_ROOT must be defined."
    sys.exit(1)
else:
    SAGE_ROOT  = os.environ['SAGE_ROOT']
    SAGE_LOCAL = SAGE_ROOT + '/local'
    SAGE_DEVEL = SAGE_ROOT + '/devel'
    SAGE_INC = SAGE_LOCAL + '/include/'

setup(
	name='Flagmatic',
	packages=['flagmatic'],
	author='Emil R. Vaughan',
	author_email='e.vaughan@qmul.ac.uk',
	version='2.0',
	cmdclass = {'build_ext': build_ext},
	ext_modules = [
		Extension('flagmatic.flag', sources=['flagmatic/flag.pyx'],
		include_dirs = [SAGE_LOCAL + '/lib/python/site-packages/numpy/core/include',
			SAGE_LOCAL + '/include/csage',
			SAGE_DEVEL + '/sage/sage/ext',
			SAGE_DEVEL + '/sage'],
		library_dirs = [SAGE_LOCAL + '/lib'])
	]
)
