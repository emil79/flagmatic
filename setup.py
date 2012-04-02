import os
import sys

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

if not 'SAGE_ROOT' in os.environ:
	print "    ERROR: The environment variable SAGE_ROOT must be defined."
	sys.exit(1)
else:
	SAGE_ROOT  = os.environ['SAGE_ROOT']

setup (

	name='Flagmatic',
	packages=['flagmatic'],
	author='Emil R. Vaughan',
	author_email='e.vaughan@qmul.ac.uk',
	version='2.0',
	cmdclass = {'build_ext': build_ext},
	ext_modules = [
		Extension('flag',
			sources=['flagmatic/flag.pyx'],
			include_dirs = [os.path.join(SAGE_ROOT, 'local/lib/python/site-packages/numpy/core/include'),
				os.path.join(SAGE_ROOT, 'local/include/csage'),
				os.path.join(SAGE_ROOT, 'devel/sage/sage/ext'),
				os.path.join(SAGE_ROOT, 'devel/sage')],
			library_dirs = [os.path.join(SAGE_ROOT, 'local/lib')]
		)
	]
)
