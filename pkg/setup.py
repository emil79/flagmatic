"""

flagmatic 2

Copyright (c) 2012, E. R. Vaughan. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import os
import sys

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

if not 'SAGE_ROOT' in os.environ:
	print "	   ERROR: The environment variable SAGE_ROOT must be defined."
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
		Extension('flagmatic.flag',
			sources=['flagmatic/flag.pyx'],
			include_dirs = [
				os.path.join(SAGE_ROOT, 'local/include'),
				os.path.join(SAGE_ROOT, 'local/include/csage'),
				os.path.join(SAGE_ROOT, 'devel/sage/sage/ext'),
	        		os.path.join(SAGE_ROOT, 'src/sage/ext'),
				os.path.join(SAGE_ROOT, 'devel/sage'),
        			os.path.join(SAGE_ROOT, 'src')],
			library_dirs = [os.path.join(SAGE_ROOT, 'local/lib')],
			extra_compile_args = ["-O3", "-Wall", "-Wno-strict-prototypes"]
		),
		Extension('flagmatic.hypergraph_flag',
			sources=['flagmatic/hypergraph_flag.pyx'],
			include_dirs = [
				os.path.join(SAGE_ROOT, 'local/include'),
				os.path.join(SAGE_ROOT, 'local/lib/python/site-packages/numpy/core/include'),
				os.path.join(SAGE_ROOT, 'local/include/csage'),
				os.path.join(SAGE_ROOT, 'devel/sage/sage/ext'),
        			os.path.join(SAGE_ROOT, 'src/sage/ext'),
				os.path.join(SAGE_ROOT, 'devel/sage'),
        			os.path.join(SAGE_ROOT, 'src')],
			library_dirs = [os.path.join(SAGE_ROOT, 'local/lib')],
			extra_compile_args = ["-O3", "-Wall", "-Wno-strict-prototypes"]
		),
		Extension('flagmatic.three_graph_flag',
			sources=['flagmatic/three_graph_flag.pyx'],
			extra_compile_args = ["-O3", "-Wall", "-Wno-strict-prototypes"]
		),
		Extension('flagmatic.graph_flag',
			sources=['flagmatic/graph_flag.pyx'],
			extra_compile_args = ["-O3", "-Wall", "-Wno-strict-prototypes"]
		),
		Extension('flagmatic.oriented_graph_flag',
			sources=['flagmatic/oriented_graph_flag.pyx'],
			extra_compile_args = ["-O3", "-Wall", "-Wno-strict-prototypes"]
		),
		Extension('flagmatic.multigraph_flag',
			sources=['flagmatic/multigraph_flag.pyx'],
			extra_compile_args = ["-O3", "-Wall", "-Wno-strict-prototypes"]
		)
	]
)
