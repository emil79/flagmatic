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

from sage.structure.sage_object import SageObject
from sage.rings.all import RationalField
from sage.matrix.all import matrix
from copy import copy


def matrix_of_independent_rows(field, rows, width):

	M = matrix(field, rows, sparse=True)
	N = matrix(field, 0, width, sparse=True)
	NE = copy(N)
	
	for i in range(M.nrows()):
		NE2 = NE.stack(M[i, :])
		NE2.echelonize()
		if not NE2[-1, :].is_zero():
			NE = NE2
			N = N.stack(M[i, :])

	return N


class Construction(SageObject):

	def __init__(self):
		self._field = RationalField()

	@property
	def field(self):
		return self._field

	def subgraph_densities(self, n):
		return None

	def zero_eigenvectors(self, tg, flags, flag_basis=None):
		return None
