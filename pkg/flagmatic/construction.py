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

import sys

from sage.structure.sage_object import SageObject
from sage.rings.all import Integer, QQ, polygen, RationalField
from sage.matrix.all import matrix
from sage.misc.all import sage_eval

from three_graph_flag import *
from graph_flag import *
from oriented_graph_flag import *


def matrix_of_independent_rows(field, rows, width):

	M = matrix(field, rows, sparse=True)

	if len(rows) == 0 or M.is_zero():
		return matrix(field, 0, width, sparse=True)

	N = M[0, :]
	NE = N
	for i in range(1, M.nrows()):
		NE2 = NE.stack(M[i, :])
		NE2.echelonize()
		if not NE2[-1,:].is_zero():
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

	def target_bound(self):
		return None

	def zero_eigenvectors(self, tg, flags, flag_basis=None):
		return None


class AdHocConstruction(Construction):

	def __init__(self, name):

		self._name = name

		if name == "maxs3":
	
			self._field_str = "NumberField(x^2+x-1/2, 'x', embedding=0.5)"
	
		elif name == "maxs4":

			self._field_str = "NumberField(x^3+x^2+x-1/3, 'x', embedding=0.5)"

		elif name == "maxs5":

			self._field_str = "NumberField(x^4+x^3+x^2+x-1/4, 'x', embedding=0.5)"

		else:
		
			self._field_str = "RationalField()"
		
		x = polygen(QQ)
		self._field = sage_eval(self._field_str, locals={'x': x})


	@property
	def field(self):
		return self._field
	

	def target_bound(self):

		K = self._field
		x = K.gen()

		if self._name == "maxs3":

			return (4 * x - 1, [0, 2, 4, 5])

		elif self._name == "maxs4":
		
			return (1 - 9 * x**2, [0, 5, 8, 24, 27, 31, 37, 38])

		elif self._name == "max42":

			return (Integer(3)/4, [0, 4, 8, 23, 24, 27, 33])


	def zero_eigenvectors(self, tg, flags):
	
		x = self._field.gen()
		rows = []

		# TODO: should check flags are in the right order!
		
		if self._name == "maxs3":

			if tg == OrientedGraphFlag("1:"):
		
				rows = [
					[1 - x,       0,     x],
					[x * (1 - x), 1 - x, x**2]
				]
	
		elif self._name == "maxs4":

			if tg == OrientedGraphFlag("2:"):
		
				rows = [
					[1 - x,       0, 0, 0, 0, 0,   0, 0, x],
					[x * (1 - x), 0, 0, 0, 0, 1-x, 0, 0, x**2]
				]
				
			elif tg == OrientedGraphFlag("2:12"):

				rows = [
					[0, 1 - x, 0, 0, 0, 0, 0, 0, x],
					[0, 0,     0, 0, 1, 0, 0, 0, -1],
					[0, 0,     0, 0, 0, 1, 0, 0, 0],
					[0, 0,     0, 0, 0, 0, 1, 0, -1]
				]

		elif self._name == "max42":

			if tg == ThreeGraphFlag("3:"):

				rows = [[1, 0, 0, 0, 1, 1, 1, 0]]
			
			elif tg == ThreeGraphFlag("3:123"):

				rows = [[0, 1, 1, 1, 0, 0, 0, 1]]


		return matrix_of_independent_rows(self._field, rows, len(flags))
