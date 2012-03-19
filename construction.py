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

from sage.all import Integer, QQ, matrix, factorial, identity_matrix
from sage.structure.sage_object import SageObject


class Construction(SageObject):

	def __init__(self):
		pass

	def edge_density(self):
		return 0

	def subgraph_densities(self, n):
		return []

	def zero_eigenvectors(self, tg, flags):
		return None


class BlowupConstruction(Construction):

	def __init__(self, g):
	
		if g.oriented:
			raise NotImplementedException("oriented graphs not supported.")
	
		self._graph = g
	
	
	def edge_density(self):
	
		den_pairs = self.subgraph_densities(self._graph.r)

		for pair in den_pairs:
			g, den = pair
			if g.ne == 1:
				return den		
	
		
	def subgraph_densities(self, n):

		cn = self._graph.n
		total = Integer(0)
		sharp_graph_counts = {}
		sharp_graphs = []

# An "inefficient" alternative (sometimes it is faster...UnorderedTuples might have overhead)
# 		for P in Tuples(range(1, cn + 1), n):
# 			factor = 1
		
		for P in UnorderedTuples(range(1, cn + 1), n):
		
			factor = factorial(n)
			
			for i in range(1, cn + 1):
				factor /= factorial(P.count(i))
				
			ig = self._graph.degenerate_induced_subgraph(P)
			ig.make_minimal_isomorph()
			
			ghash = hash(ig)
			if ghash in sharp_graph_counts:
				sharp_graph_counts[ghash] += factor
			else:
				sharp_graphs.append(ig)
				sharp_graph_counts[ghash] = factor

			total += factor
		
		return [(g, sharp_graph_counts[hash(g)] / total) for g in sharp_graphs]


	def zero_eigenvectors(self, tg, flags, flag_basis=None):
	
		if flag_basis is None:
			flag_basis = identity_matrix(QQ, len(flags))
	
		rows = set()
		for tv in Tuples(range(1, self._graph.n + 1), tg.n):
			rows.add(tuple(self._graph.degenerate_flag_density(tg, flags, tv)))
	
		M = matrix(QQ, list(rows), sparse=True) * flag_basis.T
		
		if M.rank() == 0:
			return matrix(QQ, 0, flag_basis.nrows(), sparse=True)
		
		M = M.echelon_form()
		M = M[:M.rank(),:]

		return M



