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

class flagmatic_construction(object):

	def __init__(self):
		pass

	def induced_subgraphs(self, n):
		return ([], [])

	def zero_eigenvectors(self, tg, flags):
		return None


class blowup_construction(flagmatic_construction):

	def __init__(self, g, vertex_transitive=False):

		self._graph = g
		self._vertex_transitive = vertex_transitive
		
	def induced_subgraphs(self, n):

		cn = self._graph[0]
		total = 0
		sharp_graphs = {}
		
		for P in UnorderedTuples(range(1, cn + 1), n):
		
			factor = factorial(n)
			
			for i in range(1, cn + 1):
				factor /= factorial(P.count(i))
		
			ig = induced_subgraph(self._graph, P)
			try:
				sharp_graphs[ig] += factor
			except KeyError:
				sharp_graphs[ig] = factor
			total += factor
		
		min_sharp_graphs = {}
		
		for ig, N in sharp_graphs.iteritems():
			mig = degenerate_minimal_isomorph(ig)
			try:
				min_sharp_graphs[mig] += N
			except KeyError:
				min_sharp_graphs[mig] = N
	
		for g in sorted(min_sharp_graphs.keys(), key = lambda G : len(G[1])):
			density = min_sharp_graphs[g] / Integer(total)
			sys.stdout.write("%s has density %s (%g).\n" % (graph_to_string(g),
				density, density))
	
		return min_sharp_graphs.keys()
	


	def zero_eigenvectors(self, tg, flags, flag_basis=None):
	
		if flag_basis == None:
			flag_basis = identity_matrix(QQ, len(flags))
	
		num_cols = flag_basis.nrows()
		cn = self._graph[0]
		M = matrix(QQ, 0, num_cols, sparse=True)
		
		for tv in Tuples(range(1, cn + 1), tg[0]):
	
			row = matrix(QQ, [asymptotic_flag_density_fixed(self._graph, tg,
					flags[i], tv) for i in range(len(flags))]) * flag_basis.T
			
			# Quick check to see if it is the zero vector!
			if row.is_zero():
				continue
			
			# See if adding row increases the rank
			M = M.stack(matrix(QQ, row))
			if M.rank() < M.nrows():
				M = M[:-1,:]
		
		return M.echelon_form()
