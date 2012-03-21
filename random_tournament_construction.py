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

from sage.all import Integer, QQ, matrix, factorial, identity_matrix, binomial, Combinations
from sage.structure.sage_object import SageObject


class RandomTournamentConstruction(Construction):

		
	def induced_subgraphs(self, n):

		tg = Flag()

		graphs = self.induced_flags(n, tg, [])

		for pair in graphs:
			g, den = pair
			sys.stdout.write("%s has density %s (%g).\n" % (g, den, den))

		return [pair[0] for pair in graphs]
		

	def zero_eigenvectors(self, tg, flags, flag_basis=None):
		
		rows = set()
		for p in Tuples([0, 1], binomial(tg.n, 2)):
			edges = []
			c = 0
			for i in range(1, tg.n + 1):
				for j in range(1, i):
					if p[c] == 0:
						edges.append((i, j))
					else:	
						edges.append((j, i))
					c += 1
			graphs = self.induced_flags(flags[0].n, tg, edges)
			row = [0 for f in flags]
			for pair in graphs:
				g, den = pair
				for i in range(len(flags)):
					if g.is_equal(flags[i]):
						row[i] = den
						break
			rows.add(tuple(row))

		if flag_basis == None:
			flag_basis = identity_matrix(QQ, len(flags), sparse=True)
		
		if len(rows) == 0:
			return matrix(K, 0, flag_basis.nrows(), sparse=True)

		M = matrix(QQ, list(rows), sparse=True) * flag_basis.T

		if M.rank() == 0:
			return matrix(QQ, 0, flag_basis.nrows(), sparse=True)
		
		M = M.echelon_form()
		M = M[:M.rank(),:]

		return M


	def induced_flags(self, n, tg, type_edges):
	
		flag_counts = {}
		flags = []
		total = 0
		
		for p in Tuples([0, 1], binomial(n, 2) - binomial(tg.n, 2)):
			
			edges = list(type_edges)
			
			c = 0
			for i in range(tg.n + 1, n + 1):
				for j in range(1, i):
					if p[c] == 0:
						edges.append((i, j))
					else:	
						edges.append((j, i))
					c += 1

			ig = Flag()
			ig.n = n
			ig.t = tg.n
			
			for s in Combinations(range(1, n + 1), 3):
				if ((s[0], s[1]) in edges and (s[1], s[2]) in edges and (s[2], s[0]) in edges) or (
					(s[0], s[2]) in edges and (s[2], s[1]) in edges and (s[1], s[0]) in edges):
					ig.add_edge(s)
			
			it = ig.induced_subgraph(range(1, tg.n + 1))
			if tg.is_equal(it):
				ig.make_minimal_isomorph()
				
				ghash = hash(ig)
				if ghash in flag_counts:
					flag_counts[ghash] += 1
				else:
					flags.append(ig)
					flag_counts[ghash] = 1
	
			total += 1
		
		return [(f, flag_counts[hash(f)] / Integer(total)) for f in flags]
		