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

def ClebschGraph():

	cleb = Flag("0:",2)
	cleb.n = 16
	edges = [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 7), (2, 8), (2, 9), (2, 10),
		(3, 7), (3, 11), (3, 12), (3, 13), (4, 8), (4, 11), (4, 14), (4,15), (5, 9),
		(5, 12), (5, 14), (5, 16), (6, 10), (6, 13), (6, 15), (6, 16), (7, 14), (7, 15),
		(7, 16), (8, 12), (8, 13), (8, 16), (9, 11), (9, 13), (9, 15), (10, 11), (10, 12),
		(10, 14), (11, 16), (12, 15), (13, 14)]
	for e in edges:
		cleb.add_edge(e)
	return cleb


class SymmetricBlowupConstruction (BlowupConstruction):


	def tuple_orbits(self, k):
		
		SG = self._graph.Graph()
		G, d = SG.automorphism_group(translation=True)

		# Sage gives the graph new labels! Get a translation dictionary.
		rd = dict((v,k) for (k,v) in d.iteritems())

		gen_str = ", ".join(str(t) for t in G.gens())
		gap_str = "g := Group(%s);" % gen_str
		gap.eval(gap_str)

		orbs = gap.new("Orbits(g, Tuples([1..%d], %d), OnTuples);" % (self._graph.n, k)).sage()

		total = 0
		orb_reprs = {}
		for o in orbs:
			sys.stdout.write("Orbit %d:\n" % len(orb_reprs))
			orb_reprs[tuple(o[0])] = len(o)
			for t in o:
				ig = self._graph.degenerate_induced_subgraph(map(lambda x : rd[x], t))
				ig.make_minimal_isomorph()
				sys.stdout.write("%s " % ig)
			sys.stdout.write("\n")

		return orb_reprs
		

		

	# NOTE: This computes orbits on *sets* and then expands these sets in different
	# ways to form k-tuples. This results in more representatives than is strictly
	# necessary, but it is much faster than doing otherwise.

	def tuple_orbit_reps(self, k, prefix=[]):
		
		s = len(prefix)
		tp = tuple(prefix)
		if s > k:
			raise ValueError
		
		SG = self._graph.Graph()
		G, d = SG.automorphism_group(translation=True)

		# Sage gives the graph new labels! Get a translation dictionary, and
		# relabel the generators back to how they should be.

		rd = dict((v,k) for (k,v) in d.iteritems())
		trans_gens = [gen.cycle_tuples() for gen in G.gens()]
		gens = [tuple(tuple(map(lambda x : rd[x], cy)) for cy in gen) for gen in trans_gens]

		# Pass generators to GAP to create a group for us.
		
		gen_str = ",".join("(" + "".join(str(cy) for cy in cys) + ")" for cys in gens)
		gap.eval("g := Group(%s);" % gen_str)
		if len(prefix) > 0:
			gap.eval("g := Stabilizer(g, %s, OnTuples);" % list(set(prefix)))

		S = []
		for i in range(1, k - s + 1):
			S.extend([tuple(sorted(list(x))) for x in Subsets(self._graph.n, i)])
		
		set_orb_reps = {}

		#sys.stdout.write("Calculating orbits")

		while len(S) > 0:

			rep = list(S[0])

			o = gap.new("Orbit(g, %s, OnSets);" % (rep,)).sage()
			o = list(set([tuple(sorted(t)) for t in o]))
			ot = o[0]
			set_orb_reps[ot] = len(o)
			for t in o:
				S.remove(t)
			#sys.stdout.write(".")
			#sys.stdout.flush()

		#sys.stdout.write("\n")

		combs = [tuple(c) for c in Compositions(k - s)]
		factors = []
		for c in combs:
			factor = factorial(k - s)
			for x in c:
				factor /= factorial(x)
			factors.append(factor)

		orb_reps = {}
		total = 0
		
		for ot, length in set_orb_reps.iteritems():

			ne = len(ot)
			for ci in range(len(combs)):
				c = combs[ci]
				if len(c) == ne:
					t = tp
					for i in range(ne):
						t += c[i] * (ot[i],)
					weight = factors[ci] * length
					orb_reps[t] = weight
					total += weight

		return (total, orb_reps)
	

	def zero_eigenvectors(self, tg, flags, flag_basis=None):

		s = tg.n
		k = flags[0].n # assume all flags the same order

		rows = []

		t_total, t_orb_reps = self.tuple_orbit_reps(s)

		for t_rep, t_factor in t_orb_reps.iteritems():
			
			for tp in Permutations(t_rep):

				it = self._graph.degenerate_induced_subgraph(tp)
				if not it.is_equal(tg):
					continue

				total, orb_reps = self.tuple_orbit_reps(k, prefix=tp)
				
				row = [0] * len(flags)
				
				for P, factor in orb_reps.iteritems():
				
					ig = self._graph.degenerate_induced_subgraph(P)
					ig.t = s
					ig.make_minimal_isomorph()
					
					for j in range(len(flags)):
						if ig.is_equal(flags[j]):
							row[j] += Integer(factor) / total
							break
				
				rows.append(row)

		if flag_basis == None:
			flag_basis = identity_matrix(QQ, len(flags), sparse=True)

		if len(rows) == 0:
			return matrix(QQ, 0, flag_basis.nrows(), sparse=True)

		M = matrix(QQ, rows, sparse=True) * flag_basis.T
		
		if M.rank() == 0:
			return matrix(QQ, 0, flag_basis.nrows(), sparse=True)
		
		M = M.echelon_form()
		M = M[:M.rank(),:]

		return M
	
	
	def subgraph_densities(self, n):

		sharp_graph_counts = {}
		sharp_graphs = []
		
		total, orb_reps = self.tuple_orbit_reps(n)
		
		sys.stdout.write("Found %d orbits.\n" % len(orb_reps))
		
		for P, factor in orb_reps.iteritems():
		
			ig = self._graph.degenerate_induced_subgraph(P)
			ig.make_minimal_isomorph()
			
			ghash = hash(ig)
			if ghash in sharp_graph_counts:
				sharp_graph_counts[ghash] += factor
			else:
				sharp_graphs.append(ig)
				sharp_graph_counts[ghash] = factor

		return [(g, sharp_graph_counts[hash(g)] / Integer(total)) for g in sharp_graphs]
	