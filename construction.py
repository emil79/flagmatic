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

class Construction(SageObject):

	def __init__(self):
		pass

	def induced_subgraphs(self, n):
		return ([], [])

	def zero_eigenvectors(self, tg, flags):
		return None

	def edge_density(self):
		return 0
		
	def subgraph_density(self, h):
		return 0


class BlowupConstruction(Construction):

	def __init__(self, g):
	
		self._graph = g
	
	
	def edge_density(self):
	
		return self._graph.degenerate_edge_density()
	
		
	def subgraph_density(self, h):
	
		return self._graph.degenerate_subgraph_density(h)
	
		
	def induced_subgraphs(self, n):

		cn = self._graph.n
		total = 0
		sharp_graph_counts = {}
		sharp_graphs = []
		
		for P in UnorderedTuples(range(1, cn + 1), n):
		
			factor = factorial(n)
			
			for i in range(1, cn + 1):
				factor /= factorial(P.count(i))
		
# 		for P in Tuples(range(1, cn + 1), n):
# 			factor = 1
		
			ig = self._graph.degenerate_induced_subgraph(P)
			ig.make_minimal_isomorph()
			
			ghash = hash(ig)
			if ghash in sharp_graph_counts:
				sharp_graph_counts[ghash] += factor
			else:
				sharp_graphs.append(ig)
				sharp_graph_counts[ghash] = factor

			total += factor
		
		sys.stdout.write("The following %d graphs appear in the construction:\n" %
			len(sharp_graphs))
		
		for gs in sorted(sharp_graphs, key = lambda g : g.ne):
			density = sharp_graph_counts[hash(gs)] / Integer(total)
			sys.stdout.write("%s has density %s (%g).\n" % (gs,
				density, density))
	
		return sharp_graphs


	def zero_eigenvectors(self, tg, flags, flag_basis=None):
	
		if flag_basis == None:
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
		

	def tuple_orbit_reps(self, k):
		
		SG = self._graph.Graph()
		G, d = SG.automorphism_group(translation=True)

		# Sage gives the graph new labels! Get a translation dictionary.
		rd = dict((v,k) for (k,v) in d.iteritems())

		gen_str = ", ".join(str(t) for t in G.gens())
		gap_str = "g := Group(%s);" % gen_str
		gap.eval(gap_str)

		T = [tuple(t) for t in UnorderedTuples(range(1, self._graph.n + 1), k)]
		total = 0
		orb_reps = {}

		sys.stdout.write("Calculating orbits:\n")

		while len(T) > 0:

			rep = T[0]

			factor = factorial(k)
			for i in range(1, self._graph.n + 1):
				factor /= factorial(rep.count(i))

			o = gap.new("Orbit(g, %s, OnTuples);" % (list(rep),)).sage()
			o = list(set(tuple(sorted(t)) for t in o))

			ot = tuple(map(lambda x : rd[x], o[0]))
			orb_reps[ot] = len(o) * factor
			total += len(o) * factor
			for t in o:
				#ig = self._graph.degenerate_induced_subgraph(ot)
				#ig.make_minimal_isomorph()
				#sys.stdout.write("%s " % ig)
				T.remove(t)
			sys.stdout.write(".")
			sys.stdout.flush()

		return (total, orb_reps)
		

# 	def tuple_orbit_reps(self, k):
# 		
# 		SG = self._graph.Graph()
# 		G, d = SG.automorphism_group(translation=True)
# 
# 		# Sage gives the graph new labels! Get a translation dictionary.
# 		rd = dict((v,k) for (k,v) in d.iteritems())
# 
# 		gen_str = ", ".join(str(t) for t in G.gens())
# 		gap_str = "g := Group(%s);" % gen_str
# 		gap.eval(gap_str)
# 
# 		orbs = gap.new("Orbits(g, Tuples([1..%d], %d), OnTuples);" % (self._graph.n, k)).sage()
# 
# 		total = 0
# 		orb_reps = {}
# 		for o in orbs:
# 			rep = tuple(map(lambda x : rd[x], o[0]))
# 			orb_reps[rep] = len(o)
# 			total += len(o)
# 	
# 		return (total, orb_reps)
		
	
	
	def induced_subgraphs(self, n):

		sharp_graph_counts = {}
		sharp_graphs = []
		
		total, orb_reps = self.tuple_orbit_reps(n)
		
		sys.stdout.write("Found %d orbits.\n" % len(orb_reps))
		
		for P, factor in orb_reps.iteritems():
		
			ig = self._graph.degenerate_induced_subgraph(P)
			ig.make_minimal_isomorph()
			
			#print P, ig, factor
			
			ghash = hash(ig)
			if ghash in sharp_graph_counts:
				sharp_graph_counts[ghash] += factor
			else:
				sharp_graphs.append(ig)
				sharp_graph_counts[ghash] = factor

		sys.stdout.write("The following %d graphs appear in the construction:\n" %
			len(sharp_graphs))
		
		for gs in sorted(sharp_graphs, key = lambda g : g.ne):
			density = sharp_graph_counts[hash(gs)] / Integer(total)
			sys.stdout.write("%s has density %s (%g).\n" % (gs,
				density, density))
	
		return sharp_graphs			
		