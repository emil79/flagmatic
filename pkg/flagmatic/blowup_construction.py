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

from sage.rings.arith import factorial
from sage.combinat.all import UnorderedTuples, Tuples, Combinations, Permutations, Compositions, Subsets
from sage.rings.all import Integer, RationalField
from sage.interfaces.gap import gap
from copy import copy

from three_graph_flag import *
from graph_flag import *
from oriented_graph_flag import *
from construction import *


class BlowupConstruction(Construction):


	def __init__(self, g, weights=None, field=None, phantom_edge=None, no_symmetry=False):
	
		if g.oriented and g.is_degenerate:
			raise NotImplementedError("degenerate oriented graphs not supported.")
	
		self._graph = copy(g)

		if weights is None:
			self._weights = None
		else:
			if len(weights) != g.n:
				raise ValueError
			self._weights = weights
	
		if field is None:
			self._field = RationalField()
		else:
			self._field = field

		if not phantom_edge is None:
			# check edge is valid; will get an Exception if not.
			h = copy(g)
			h.add_edge(phantom_edge)
			self._phantom_edge = phantom_edge

		# Only make use of symmetry when all the conditions are right...
		# Should probably allow OrientedGraphFlag
	
		if (field is None and weights is None and g.n > 4 and not no_symmetry
			and (type(g) is GraphFlag or type(g) is OrientedGraphFlag)
			and phantom_edge is None):
			self._use_symmetry = True
		else:
			self._use_symmetry = False


	@property
	def graph(self):
		return self._graph

		
	@property
	def weights(self):
		return self._weights


	@property
	def field(self):
		return self._field
	

	def subgraph_densities(self, n):

		if self._use_symmetry:
			return self.symm_subgraph_densities(n)

		cn = self._graph.n
		total = Integer(0)
		sharp_graph_counts = {}
		sharp_graphs = []

		for P in UnorderedTuples(range(1, cn + 1), n):
		
			factor = factorial(n)
			for i in range(1, cn + 1):
				factor /= factorial(P.count(i))
			
			if self._weights:
				for v in P:
					factor *= self._weights[v - 1]
			
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


	def zero_eigenvectors(self, tg, flags):

		if self._use_symmetry:
			return self.symm_zero_eigenvectors(tg, flags)

		cn = self._graph.n
		s = tg.n
		k = flags[0].n # assume all flags the same order

		rows = []

		for tv in Tuples(range(1, cn + 1), s):

			it = self._graph.degenerate_induced_subgraph(tv)
			
			using_phantom_edge = False

			if hasattr(self, "_phantom_edge") and it.ne == tg.ne - 1:
				extra_edges = [e for e in tg if not e in it]
				if len(extra_edges) == 1:
					phantom_edge = extra_edges[0]
					if all(tv[phantom_edge[i] - 1] == self._phantom_edge[i] for i in range(tg.r)):
						it.add_edge(phantom_edge)
						using_phantom_edge = True
		
			if not (using_phantom_edge or it.is_labelled_isomorphic(tg)):
				continue

			total = Integer(0)
			row = [0] * len(flags)
		
			for ov in UnorderedTuples(range(1, cn + 1), k - s):
		
				factor = factorial(k - s)
				for i in range(1, cn + 1):
					factor /= factorial(ov.count(i))

				if self._weights:
					for v in ov:
						factor *= self._weights[v - 1]
				
				ig = self._graph.degenerate_induced_subgraph(tv + ov)
				if using_phantom_edge:
					ig.add_edge(phantom_edge)
				ig.t = s
				ig.make_minimal_isomorph()
				
				for j in range(len(flags)):
					if ig.is_labelled_isomorphic(flags[j]):
						row[j] += factor
						total += factor
						break
						
			for j in range(len(flags)):
				row[j] /= total	
			rows.append(row)

		return matrix_of_independent_rows(self._field, rows, len(flags))


	#
	# "Symmetric" versions follow.
	#
		

	# NOTE: This computes orbits on *sets* and then expands these sets in different
	# ways to form k-tuples. This results in more representatives than is strictly
	# necessary, but it is much faster than doing otherwise.

	def tuple_orbit_reps(self, k, prefix=[]):
		
		s = len(prefix)
		tp = tuple(prefix)
		if s > k:
			raise ValueError

		gens = self._graph.automorphism_group_gens()
		
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
	

	def symm_zero_eigenvectors(self, tg, flags, flag_basis=None):

		s = tg.n
		k = flags[0].n # assume all flags the same order

		rows = []

		t_total, t_orb_reps = self.tuple_orbit_reps(s)

		for t_rep, t_factor in t_orb_reps.iteritems():
			
			for tp in Permutations(t_rep):

				it = self._graph.degenerate_induced_subgraph(tp)
				if not it.is_labelled_isomorphic(tg):
					continue

				total, orb_reps = self.tuple_orbit_reps(k, prefix=tp)
				
				row = [0] * len(flags)
				
				for P, factor in orb_reps.iteritems():
				
					ig = self._graph.degenerate_induced_subgraph(P)
					ig.t = s
					ig.make_minimal_isomorph()
					
					for j in range(len(flags)):
						if ig.is_labelled_isomorphic(flags[j]):
							row[j] += Integer(factor) / total
							break
				
				rows.append(row)

		return matrix_of_independent_rows(self._field, rows, len(flags))
	
	
	def symm_subgraph_densities(self, n):

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
