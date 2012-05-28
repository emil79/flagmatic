

class Construction17(BlowupConstruction):


	def __init__(self):
	
		self._graph = ThreeGraphFlag("3:112223331123")
		self._field = RationalField()
		self._weights = None
		
		self._phantom_edges = []
	

	def subgraph_densities(self, n):

		cn = self._graph.n
		total = Integer(0)
		sharp_graph_counts = {}
		sharp_graphs = []

		for P in Tuples(range(1, cn + 1), n):
		
			factor = 1
			
			ig = self._graph.degenerate_induced_subgraph(P)
			igc = copy(ig) # copy for phantom edge
			ig.make_minimal_isomorph()
			
			ghash = hash(ig)
			if ghash in sharp_graph_counts:
				sharp_graph_counts[ghash] += factor
			else:
				sharp_graphs.append(ig)
				sharp_graph_counts[ghash] = factor

			total += factor
		
			if P.count(1) >= 2:
				ig = igc
				one_indices = [i for i in range(n) if P[i] == 1]
				two_indices = [i for i in range(n) if P[i] == 2]
				a1 = one_indices[0] + 1
				b1 = one_indices[1] + 1
				for x in two_indices:
					ig.delete_edge((a1, b1, x + 1))
				for x in one_indices[2:]:
					ig.add_edge((a1, b1, x + 1))
				ig.make_minimal_isomorph()
				ghash = hash(ig)
				if not ghash in sharp_graph_counts:
					sharp_graphs.append(ig)
					sharp_graph_counts[ghash] = Integer(0)
				
		return [(g, sharp_graph_counts[hash(g)] / total) for g in sharp_graphs]


	def zero_eigenvectors(self, tg, flags):

		cn = self._graph.n
		s = tg.n
		k = flags[0].n # assume all flags the same order

		rows = []

		for tv in Tuples(range(1, cn + 1), s):

			it = self._graph.degenerate_induced_subgraph(tv)
			
			if not it.is_labelled_isomorphic(tg):
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

		for tv in Tuples(range(1, cn + 1), s):

			if tv.count(1) < 2:
				continue
			
			it = self._graph.degenerate_induced_subgraph(tv)

			one_indices = [i for i in range(s) if tv[i] == 1]
			two_indices = [i for i in range(s) if tv[i] == 2]
			a1 = one_indices[0] + 1
			b1 = one_indices[1] + 1
			for x in two_indices:
				it.delete_edge((a1, b1, x + 1))
			for x in one_indices[2:]:
				it.add_edge((a1, b1, x + 1))
			
			if not it.is_labelled_isomorphic(tg):
				continue

			print tv, it

			total = Integer(0)
			row = [0] * len(flags)
		
			for ov in Tuples(range(1, cn + 1), k - s):
		
				factor = 1
				
				ig = self._graph.degenerate_induced_subgraph(tv + ov)
				ext_one_indices = one_indices + [i + s for i in range(k - s) if ov[i] == 1]
				ext_two_indices = two_indices + [i + s for i in range(k - s) if ov[i] == 2]
				for x in ext_two_indices:
					ig.delete_edge((a1, b1, x + 1))
				for x in ext_one_indices[2:]:
					ig.add_edge((a1, b1, x + 1))
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
