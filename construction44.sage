

class Construction44(BlowupConstruction):


	def __init__(self):
	
		self._graph = GraphFlag("8:122886633447755156781122334455667788")

		x = polygen(QQ)
		K = NumberField(x^3 - 2*x^2 + 2*x - 2/3, 'x', embedding=RDF(0.5))
		x = K.gen()

		self._field = K

		self._weights = [x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4]
	
		self._phantom_edges = [(1, 3), (1, 4)]
	

	def subgraph_densities(self, n):

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
			igc = copy(ig) # copy for phantom edge
			ig.make_minimal_isomorph()
			
			ghash = hash(ig)
			if ghash in sharp_graph_counts:
				sharp_graph_counts[ghash] += factor
			else:
				sharp_graphs.append(ig)
				sharp_graph_counts[ghash] = factor

			total += factor
		
			for pes in Subsets(self._phantom_edges, 1):
				if all(all(x in P for x in pe) for pe in pes):
					ig = copy(igc)
					for pe in pes:
						ig.add_edge([P.index(x) + 1 for x in pe])
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


