

class Construction44(BlowupConstruction):


	def __init__(self):
	
		self._graph = GraphFlag("8:122886633447755156781122334455667788")

		x = polygen(QQ)
		K = NumberField(x^3 - 2*x^2 + 2*x - 2/3, 'x', embedding=RDF(0.5))
		x = K.gen()

		self._field = K

		self._weights = [x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4]
	
		#self._phantom_edges = [(1, 3), (1, 4)]
		self._phantom_edges = [(1, 4)]



	def k4_density_epsilon(self):

		cn = self._graph.n

		pr = PolynomialRing(self._field, ['epsilon','c1','c2','c3','c4','c5','c6','c7'])
		pr_gens = pr.gens()
		epsilon_weights = [pr_gens[0] * pr_gens[v] for v in range(1, cn)]
		epsilon_weights.append(sum(-x for x in epsilon_weights))

		#pr = PolynomialRing(self._field, ['epsilon','c1','c2','c3','c4','c5','c6','c7','c8'])
		#pr_gens = pr.gens()
		#epsilon_weights = [pr_gens[0] * pr_gens[v] for v in range(1, cn + 1)]

		found = pr(0)
		
		for P in Tuples(range(1, cn + 1), 4):
		
			factor = pr(1)
			for v in P:
				factor *= self._weights[v - 1] + epsilon_weights[v - 1]

			ig = self._graph.degenerate_induced_subgraph(P)
			if ig.ne == 6:
				found += factor
		
		return found


	def subgraph_densities_epsilon(self, n):

		cn = self._graph.n
		sharp_graph_counts = {}
		sharp_graphs = []

		pr = PolynomialRing(self._field, ['epsilon','c1','c2','c3','c4','c5','c6','c7'])
		pr_gens = pr.gens()
		epsilon_weights = [pr_gens[0] * pr_gens[v] for v in range(1, cn)]
		epsilon_weights.append(sum(-x for x in epsilon_weights))
		total = pr(0)
		
		for P in UnorderedTuples(range(1, cn + 1), n):
		
			factor = pr(factorial(n))
			for i in range(1, cn + 1):
				factor /= factorial(P.count(i))
			
			for v in P:
				factor *= self._weights[v - 1] + epsilon_weights[v - 1]

			total += factor
			
			ig = self._graph.degenerate_induced_subgraph(P)
			ig.make_minimal_isomorph()
			
			ghash = hash(ig)
			if ghash in sharp_graph_counts:
				sharp_graph_counts[ghash] += factor
			else:
				sharp_graphs.append(ig)
				sharp_graph_counts[ghash] = factor
		
			
		return [(g, sharp_graph_counts[hash(g)] / total) for g in sharp_graphs]


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

	
	def subgraphs(self, n):
		return [p[0] for p in self.subgraph_densities(n)]


	def zero_eigenvectors(self, tg, flags):

		cn = self._graph.n
		s = tg.n
		k = flags[0].n # assume all flags the same order

		pr = PolynomialRing(self._field, "e")

		rows = []

		if True: # no special edges
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


		if True: # one special edge
			for tv in Tuples(range(1, cn + 1), s):

				one_indices = [i for i in range(s) if tv[i] == 1]
				if len(one_indices) == 0:
					continue
					
				four_indices = [i for i in range(s) if tv[i] == 4]
				if len(four_indices) == 0:
					continue
	
				it = self._graph.degenerate_induced_subgraph(tv)
				
				if it.ne != tg.ne - 1:
					continue
				a = one_indices[0] + 1
				b = four_indices[0] + 1
				it.add_edge((a,b))
				
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
					ig.add_edge((a,b))
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

		if True: # two special edges 1-4, 1-4
			for tv in Tuples(range(1, cn + 1), s):

				one_indices = [i for i in range(s) if tv[i] == 1]
				if len(one_indices) == 0:
					continue
				
				four_indices = [i for i in range(s) if tv[i] == 4]
				if len(four_indices) < 2:
					continue
	
				it = self._graph.degenerate_induced_subgraph(tv)
				
				if it.ne != tg.ne - 2:
					continue
				a = one_indices[0] + 1
				b = four_indices[0] + 1
				c = four_indices[1] + 1
				it.add_edge((a,b))
				it.add_edge((a,c))
				
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
					ig.add_edge((a,b))
					ig.add_edge((a,c))
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


		if True: # two special edges 1-3, 1-4
			for tv in Tuples(range(1, cn + 1), s):

				one_indices = [i for i in range(s) if tv[i] == 1]
				if len(one_indices) == 0:
					continue

				three_indices = [i for i in range(s) if tv[i] == 3]
				if len(three_indices) == 0:
					continue
				
				four_indices = [i for i in range(s) if tv[i] == 4]
				if len(four_indices) == 0:
					continue
	
				it = self._graph.degenerate_induced_subgraph(tv)
				
				if it.ne != tg.ne - 2:
					continue
				a = one_indices[0] + 1
				b = three_indices[0] + 1
				c = four_indices[0] + 1
				it.add_edge((a,b))
				it.add_edge((a,c))
				
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
					ig.add_edge((a,b))
					ig.add_edge((a,c))
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


		if False: # old code for s==6
			for tv in Tuples(range(1, cn + 1), s):
	
				one_indices = [i for i in range(s) if tv[i] == 1]
				
				if len(one_indices) == 0:
					continue
				
				it = self._graph.degenerate_induced_subgraph(tv)
				
				if not it.is_labelled_isomorphic(tg):
					continue
	
				row = [0] * len(flags)
				ov = [4,]
				ig = self._graph.degenerate_induced_subgraph(tv + ov)
				igp = copy(ig)
				igp.add_edge((one_indices[0] + 1, s + 1))
				ig.t = s
				igp.t = s
				ig.make_minimal_isomorph()
				igp.make_minimal_isomorph()
				for j in range(len(flags)):
					if ig.is_labelled_isomorphic(flags[j]):
						row[j] = -1
						break
				for j in range(len(flags)):
					if igp.is_labelled_isomorphic(flags[j]):
						row[j] = 1
						break
						
				rows.append(row)

		if True:
			for tv in Tuples(range(1, cn + 1), s):
	
				one_indices = [i for i in range(s) if tv[i] == 1]
				
				if len(one_indices) == 0:
					continue
				
				a = one_indices[0] + 1
				
				it = self._graph.degenerate_induced_subgraph(tv)
				
				if not it.is_labelled_isomorphic(tg):
					continue
	
				row = [pr(0)] * len(flags)
				epsilon = pr.gen()

				for ov in Tuples(range(1, cn + 2), k - s):
			
					factor = pr(1)
					
					for v in ov:
						if v == cn + 1:
							factor *= epsilon
						elif v == 4:
							factor *= self._weights[4 - 1] - epsilon
						else:
							factor *= self._weights[v - 1]
					
					ov2 = [v if v != cn + 1 else 4 for v in ov]
					ig = self._graph.degenerate_induced_subgraph(tv + ov2)
					for i in range(k - s):
						if ov[i] == cn + 1:
							ig.add_edge((a, s + i + 1))
					ig.t = s
					ig.make_minimal_isomorph()
					
					for j in range(len(flags)):
						if ig.is_labelled_isomorphic(flags[j]):
							row[j] += factor
							break
					else:
						raise ValueError
				
				nrow = [entry[1] for entry in row]
				rows.append(nrow)
			
		#return rows	
		return matrix_of_independent_rows(self._field, rows, len(flags))
