import itertools

def generate_graphs(n):

	def degrees(g):
		return [len([e for e in g[1] if x in e]) for x in range(1, g[0] + 1)] 

	if n == 0:
	
		return [(0, ())]

	max_ne = (n - 1) * (n - 2) / 2
	max_e = n * max_ne / 3
	
	new_graphs = [set() for i in range(max_e + 1)]
	
	smaller_graphs = generate_graphs(n - 1)
	
	possible_nbrs = [p for p in itertools.combinations(range(1, n), 2)]

	for sg in smaller_graphs:
	
		pe = len(sg[1])
		ds = degrees(sg)
		maxd = max(ds + [0])
				
		for ne in range(maxd, max_ne + 1):
		
			for nb in itertools.combinations(possible_nbrs, ne):
			
				ng = (n, sg[1] + tuple([(v[0], v[1], n) for v in nb]))
				mng = minimal_isomorph(ng)
				new_graphs[pe + ne].add(mng)

	#print new_graphs

	return [g for graphs in new_graphs for g in graphs]

