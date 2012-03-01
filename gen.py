import itertools

def degrees(g):

	n = g[0]
	edges = g[1]
	return [len([e for e in edges if x in e]) for x in range(1, n + 1)]

def generate_flags(n, tg, forbidden_edge_numbers={}):

	check_forbidden_edge_numbers = len(forbidden_edge_numbers) > 0

	t = tg[0]

	if n < t:
		return []

	if n == t:
		return [tg]

	max_ne = (n - 1) * (n - 2) / 2
	max_e = n * max_ne / 3
	
	new_graphs = [set() for i in range(max_e + 1)]
	
	smaller_graphs = generate_flags(n - 1, tg, forbidden_edge_numbers=forbidden_edge_numbers)
	
	possible_nbrs = [p for p in itertools.combinations(range(1, n), 2)]

	for sg in smaller_graphs:
	
		pe = len(sg[1])
		ds = degrees(sg)
		maxd = max(ds[t:] + [0])
			
		for ne in range(maxd, max_ne + 1):
		
			for nb in itertools.combinations(possible_nbrs, ne):
			
				ng = (n, sg[1] + tuple([(v[0], v[1], n) for v in nb]))

				if check_forbidden_edge_numbers:
					if has_forbidden_edge_numbers(ng, forbidden_edge_numbers):
						continue

				mng = minimal_isomorph(ng, tg)
				new_graphs[pe + ne].add(mng)

	return [g for graphs in new_graphs for g in graphs]


def generate_graphs(n, forbidden_edge_numbers={}):
	
	t = (0, ())
	return generate_flags(n, t, forbidden_edge_numbers=forbidden_edge_numbers)
	