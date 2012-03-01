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

import itertools

def degrees(g):

	n = g[0]
	edges = g[1]
	return [len([e for e in edges if x in e]) for x in range(1, n + 1)]

def generate_flags(n, tg, forbidden_edge_numbers={}):

	check_forbidden_edge_numbers = len(forbidden_edge_numbers) > 0

	s = tg[0]

	if n < s:
		return []

	if n == s:
		return [tg]

	max_ne = (n - 1) * (n - 2) / 2
	max_e = n * max_ne / 3
	
	new_graphs = [set() for i in range(max_e + 1)]
	
	smaller_graphs = generate_flags(n - 1, tg, forbidden_edge_numbers=forbidden_edge_numbers)
	
	possible_nbrs = [p for p in itertools.combinations(range(1, n), 2)]

	for sg in smaller_graphs:
	
		pe = len(sg[1])
		ds = degrees(sg)
		maxd = max(ds[s:] + [0])
			
		for ne in range(maxd, max_ne + 1):
		
			for nb in itertools.combinations(possible_nbrs, ne):
			
				ng = (n, sg[1] + tuple([(v[0], v[1], n) for v in nb]))

				if check_forbidden_edge_numbers:
					if has_forbidden_edge_numbers(ng, forbidden_edge_numbers, must_have_highest=True):
						continue

				mng = minimal_isomorph(ng, tg)
				new_graphs[pe + ne].add(mng)

	return [g for graphs in new_graphs for g in graphs]


def generate_graphs(n, forbidden_edge_numbers={}):
	
	tg = (0, ())
	return generate_flags(n, tg, forbidden_edge_numbers=forbidden_edge_numbers)


def flag_products (g, s, m, typs, flags):

	n = g[0]

	vertices = range(1, n + 1)

	num_typs = len(typs)
	num_flags = [len(fl) for fl in flags]
	
	pair_densities = dict(((i, j, k), 0) for i in range(num_typs)
		for j in range(num_flags[i]) for k in range(num_flags[i]))
	
	for tv in Permutations(vertices, s):
	
		tg = induced_subgraph(g, tv)
		
		if not tg in typs:
			continue
		
		tindex = typs.index(tg)
	
		non_typ_verts = [x for x in vertices if not x in tv]
	
		for fav in Combinations(non_typ_verts, m - s):

			fag = minimal_isomorph(induced_subgraph(g, tv + fav), tg)
			faindex = flags[tindex].index(fag)
			
			remaining_verts = [x for x in non_typ_verts if not x in fav]
			
			for fbv in Combinations(remaining_verts, m - s):

				fbg = minimal_isomorph(induced_subgraph(g, tv + fbv), tg)
				fbindex = flags[tindex].index(fbg)
				
				pair_densities[(tindex, faindex, fbindex)] += 1

	total = falling_factorial(n, s) * binomial(n - s, m - s) * binomial(n - m, m - s)

	for key in pair_densities.iterkeys():
		pair_densities[key] = Integer(pair_densities[key]) / total

	return pair_densities


