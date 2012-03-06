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

def degrees(g):

	n = g[0]
	edges = g[1]
	return [len([e for e in edges if x in e]) for x in range(1, n + 1)]


def generate_flags(n, tg, forbidden_edge_numbers={}, forbidden_graphs = [], forbidden_induced_graphs=[]):

	# TODO: implement forbidden_induced_graphs - currently ignored.
	# TODO: turn some forbidden graphs into forbidden edge numbers.

	check_forbidden_edge_numbers = len(forbidden_edge_numbers) > 0
	check_forbidden_graphs = len(forbidden_graphs) > 0

	if check_forbidden_graphs:
		fg_block = make_graph_block(forbidden_graphs, 0)

	s = tg[0]

	if n < s:
		return []

	if n == s:
		return [tg]

	max_ne = (n - 1) * (n - 2) / 2
	max_e = n * max_ne / 3
	
	new_graphs = [set() for i in range(max_e + 1)]
	
	smaller_graphs = generate_flags(n - 1, tg, forbidden_edge_numbers=forbidden_edge_numbers,
		forbidden_graphs=forbidden_graphs, forbidden_induced_graphs=forbidden_induced_graphs)
	
	possible_nbrs = [p for p in Combinations(range(1, n), 2)]

	for sg in smaller_graphs:
	
		pe = len(sg[1])
		ds = degrees(sg)
		maxd = max(ds[s:] + [0])
			
		for ne in range(maxd, max_ne + 1):
		
			for nb in Combinations(possible_nbrs, ne):
			
				ng = (n, sg[1] + tuple([(v[0], v[1], n) for v in nb]))

				if check_forbidden_edge_numbers:
					if has_forbidden_edge_numbers(ng, forbidden_edge_numbers, must_have_highest=True):
						continue

				if check_forbidden_graphs:
					if has_forbidden_graphs(ng, fg_block, must_have_highest=True):
						continue

				mng = minimal_isomorph(ng, tg)
				new_graphs[pe + ne].add(mng)

	return [g for graphs in new_graphs for g in graphs]


def generate_graphs(n, forbidden_edge_numbers={}, forbidden_graphs = [], forbidden_induced_graphs=[]):
	
	tg = (0, ())
	return generate_flags(n, tg, forbidden_edge_numbers=forbidden_edge_numbers,
		forbidden_graphs=forbidden_graphs, forbidden_induced_graphs=forbidden_induced_graphs)


def flag_orbits(tg, flags):

	s = tg[0]
	min_flags = []

	for fg in flags:
		mfg = fg
		for perm in Permutations(range(1, s + 1)):
			permplus = perm + range(s + 1, fg[0] + 1)
			ntg = (tg[0], tuple(sorted([tuple(sorted([perm[e[i] - 1] for i in range(3)])) for e in tg[1]])))
			nfg = (fg[0], tuple(sorted([tuple(sorted([permplus[e[i] - 1] for i in range(3)])) for e in fg[1]])))
			mnfg = minimal_isomorph(nfg, ntg)
			if mnfg < mfg:
				mfg = mnfg
		min_flags.append(mfg)

	orbs = []
	for mfg in set(min_flags):
		orbs.append(tuple([i for i in range(len(min_flags)) if min_flags[i] == mfg]))

	return sorted(orbs)
	

from sage.modules.misc import gram_schmidt

def flag_basis(tg, flags, orthogonalize=True):

	orbs = flag_orbits(tg, flags)
	
	Inv = matrix(QQ, len(orbs), len(flags), sparse=True)
	row = 0
	for orb in orbs:
		for j in orb:
			Inv[row, j] = 1
		row += 1
	
	# There might be no anti-invariant part.
	if len(orbs) == len(flags):
		return Inv
	
	AntiInv = matrix(QQ, len(flags) - len(orbs), len(flags), sparse=True)
	row = 0
	for orb in orbs:
		for j in orb[1:]:
			AntiInv[row, orb[0]] = 1
			AntiInv[row, j] = -1
			row += 1

	sys.stdout.write("Inv-AntiInv: %d + %d = %d\n" % (Inv.nrows(), AntiInv.nrows(),
		len(flags)))
	
	if orthogonalize:
	
		# Note: the following does not preserve sparsity
		#AntiInv, mu = AntiInv.gram_schmidt()
	
		AntiInvRows, mu = gram_schmidt(AntiInv.rows())
		AntiInv = matrix(QQ, AntiInvRows, sparse=True)

	return block_matrix([[Inv],[AntiInv]])


# Deprecated: use flag_products instead

def slow_flag_products (g, s, m, typs, flags):

	n = g[0]

	vertices = range(1, n + 1)

	num_typs = len(typs)
	num_flags = [len(fl) for fl in flags]
	
	#pair_densities = dict(((i, j, k), 0) for i in range(num_typs)
	#	for j in range(num_flags[i]) for k in range(num_flags[i]))
	pair_densities = {}
	
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
				
				key = (tindex, faindex, fbindex)
				pair_densities[key] = pair_densities.setdefault(key, 0) + 1

	total = falling_factorial(n, s) * binomial(n - s, m - s) * binomial(n - m, m - s)

	for key in pair_densities.iterkeys():
		pair_densities[key] = Integer(pair_densities[key]) / total

	return pair_densities


def induced_subgraph (g, S):

	good_edges = [e for e in g[1] if all(x in S for x in e)]
	p = [0 for i in range(g[0] + 1)]
	for i in range(len(S)):
		p[S[i]] = i + 1

	edges = sorted([tuple(sorted([p[x] for x in e])) for e in good_edges])

	return (len(S), tuple(edges))


# Deprecated: use minimal_isomorph instead

def slow_minimal_isomorph (g):
	
	n = g[0]
	min_edges = g[1]
	
	for p in Permutations(range(1, n + 1)):
		
		edges = tuple(sorted([tuple(sorted([p[e[i] - 1] for i in range(2)]))
			for e in g[1]]))
		
		if edges < min_edges:
			min_edges = edges
			
	return (n, min_edges)


def sparse_symm_matrix_to_compact_repr(M):

	ed = {}
	for key, value in M.dict().iteritems():
		x, y = key
		if x <= y:
			ed[key] = repr(value)

	d = {
		"n" : M.nrows(),
		"blocks" : M.subdivisions()[0],
		"entries" : ed
	}

	return repr(d)
	

def sparse_symm_matrix_from_compact_repr(ds):

	d = eval(ds)
	n = d["n"]
	M = matrix(QQ, n, n, sparse=True)
	for key, value in d["entries"].iteritems():
		x, y = key
		M[x, y] = sage_eval(value)
		M[y, x] = M[x, y]
	M.subdivide(d["blocks"], d["blocks"])
	return M
	