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
	

# TODO: implement a forbidden_induced_edge_numbers.
# TODO: turn some forbidden graphs into forbidden edge numbers.

from sage.rings.arith import binomial
from sage.combinat.all import Combinations, Permutations
from sage.rings.all import Integer, QQ
from sage.matrix.all import matrix, block_matrix
from sage.modules.misc import gram_schmidt

from flag import *

def generate_flags(n, tg, r=3, oriented=False, forbidden_edge_numbers={}, forbidden_graphs = [], forbidden_induced_graphs=[]):
	"""
	For an integer n, and a type tg, returns a list of all tg-flags on n
	vertices, that satisfy certain constraints.
	
	forbidden_edge_numbers should be a dictionary whose keys and values are integers,
	where an item (k, v) specifies that k-sets are to span fewer than v edges.
		
	forbidden_graphs should be a list of graphs that are forbidden as subgraphs.
	
	forbidden_induced_subgraphs should be a list of graphs that are forbidden as
	_induced_ subgraphs.
	
	EXAMPLES:
	
		sage: tg = (1,())
		sage: generate_flags(4, tg, forbidden_edge_numbers={4:2})
		[(4, ()), (4, ((2, 3, 4),)), (4, ((1, 2, 3),))]
	
	"""

	if not (r == 2 or r == 3):
		raise NotImplementedError
		
	if oriented and r != 2:
		raise NotImplementedError

	if tg is None:
		tg = Flag(r=r, oriented=oriented)

	if r != tg.r or oriented != tg.oriented:
		raise ValueError

	if tg.t != 0:
		raise NotImplementedError("type must not contain labelled vertices.")

	s = tg.n

	if n < s:
		return []

	if n == s:
		ntg = tg.copy()
		ntg.t = s
		return [ntg]

	#max_ne = (n - 1) * (n - 2) / 2
	max_ne = binomial(n - 1, r - 1)

	#max_e = n * max_ne / 3
	max_e = binomial(n, r)
	
	new_graphs = []
	hashes = set()
	
	smaller_graphs = generate_flags(n - 1, tg, r, oriented, forbidden_edge_numbers=forbidden_edge_numbers,
		forbidden_graphs=forbidden_graphs, forbidden_induced_graphs=forbidden_induced_graphs)
	
	possible_edges = []

	if r == 3:
		for c in Combinations(range(1, n), 2):
			possible_edges.append((c[0], c[1], n))

	elif r == 2:
		for x in range(1, n):
			possible_edges.append((x, n))
			if oriented:
				possible_edges.append((n, x))

	for sg in smaller_graphs:
	
		pe = sg.ne
		ds = sg.degrees()
		maxd = max(ds[s:] + (0,))
			
		for ne in range(maxd, max_ne + 1):
		
			for nb in Combinations(possible_edges, ne):

				# For oriented graphs, can't have bidirected edges.
				# TODO: exclude these in a more efficient way!
				if oriented:
					if any(e in nb and (e[1], e[0]) in nb for e in possible_edges):
						continue
						
				ng = sg.copy()
				ng.n = n
				for e in nb:
					ng.add_edge(e)

				if ng.has_forbidden_edge_numbers(forbidden_edge_numbers, must_have_highest=True):
					continue

				if ng.has_forbidden_graphs(forbidden_graphs, must_have_highest=True):
					continue

				if ng.has_forbidden_graphs(forbidden_induced_graphs, must_have_highest=True, induced=True):
					continue

				ng.make_minimal_isomorph()
				ng_hash = hash(ng)
				if not ng_hash in hashes:
					new_graphs.append(ng)
					hashes.add(ng_hash)

	return new_graphs


def generate_graphs(n, r=3, oriented=False, forbidden_edge_numbers={}, forbidden_graphs = [], forbidden_induced_graphs=[]):
	"""
	For an integer n, return a list of all 3-graphs on n vertices that satisfy certain
	constraints.
	
	forbidden_edge_numbers should be a dictionary whose keys and values are integers,
	where an item (k, v) specifies that k-sets are to span fewer than v edges.
		
	forbidden_graphs should be a list of graphs that are forbidden as subgraphs.
	
	forbidden_induced_subgraphs should be a list of graphs that are forbidden as
	_induced_ subgraphs.
	
	EXAMPLES:
	
		sage: generate_graphs(4, forbidden_edge_numbers={4:3})
		[(4, ()), (4, ((1, 2, 3),)), (4, ((1, 2, 3), (1, 2, 4)))]
	
	"""
	return generate_flags(n, None, r, oriented, forbidden_edge_numbers=forbidden_edge_numbers,
		forbidden_graphs=forbidden_graphs, forbidden_induced_graphs=forbidden_induced_graphs)


def flag_orbits(tg, flags):
	"""
	flags should be a list of flags of the type tg. Returns a list of tuples.
	Each tuple contains the indices of the flags that are in the same orbit
	under the action of relabelling the vertices of tg.
	"""
	s = tg.n
	min_flags = []

	for fg in flags:
		mfgs = str(fg)
		for perm in Permutations(range(1, s + 1)):
			permplus = perm + range(s + 1, fg.n + 1)
			ntg = tg.copy()
			ntg.relabel(perm)
			nfg = fg.copy()
			nfg.relabel(permplus)
			nfg.make_minimal_isomorph()
			nfgs = str(nfg)
			if nfgs < mfgs:
				mfgs = nfgs
		min_flags.append(mfgs)

	orbs = []
	for mfgs in set(min_flags):
		orbs.append(tuple([i for i in range(len(min_flags)) if min_flags[i] == mfgs]))

	return sorted(orbs)


def flag_basis(tg, flags, orthogonalize=True):
	"""
	flags should be a list of flags of the type tg. Returns a basis for the flags
	that is a block matrix of two blocks. Uses flag orbits to create invariant-
	anti-invariant decomposition. If orthogonalize=True, perform Gram-Schmidt
	orthogonalization.
	"""
	
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

	#sys.stdout.write("Invariant-anti-invariant split: %d + %d = %d\n" % (Inv.nrows(), AntiInv.nrows(),
	#	len(flags)))
	
	if orthogonalize:
	
		# Note: the following does not preserve sparsity
		#AntiInv, mu = AntiInv.gram_schmidt()
	
		AntiInvRows, mu = gram_schmidt(AntiInv.rows())
		AntiInv = matrix(QQ, AntiInvRows, sparse=True)

	return block_matrix([[Inv],[AntiInv]])


def homomorphic_images(G):
	"""
	For an unlabelled Flag G, returns a list of Flags of order at most that of G,
	that are homomorphic images of G. The first Flag in the list will be a copy of G.
	
	"""

	if G.t != 0:
		raise ValueError

	mg = copy(G)
	mg.make_minimal_isomorph()

	if G.n <= 1:
		return [copy(G)]

	graph_hashes = set()
	graphs = []

	bad_pairs = set()
	
	for e in mg.edges:
		bad_pairs.add((e[0], e[1]))
		bad_pairs.add((e[0], e[2]))
		bad_pairs.add((e[1], e[2]))

	for i in range(1, G.n + 1):
		for j in range(i + 1, G.n + 1):
			
			if (i, j) in bad_pairs:
				continue
			
			ig = copy(mg)
			ig.identify_vertices(i, j)
			ig.make_minimal_isomorph()
	
			ghash = hash(ig)
			if not ghash in graph_hashes:
				graph_hashes.add(ghash)
				graphs.append(ig)
				s_graphs = homomorphic_images(ig)
				for sg in s_graphs:
					sghash = hash(sg)
					if not sghash in graph_hashes:
						graph_hashes.add(sghash)
						graphs.append(sg)

	if len(graphs) == 0:
		return []

	return [copy(G)] + graphs
	
