from flag cimport Flag

# 35 + 42 + 7 = 84, 84 * 3 = 252
DEF MAX_NUMBER_OF_EDGE_INTS = 256

cdef class HypergraphFlag (Flag):

	cdef int _n
	cdef int _r
	cdef bint _oriented
	cdef int _t
	cdef readonly bint is_degenerate
	cdef readonly int ne
	cdef int _edges[MAX_NUMBER_OF_EDGE_INTS]
	cpdef is_equal(self, HypergraphFlag other)
	cdef HypergraphFlag c_induced_subgraph(self, int *verts, int num_verts)
	cdef int c_has_subgraph (self, HypergraphFlag h)

	
cdef class combinatorial_info_block:
	cdef int np
	cdef int *p

cdef int *generate_permutations_fixing(int n, int s, int *number_of)
cdef int *generate_combinations(int n, int s, int *number_of)
cdef int *generate_combinations_plus(int n, int s, int *number_of)
cdef int *generate_pair_combinations(int n, int s, int m1, int m2, int *number_of)
cdef int *generate_equal_pair_combinations(int n, int s, int m, int *number_of)

cdef class graph_block:
	cdef int n, len
	cdef void **graphs
