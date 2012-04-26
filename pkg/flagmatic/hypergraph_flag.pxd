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

from flag cimport Flag

# 35 + 42 + 7 = 84, 84 * 3 = 252
DEF MAX_NUMBER_OF_EDGE_INTS = 256
DEF MAX_NUMBER_OF_VERTICES = 35

cdef class HypergraphFlag (Flag):

	cdef int _n
	cdef int _r
	cdef bint _oriented
	cdef int _t
	cdef readonly bint is_degenerate
	cdef readonly int ne
	cdef int _edges[MAX_NUMBER_OF_EDGE_INTS]
	cpdef is_labelled_isomorphic(self, HypergraphFlag other)
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
