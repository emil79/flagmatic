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

from sage.arith.misc import binomial
from hypergraph_flag cimport HypergraphFlag

cdef class ThreeGraphFlag (HypergraphFlag):


	def __init__(self, representation=None):
		super(ThreeGraphFlag, self).__init__(representation=representation, r=3, oriented=False)


	def __reduce__(self):
		return (type(self), (self._repr_(),))


	@classmethod
	def description(cls):
		return "3-graph"


	@classmethod
	def default_density_graph(cls):
		return cls("3:123")


	@classmethod
	def max_number_edges(cls, n):
		return binomial(n, 3)


	@classmethod
	def generate_flags(cls, n, tg, forbidden_edge_numbers=None, forbidden_graphs=None, forbidden_induced_graphs=None):
		return HypergraphFlag.generate_flags(n, tg, r=3, oriented=False, forbidden_edge_numbers=forbidden_edge_numbers,
			forbidden_graphs=forbidden_graphs, forbidden_induced_graphs=forbidden_induced_graphs)


	@classmethod
	def generate_graphs(cls, n, forbidden_edge_numbers=None, forbidden_graphs=None, forbidden_induced_graphs=None):
		return HypergraphFlag.generate_flags(n, cls(), r=3, oriented=False, forbidden_edge_numbers=forbidden_edge_numbers,
			forbidden_graphs=forbidden_graphs, forbidden_induced_graphs=forbidden_induced_graphs)
