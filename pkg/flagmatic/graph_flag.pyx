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

from hypergraph_flag cimport HypergraphFlag

cdef class GraphFlag (HypergraphFlag):

	def __init__(self, string_rep=None, r=3, oriented=False):
		super(GraphFlag, self).__init__(string_rep=string_rep, r=2, oriented=False)

	@classmethod
	def generate_flags(cls, n, tg, forbidden_edge_numbers={}, forbidden_graphs=[], forbidden_induced_graphs=[]):
		return HypergraphFlag.generate_flags(n, tg, r=2, oriented=False, forbidden_edge_numbers={},
			forbidden_graphs=[], forbidden_induced_graphs=[])

	@classmethod
	def generate_graphs(cls, n, forbidden_edge_numbers={}, forbidden_graphs=[], forbidden_induced_graphs=[]):
		return HypergraphFlag.generate_flags(n, cls(), r=2, oriented=False, forbidden_edge_numbers={},
			forbidden_graphs=[], forbidden_induced_graphs=[])
