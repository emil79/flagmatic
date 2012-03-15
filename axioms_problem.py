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


class AxiomsProblem(Problem):
	
	def __init__(self, r=3, oriented=False):
	
		Problem.__init__(self, r, oriented)
	

	def save_more_json(self, d):
		pass

	
	@classmethod
	def load_more_json(cls, d, obj):
		pass	
	

	def clear_axioms(self):
		
		self._densities = []


	def add_axiom(self, tg, flags, flag_coeffs):

		m = self.n - max([f.n for f in flags]) + tg.n

		axiom_flags = generate_flags(m, tg,
			forbidden_edge_numbers=self.forbidden_edge_numbers,
			forbidden_graphs=self.forbidden_graphs,
			forbidden_induced_graphs=self.forbidden_induced_graphs)
		
		num_densities = len(axiom_flags)
		sys.stdout.write("Added %d quantum graphs.\n" % num_densities)
		
		num_graphs = len(self._graphs)
		quantum_graphs = [[Integer(0) for i in range(num_graphs)] for j in range(num_densities)]
		
		axiom_flags_block = make_graph_block(axiom_flags, m)
		graph_block = make_graph_block(self._graphs, self.n)

		for i in range(len(flags)):
			flags_block = make_graph_block([flags[i]], flags[i].n)
			DL = flag_products(graph_block, tg, flags_block, axiom_flags_block)
			for j in range(num_graphs):
				M = DL[j]
				for k in range(num_densities):
					quantum_graphs[k][j] += M[0, k] * flag_coeffs[i]
		
		self._densities.extend(quantum_graphs)		
	
	
	def add_codegree_axiom(self, value):
	
		tg = Flag("2:")
		f1 = Flag("3:123(2)")
		f2 = Flag("2:(2)")
		self.add_axiom(tg, [f1, f2], [Integer(1), -value])


	def add_degree_axiom(self, value):
	
		tg = Flag("1:")
		f1 = Flag("3:123(1)")
		f2 = Flag("1:(1)")
		self.add_axiom(tg, [f1, f2], [Integer(1), -value])
