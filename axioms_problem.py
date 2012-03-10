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

	
	def __init__(self):
	
		Problem.__init__(self)
		self._quantum_graphs = []


	def clear_axioms(self):
		
		self._quantum_graphs = []


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
		
		self._quantum_graphs.extend(quantum_graphs)		
	
	
	def add_codegree_axiom(self, value):
	
		tg = Flag("2:")
		f1 = Flag("3:123", tg)
		f2 = Flag("2:", tg)
		self.add_axiom(tg, [f1, f2], [Integer(1), -value])


	def add_degree_axiom(self, value):
	
		tg = Flag("1:")
		f1 = Flag("3:123", tg)
		f2 = Flag("1:", tg)
		self.add_axiom(tg, [f1, f2], [Integer(1), -value])


	def process_sdp_output_file(self, f):

		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._quantum_graphs)
		
		inv_block_sizes = []
		anti_inv_block_sizes = []
		for i in range(num_types):
			row_div = self._flag_bases[i].subdivisions()[0]
			if len(row_div) > 0:
				inv_block_sizes.append(row_div[0])
				anti_inv_block_sizes.append(self._flag_bases[i].nrows() - row_div[0])
			else:
				inv_block_sizes.append(self._flag_bases[i].nrows())
				anti_inv_block_sizes.append(1)	
		
		self._quantum_graphs_fcoeffs = [0.0 for i in range(num_densities)]
		
		for line in f:
			numbers = line.split()
			if numbers[0] != "2":
				continue
			ti = int(numbers[1]) - 2
			if ti == 2 * num_types + 1:
				j = int(numbers[2]) - 1
				self._quantum_graphs_fcoeffs[j] = RDF(numbers[4])
				continue
			if ti < 0 or ti >= 2 * num_types:
				continue
			j = int(numbers[2]) - 1
			k = int(numbers[3]) - 1
			if ti % 2:
				ti = (ti - 1) / 2
				j += inv_block_sizes[ti]
				k += inv_block_sizes[ti]
				if j >= self._flag_bases[ti].nrows():
					continue # Might be 'dummy' block for types with no anti-inv space
			else:
				ti /= 2
			self._sdp_Q_matrices[ti][j, k] = numbers[4]
			self._sdp_Q_matrices[ti][k, j] = self._sdp_Q_matrices[ti][j, k]


	
	def write_sdp_input_file(self):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._quantum_graphs)
		
		self._obj_value_factor = 1
		
		if num_densities < 1:
			raise NotImplementedError("at least one axiom must be provided.")
		
		self._sdp_input_filename = os.path.join(SAGE_TMP, "sdp.dat-s")
	
		with open(self._sdp_input_filename, "w") as f:
	
			f.write("%d\n" % (num_graphs + 1,))
			f.write("%d\n" % (2 * num_types + 3,))
			
			f.write("1 ")
			self.write_block_sizes(f)
			f.write("-%d -%d\n" % (num_graphs, num_densities))
			f.write("%s1.0\n" % ("0.0 " * num_graphs,))
			f.write("0 1 1 1 -1.0\n")
	
			for i in range(num_graphs):
				f.write("%d 1 1 1 -1.0\n" % (i + 1,))
				# TODO: omit for sharp graphs
				f.write("%d %d %d %d 1.0\n" % (i + 1, 2 * num_types + 2, i + 1, i + 1))
	
			for i in range(num_graphs):
				for j in range(num_densities):
					d = self._quantum_graphs[j][i]
					if d != 0:
						f.write("%d %d %d %d %s\n" % (i + 1, 2 * num_types + 3, j + 1, j + 1, d.n(digits=64)))
			
			for j in range(num_densities):
				f.write("%d %d %d %d 1.0\n" % (num_graphs + 1, 2 * num_types + 3, j + 1, j + 1))
		
			self.write_blocks(f)


	def write_alternate_sdp_input_file(self):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._quantum_graphs)
		
		self._obj_value_factor = 1
		
		if num_densities < 1:
			raise NotImplementedError("at least one axiom must be provided.")
		
		self._sdp_input_filename = os.path.join(SAGE_TMP, "sdp.dat-s")
	
		with open(self._sdp_input_filename, "w") as f:
	
			f.write("%d\n" % (num_graphs + num_densities * 2,))
			f.write("%d\n" % (2 * num_types + 5,))
			
			f.write("1 ")
			self.write_block_sizes(f)
			f.write("-%d -%d -%d -%d\n" % (num_graphs, num_densities, num_densities, num_densities))
			f.write("%s%s%s\n" % ("0.0 " * num_graphs, "0.0 " * num_densities, "1.0 " * num_densities))
			f.write("0 1 1 1 1.0\n")
	
			for i in range(num_graphs):
				# TODO: omit for sharp graphs
				f.write("%d %d %d %d 1.0\n" % (i + 1, 2 * num_types + 2, i + 1, i + 1))
	
			for i in range(num_graphs):
				for j in range(num_densities):
					d = self._quantum_graphs[j][i]
					if d != 0:
						f.write("%d %d %d %d %s\n" % (i + 1, 2 * num_types + 3, j + 1, j + 1, d.n(digits=64)))
			
			for j in range(num_densities):
				f.write("%d 1 1 1 -1.0\n" % (num_graphs + 1 + j,))
				f.write("%d %d %d %d 1.0\n" % (num_graphs + 1 + j, 2 * num_types + 3, j + 1, j + 1))
				f.write("%d %d %d %d -1.0\n" % (num_graphs + 1 + j, 2 * num_types + 4, j + 1, j + 1))
		
			for j in range(num_densities):
				f.write("%d %d %d %d 1.0\n" % (num_graphs + num_densities + 1 + j, 2 * num_types + 3, j + 1, j + 1))
				f.write("%d %d %d %d 1.0\n" % (num_graphs + num_densities + 1 + j, 2 * num_types + 5, j + 1, j + 1))
	
			self.write_blocks(f)
	
			
	def find_sharps(self, tolerance = 0.00001):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._quantum_graphs)
		
		fbounds = [0.0 for i in range(num_graphs)]

		for i in range(num_graphs):
			for j in range(num_densities):
				fbounds[i] += self._quantum_graphs_fcoeffs[j] * self._quantum_graphs[j][i]
				
		for ti in range(num_types):
			num_flags = len(self._flags[ti])
			for gi in range(num_graphs):
				D = sparse_symm_matrix_from_compact_repr(self._product_densities[(gi, ti)])
				for j in range(D.nrows()):
					for k in range(j, D.nrows()):
						value = self._sdp_Q_matrices[ti][j, k]
						if j != k:
							value *= 2
						d = D[j, k]
						if d != 0:
							fbounds[gi] += d * value

		bound = max(fbounds)
		sharp_indices = [gi for gi in range(num_graphs)
			if abs(fbounds[gi] - bound) < tolerance]
		
		sys.stdout.write("The following %d graphs are sharp:\n" % len(sharp_indices))
		for gi in sharp_indices:
			sys.stdout.write("%s : graph %d (%s)\n" % (fbounds[gi], gi + 1, self._graphs[gi]))
