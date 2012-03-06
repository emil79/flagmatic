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

import os
import pexpect
import sys

# pexpect in Sage has a bug, which prevents it using commands with full paths.
# So for now, CSDP has to be in a directory in $PATH.

cdsp_cmd = "csdp"

class flagmatic_problem(object):

	forbidden_edge_numbers = {}
	forbidden_graphs = []
	forbidden_induced_graphs = []

	def __init__(self):
		self._n = 0
		self._graphs = []
		self._types = []
		self._flags = []
		
	@property
	def n(self):
		return self._n

	@n.setter
	def n(self, n):
		self._n = n

		sys.stdout.write("Generating graphs...\n")
		self._graphs = generate_graphs(n, forbidden_edge_numbers=self.forbidden_edge_numbers)
		self._graph_block = make_graph_block(self._graphs, n)
		sys.stdout.write("Generated %d graphs.\n" % len(self._graphs))
	
		self._graph_densities = []
		for g in self._graphs:
			self._graph_densities.append(len(g[1]) / binomial(n, 3))
	
		sys.stdout.write("Generating types and flags...\n")
		self._types = []
		self._flags = []
	
		for s in range(n % 2, n - 1, 2):
			
			these_types = generate_graphs(s, forbidden_edge_numbers=self.forbidden_edge_numbers)
			sys.stdout.write("Generated %d types of order %d, " % (
				len(these_types), s))

			m = (n + s) / 2
			these_flags = []
			for tg in these_types:
				these_flags.append(generate_flags(m, tg, forbidden_edge_numbers=self.forbidden_edge_numbers))
			sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))
						
			self._types.extend(these_types)
			self._flags.extend(these_flags)

		self._flag_bases = []
		for ti in range(len(self._types)):
			self._flag_bases.append(identity_matrix(QQ, len(self._flags[ti]), sparse=True))

	@property
	def graphs(self):
		return self._graphs

	@property
	def types(self):
		return self._types

	@property
	def flags(self):
		return self._flags
	
	def set_inv_anti_inv_bases(self):

		for ti in range(len(self._types)):
			B = flag_basis(self._types[ti], self._flags[ti])
			self._flag_bases[ti] = B

 	def calculate_flag_products(self):
 	
		self._flag_products = []
 		sys.stdout.write("Averaging flag products...\n")

 		for ti in range(len(self._types)):

			B = self._flag_bases[ti]
			row_div = B.subdivisions()[0]
			try:
				inv_rows = row_div[0]
				is_subdivided = True
			except IndexError:
				is_subdivided = False

 			tg = self._types[ti]
 			s = tg[0]
 			m = (self._n + s) / 2

 			sys.stdout.write("Doing type %d (order %d; flags %d).\n" % (ti + 1, s, m))
 			flags_block = make_graph_block(self._flags[ti], m)
			DL = flag_products(self._graph_block, make_graph_block([tg], s), flags_block, None)
		
			nfp = []
			for gi in range(len(self._graphs)):
				D = DL[gi]
				ND = B * D * B.T
				if is_subdivided:
					ND.subdivide(inv_rows, inv_rows)
				nfp.append(ND)
			self._flag_products.append(nfp)
		
	def write_sdp_input_file(self):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		
		inv_block_sizes = []
		anti_inv_block_sizes = []
		for i in range(num_types):
			row_div = self._flag_bases[i].subdivisions()[0]
			if len(row_div) > 0:
				inv_block_sizes.append(row_div[0])
				anti_inv_block_sizes.append(len(self._flags[i]) - row_div[0])
			else:
				inv_block_sizes.append(len(self._flags[i]))
				anti_inv_block_sizes.append(1)

		self._sdp_input_filename = os.path.join(SAGE_TMP, "sdp.dat-s")
	
		with open(self._sdp_input_filename, "w") as f:
	
			f.write("%d\n" % (num_graphs + 1,))
			f.write("%d\n" % (2 * num_types + 3,))
			f.write("1 ")
			for i in range(num_types):
				f.write("%d " % inv_block_sizes[i])
				f.write("%d " % anti_inv_block_sizes[i])
			f.write("-%d -1\n" % num_graphs)
			f.write("%s1.0\n" % ("0.0 " * num_graphs,))
			f.write("0 1 1 1 -1.0\n")
	
			for i in range(num_graphs):
				f.write("%d 1 1 1 -1.0\n" % (i + 1,))
				f.write("%d %d %d %d 1.0\n" % (i + 1, 2 * num_types + 2, i + 1, i + 1))
	
			for i in range(num_graphs):
				d = self._graph_densities[i]
				if d != 0:
					f.write("%d %d 1 1 %s\n" % (i + 1, 2 * num_types + 3, d.n(digits=64)))
			
			f.write("%d %d 1 1 1.0\n" % (num_graphs + 1, 2 * num_types + 3))
		
			for i in range(num_graphs):
				for j in range(num_types):
					for key, value in self._flag_products[j][i].dict().iteritems():
						row, col = key
						if row < col: # only print upper triangle
							continue
						if row < inv_block_sizes[j]:
							block = 2 * j + 2
						else:
							block = 2 * j + 3
							row -= inv_block_sizes[j]
							col -= inv_block_sizes[j]
						f.write("%d %d %d %d %s\n" % (i + 1, block, row + 1, col + 1,
							value.n(digits=64)))

	def run_csdp(self, show_output=False):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		
		inv_block_sizes = []
		anti_inv_block_sizes = []
		for i in range(num_types):
			row_div = self._flag_bases[i].subdivisions()[0]
			if len(row_div) > 0:
				inv_block_sizes.append(row_div[0])
				anti_inv_block_sizes.append(len(self._flags[i]) - row_div[0])
			else:
				inv_block_sizes.append(len(self._flags[i]))
				anti_inv_block_sizes.append(1)
	
		self._sdp_output_filename = os.path.join(SAGE_TMP, "sdp.out")
	
		cmd = "%s %s %s" % (cdsp_cmd, self._sdp_input_filename, self._sdp_output_filename)
		p = pexpect.spawn(cmd, timeout=60*60*24*7)
		obj_val = None
		while True:
			if p.eof():
				break
			try:
				p.expect("\r\n")
				line = p.before.strip() + "\n"
				if "Primal objective value:" in line:
					obj_val = -RDF(line.split()[-1])
				if show_output:
					sys.stdout.write(line)
			except pexpect.EOF:
				break
		p.close()
		returncode = p.exitstatus
		
		sys.stdout.write("Returncode is %d. Objective value is %s.\n" % (returncode,
			obj_val))
		
		self._sdp_Q_matrices = [matrix(RDF, len(self._flags[i]), len(self._flags[i]))
			for i in range(num_types)]
		
		with open(self._sdp_output_filename, "r") as f:
			
			for line in f:
				numbers = line.split()
				if numbers[0] != "2":
					continue
				ti = int(numbers[1]) - 2
				if ti < 0 or ti >= 2 * num_types:
					continue
				j = int(numbers[2]) - 1
				k = int(numbers[3]) - 1
				if ti % 2:
					ti = (ti - 1) / 2
					j += inv_block_sizes[ti]
					k += inv_block_sizes[ti]
					if j >= len(self._flags[ti]):
						continue # Might be 'dummy' block for types with no anti-inv space
				else:
					ti /= 2
				self._sdp_Q_matrices[ti][j, k] = numbers[4]
				self._sdp_Q_matrices[ti][k, j] = self._sdp_Q_matrices[ti][j, k]

	def find_sharps(self):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		fbounds = [RDF(self._graph_densities[i]) for i in range(num_graphs)]
		
		for ti in range(num_types):
			num_flags = len(self._flags[ti])
			for j in range(num_flags):
				for k in range(j, num_flags):
					value = self._sdp_Q_matrices[ti][j, k]
					if j != k:
						value *= 2
					for gi in range(num_graphs):
						key = (j, k)
						d = self._flag_products[ti][gi].dict()
						if key in d:
							fbounds[gi] += d[key] * value

		tolerance = 0.00001
		bound = max(fbounds)
		sharp_indices = [gi for gi in range(num_graphs)
			if abs(fbounds[gi] - bound) < tolerance]
		
		sys.stdout.write("The following %d graphs are sharp:\n" % len(sharp_indices))
		for gi in sharp_indices:
			sys.stdout.write("%s : graph %d (%s)\n" % (fbounds[gi], gi + 1, self._graphs[gi]))
		

def test_sdp(n, show_output=False):
	P = flagmatic_problem()
	P.forbidden_edge_numbers={4:3}
	P.n = n
	P.set_inv_anti_inv_bases()
	P.calculate_flag_products()
	P.write_sdp_input_file()
	P.run_csdp(show_output=show_output)
	return P
