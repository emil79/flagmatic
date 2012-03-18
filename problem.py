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

import base64
import gzip
import json
import numpy
import os
import pexpect
import sys

from sage.all import Integer, QQ, matrix
from sage.structure.sage_object import SageObject


# pexpect in Sage has a bug, which prevents it using commands with full paths.
# So for now, CSDP has to be in a directory in $PATH.

cdsp_cmd = "csdp"

    
class Problem(SageObject):

	def __init__(self, r=3, oriented=False):
	
		self._n = 0
		self._r = r
		self._oriented = oriented

		self.forbidden_edge_numbers = {}
		self.forbidden_graphs = []
		self.forbidden_induced_graphs = []
		self._graphs = []
		self._densities = []
		self._density_graph = None

		self._types = []
		self._flags = []
		self._block_bases = []
		self._flag_bases = []
		self._target_bound = None
		self._sharp_graphs = []
		self._zero_eigenvectors = []
		self._product_densities_dumps = {}
	
		self._obj_value_factor = -1
		self._minimize = False
		self._force_sharps = False
		self._sdp_input_filename = None
		self._sdp_output_filename = None
	
		self._sdp_Q_matrices = []
		self._sdp_density_coeffs = []
		self._exact_Q_matrices = []
		self._exact_density_coeffs = []
		self._bounds = []


	def save_json(self, filename):
		
		d = {}
		
		d["n"] = self._n
		d["r"] = self._r
		d["oriented"] = self._oriented
		d["forbidden_edge_numbers"] = self.forbidden_edge_numbers
		d["forbidden_graphs"] = [repr(g) for g in self.forbidden_graphs]
		d["forbidden_induced_graphs"] = [repr(g) for g in self.forbidden_induced_graphs]
		
		d["graphs"] = [repr(g) for g in self._graphs]
		d["densities"] = [[repr(r) for r in x] for x in self._densities]
		if not self._density_graph is None:
			d["density_graph"] = repr(self._density_graph)
			
		d["types"] = [repr(g) for g in self._types]
		d["flags"] = [[repr(g) for g in x] for x in self._flags]
		d["block_bases"] = [base64.b64encode(dumps(M)) for M in self._block_bases]
		d["flag_bases"] = [base64.b64encode(dumps(M)) for M in self._flag_bases]
		if not self._target_bound is None:
			d["target_bound"] = repr(self._target_bound)
		d["sharp_graphs"] = self._sharp_graphs
		d["zero_eigenvectors"] = [base64.b64encode(dumps(M)) for M in self._zero_eigenvectors]
		d["product_densities_dumps"] = [[base64.b64encode(s) for s in x] for x in self._product_densities_dumps]

		d["obj_value_factor"] = repr(self._obj_value_factor)
		d["minimize"] = self._minimize
		d["force_sharps"] = self._force_sharps
		if not self._sdp_input_filename is None:
			d["sdp_input_filename"] = self._sdp_input_filename,
		if not self._sdp_output_filename is None:
			d["sdp_output_filename"] = self._sdp_output_filename
			
		d["sdp_Q_matrices"] = [base64.b64encode(dumps(M)) for M in self._sdp_Q_matrices]
		d["sdp_density_coeffs"] = [repr(n) for n in self._sdp_density_coeffs]
		d["exact_Q_matrices"] = [base64.b64encode(dumps(M)) for M in self._exact_Q_matrices]
		d["exact_density_coeffs"] = [repr(r) for r in self._exact_density_coeffs]
		d["bounds"] = [repr(r) for r in self._bounds]
		
		self.save_more_json(d)
		
		with open(filename, "wb") as f:
			json.dump(d, f)


	def save_more_json(self, d):
		pass		


	@classmethod
	def load_json(cls, filename):
		
		obj = cls()
		
		with open(filename, "rb") as f:
			d = json.load(f)
		
		if "n" in d:
			obj._n = d["n"]
		if "r" in d:
			obj._r = d["r"]
		if "oriented" in d:
			obj._oriented = d["oriented"]
		if "forbidden_edge_numbers" in d:	
			obj.forbidden_edge_numbers = d["forbidden_edge_numbers"]
		if "forbidden_graphs" in d:
			obj.forbidden_graphs = [Flag(s) for s in d["forbidden_graphs"]]
		if "forbidden_induced_graphs" in d:
			obj.forbidden_induced_graphs = [Flag(s) for s in d["forbidden_induced_graphs"]]
		
		if "graphs" in d:
			obj._graphs = [Flag(s) for s in d["graphs"]]
		if "densities" in d:
			obj._densities = [[sage_eval(s) for s in x] for x in d["densities"]]
		if "density_graph" in d:
			obj._density_graph = Flag(d["density_graph"])

		if "types" in d:
			obj._types = [Flag(s) for s in d["types"]]
		if "flags" in d:
			obj._flags = [[Flag(s) for s in x] for x in d["flags"]]
		if "block_bases" in d:
			obj._block_bases = [loads(base64.b64decode(s)) for s in d["block_bases"]]
		if "flag_bases" in d:
			obj._flag_bases = [loads(base64.b64decode(s)) for s in d["flag_bases"]]
		if "target_bound" in d:
			obj._target_bound = sage_eval(d["target_bound"]) 
		
		if "sharp_graphs" in d:
			obj._sharp_graphs = d["sharp_graphs"]
		if "zero_eigenvectors" in d:
			obj._zero_eigenvectors = [loads(base64.b64decode(s)) for s in d["zero_eigenvectors"]]
		if "product_densities_dumps" in d:
			obj._product_densities_dumps = [[base64.b64decode(s) for s in x] for x in d["product_densities_dumps"]]

		if "obj_value_factor" in d:
			obj._obj_value_factor = RDF(d["obj_value_factor"])
		if "minimize" in d:
			obj._minimize = d["minimize"]
		if "force_sharps" in d:
			obj._force_sharps = d["force_sharps"]
		if "sdp_input_filename" in d:
			obj._sdp_input_filename = d["sdp_input_filename"]
		if "sdp_output_filename" in d:
			obj._sdp_output_filename = d["sdp_output_filename"]
			
		if "sdp_Q_matrices" in d:
			obj._sdp_Q_matrices = [loads(base64.b64decode(s)) for s in d["sdp_Q_matrices"]]
		if "sdp_density_coeffs" in d:
			obj._sdp_density_coeffs = [RDF(s) for s in d["sdp_density_coeffs"]]
		if "exact_Q_matrices" in d:
			obj._exact_Q_matrices = [loads(base64.b64decode(s)) for s in d["exact_Q_matrices"]]
		if "exact_density_coeffs" in d:
			obj._exact_density_coeffs = [sage_eval(s) for s in d["exact_density_coeffs"]]
		if "bounds" in d:
			obj._bounds = [sage_eval(s) for s in d["bounds"]]
	
		cls.load_more_json(d, obj)
	
		return obj

	@classmethod
	def load_more_json(cls, d, obj):
		pass

	@property
	def r(self):
		return self._r
		
	@property
	def oriented(self):
		return self._oriented

	@property
	def n(self):
		return self._n

	@n.setter
	def n(self, n):
		self._n = n

		sys.stdout.write("Generating graphs...\n")
		self._graphs = generate_graphs(n, self._r, self._oriented,
			forbidden_edge_numbers=self.forbidden_edge_numbers,
			forbidden_graphs=self.forbidden_graphs,
			forbidden_induced_graphs=self.forbidden_induced_graphs)
		sys.stdout.write("Generated %d graphs.\n" % len(self._graphs))
	
		self._densities = [[g.edge_density() for g in self._graphs]]
	
		sys.stdout.write("Generating types and flags...\n")
		self._types = []
		self._flags = []
	
		for s in range(n % 2, n - 1, 2):
			
			these_types = generate_graphs(s, self._r, self._oriented,
				forbidden_edge_numbers=self.forbidden_edge_numbers,
				forbidden_graphs=self.forbidden_graphs,
				forbidden_induced_graphs=self.forbidden_induced_graphs)
			sys.stdout.write("Generated %d types of order %d, " % (
				len(these_types), s))

			m = (n + s) / 2
			these_flags = []
			for tg in these_types:
				these_flags.append(generate_flags(m, tg, self._r, self._oriented,
					forbidden_edge_numbers=self.forbidden_edge_numbers,
					forbidden_graphs=self.forbidden_graphs,
					forbidden_induced_graphs=self.forbidden_induced_graphs))
			sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))
						
			self._types.extend(these_types)
			self._flags.extend(these_flags)

		self._flag_bases = [identity_matrix(QQ, len(self._flags[ti]), sparse=True) for ti in range(len(self._types))]

		for ti in range(len(self._types)):
			self._flag_bases[ti].set_immutable()


	@property
	def graphs(self):
		return self._graphs

	@property
	def types(self):
		return self._types

	@property
	def flags(self):
		return self._flags

	# TODO: default self._density_graph to Flag("3:123")
	
	@property
	def density_graph(self):
		return self._density_graph

	@density_graph.setter
	def density_graph(self, dg):
		self._density_graph = dg
		self._densities = [[g.subgraph_density(dg) for g in self._graphs]]


	def remove_types(self, indices):
	
		num_types = len(self._types)
		remaining = [i for i in range(num_types) if not i in indices]

		sys.stdout.write("%d types removed.\n" % (num_types - len(remaining)))

		if num_types > 0:
			self._types = [self._types[i] for i in remaining]
		if len(self._flags) > 0:
			self._flags = [self._flags[i] for i in remaining]
			sys.stdout.write("Remaining types have %s flags.\n" % [len(flags) for flags in self._flags])
		if len(self._flag_bases) > 0:
			self._flag_bases = [self._flag_bases[i] for i in remaining]
		if len(self._block_bases) > 0:
			self._block_bases = [self._block_bases[i] for i in remaining]


	def create_block_bases(self):

		self._block_bases = []
		for ti in range(len(self._types)):
			B = flag_basis(self._types[ti], self._flags[ti])
			self._block_bases.append(B)


	def use_construction(self, c):

		# TODO: axioms construction needs density of zero
		# maybe create a self.density(g) function

		sys.stdout.write("Determining which graphs appear in construction...\n")
		
		sharp_graphs = c.subgraph_densities(self._n)
		self._sharp_graphs = []
		for pair in sharp_graphs:
			g, den = pair
			for gi in range(len(self._graphs)):
				if g.is_equal(self._graphs[gi]):
					self._sharp_graphs.append(gi)
					break
			else:
				sys.stdout.write("Warning: non-admissible graph %s appears in construction!\n" % sg)

		if self._density_graph is None:
		
			self._target_bound = c.edge_density()
			sys.stdout.write("Edge density of construction is %s.\n" % self._target_bound)
		
		else:
		
			if self._density_graph.n == self._n:
				
				# There's no guarantee that self._density_graph is minimally isomorphic!
				mdg = self._density_graph.copy()
				mdg.make_minimal_isomorph()
				
				for pair in sharp_graphs:
					g, den = pair
					if g.is_equal(mdg):
						self._target_bound = den
						break
				else:
					# if graph does not appear in construction, density is 0.
					self._target_bound = 0		

			else:
			
				self._target_bound = 0
				
				for pair in sharp_graphs:
					g, den = pair
					self._target_bound += den * g.subgraph_density(self._density_graph)
			
			sys.stdout.write("%s density of construction is %s.\n" % (self._density_graph, self._target_bound))

		
		self._zero_eigenvectors = []
		
		for ti in range(len(self._types)):
		
			if len(self._block_bases) != 0:

				row_div = self._block_bases[ti].subdivisions()[0]
				div_sizes = row_div + [len(self._flags[ti])]
				for i in range(1, len(div_sizes)):
					div_sizes[i] -= div_sizes[i - 1]
	
				M = []
				for i in range(len(div_sizes)):
					M.append(c.zero_eigenvectors(self._types[ti], self._flags[ti],
						self._block_bases[ti].subdivision(i,0)))
				self._zero_eigenvectors.append(block_diagonal_matrix(M))

			else:

				self._zero_eigenvectors.append(c.zero_eigenvectors(self._types[ti], self._flags[ti]))
		
		for ti in range(len(self._types)):
			self._zero_eigenvectors[ti].set_immutable()

	
	# TODO: make this work if _block_bases == [].

	def add_zero_eigenvectors(self, ti, bi, M):
	
		OM = self._flag_bases[ti].subdivision(bi,0).solve_right(M.T)
		MM = OM.T * self._block_bases[ti].subdivision(bi,0).T
		
		row_div = self._block_bases[ti].subdivisions()[0]
		div_sizes = row_div + [len(self._flags[ti])]
		for i in range(1, len(div_sizes)):
			div_sizes[i] -= div_sizes[i - 1]
		
		B = []
		for i in range(len(row_div) + 1):
			BL = self._zero_eigenvectors[ti].subdivision(i, i)
			if BL.nrows() == 0:
				BL = matrix(QQ, 0, div_sizes[i])
			if i == bi:
				B.append(BL.stack(MM))
			else:
				B.append(BL)

		self._zero_eigenvectors[ti] = block_diagonal_matrix(B)



	def set_new_bases(self):

		# If there are no zero eigenvectors, just use the block bases
		if len(self._zero_eigenvectors) == 0:
			for ti in range(len(self._types)):
				self._flag_bases[ti] = self._block_bases[ti]
			return

		for ti in range(len(self._types)):
	
			col_div = self._zero_eigenvectors[ti].subdivisions()[1]
			num_M = len(col_div) + 1
			div_sizes = col_div + [self._zero_eigenvectors[ti].ncols()]
			for i in range(1, len(div_sizes)):
				div_sizes[i] -= div_sizes[i - 1]

			M = [self._zero_eigenvectors[ti].subdivision(i,i) for i in range(num_M)]
	
			for i in range(num_M):
				nzev = M[i].nrows()
				if nzev == 0:
					M[i] = identity_matrix(QQ, div_sizes[i])
				else:
					M[i] = M[i].stack(M[i].right_kernel().basis_matrix())
					M[i], mu = M[i].gram_schmidt()
					M[i] = M[i][nzev:,:] # delete rows corresponding to zero eigenvectors

			if len(self._block_bases) != 0:
				self._flag_bases[ti] = block_diagonal_matrix(M) * self._block_bases[ti]
			else:
				self._flag_bases[ti] = block_diagonal_matrix(M)
			self._flag_bases[ti].set_immutable()



 	def calculate_product_densities(self):
 	
 		graph_block = make_graph_block(self._graphs, self.n)
		
		self._product_densities_dumps = []
 		sys.stdout.write("Calculating product densities...\n")

 		for ti in range(len(self._types)):

			B = self._flag_bases[ti]
			row_div = B.subdivisions()[0]
			try:
				inv_rows = row_div[0]
				is_subdivided = True
			except IndexError:
				is_subdivided = False

 			tg = self._types[ti]
 			s = tg.n
 			m = (self._n + s) / 2

 			sys.stdout.write("Doing type %d (order %d; flags %d).\n" % (ti, s, m))
 			flags_block = make_graph_block(self._flags[ti], m)
			DL = flag_products(graph_block, tg, flags_block, None)
		
			this_type_dumps = []
			for gi in range(len(self._graphs)):
				D = DL[gi]
				ND = B * D * B.T
				if is_subdivided:
					ND.subdivide(inv_rows, inv_rows)
				this_type_dumps.append(dumps(ND))
			
			self._product_densities_dumps.append(this_type_dumps)

	
	def write_blocks(self, f, write_sizes=False):

		num_graphs = len(self._graphs)
		num_types = len(self._types)
		
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

		if write_sizes:

			for i in range(num_types):
				f.write("%d " % inv_block_sizes[i])
				f.write("%d " % anti_inv_block_sizes[i])

			return

		for i in range(num_graphs):
			for j in range(num_types):
				D = loads(self._product_densities_dumps[j][i])
				for key, value in D.dict().iteritems():
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


	def write_sdp_input_file(self):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)
		
		if not self._minimize:
			# Multiply SDP solver objective value by -1
			self._obj_value_factor = -1.0
		else:
			self._obj_value_factor = 1.0
		
		if num_densities < 1:
			raise NotImplementedError("there must be at least one density.")
		
		self._sdp_input_filename = os.path.join(SAGE_TMP, "sdp.dat-s")
	
		with open(self._sdp_input_filename, "w") as f:
	
			f.write("%d\n" % (num_graphs + 1,))
			f.write("%d\n" % (2 * num_types + 3,))
			
			f.write("1 ")
			self.write_blocks(f, True)
			f.write("-%d -%d\n" % (num_graphs, num_densities))
			f.write("%s1.0\n" % ("0.0 " * num_graphs,))
			
			if not self._minimize:
				f.write("0 1 1 1 -1.0\n")
			else:
				f.write("0 1 1 1 1.0\n")
	
			for i in range(num_graphs):
				if not self._minimize:
					f.write("%d 1 1 1 -1.0\n" % (i + 1,))
				else:
					f.write("%d 1 1 1 1.0\n" % (i + 1,))
				if not (self._force_sharps and i in self._sharp_graphs):
					f.write("%d %d %d %d 1.0\n" % (i + 1, 2 * num_types + 2, i + 1, i + 1))
	
			for i in range(num_graphs):
				for j in range(num_densities):
					d = self._densities[j][i]
					if d != 0:
						if self._minimize:
							d *= -1
						f.write("%d %d %d %d %s\n" % (i + 1, 2 * num_types + 3, j + 1, j + 1, d.n(digits=64)))
			
			for j in range(num_densities):
				f.write("%d %d %d %d 1.0\n" % (num_graphs + 1, 2 * num_types + 3, j + 1, j + 1))
		
			self.write_blocks(f)


	def write_alternate_sdp_input_file(self):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)

		self._obj_value_factor = 1.0
		
		if num_densities < 1:
			raise NotImplementedError("there must be at least one density.")
		
		self._sdp_input_filename = os.path.join(SAGE_TMP, "sdp.dat-s")
	
		with open(self._sdp_input_filename, "w") as f:
	
			f.write("%d\n" % (num_graphs + num_densities * 2,))
			f.write("%d\n" % (2 * num_types + 5,))
			
			f.write("1 ")
			self.write_blocks(f, True)
			f.write("-%d -%d -%d -%d\n" % (num_graphs, num_densities, num_densities, num_densities))
			f.write("%s%s%s\n" % ("0.0 " * num_graphs, "0.0 " * num_densities, "1.0 " * num_densities))
			f.write("0 1 1 1 1.0\n")
	
			for i in range(num_graphs):
				if not (self._force_sharps and i in self._sharp_graphs):
					f.write("%d %d %d %d 1.0\n" % (i + 1, 2 * num_types + 2, i + 1, i + 1))
	
			for i in range(num_graphs):
				for j in range(num_densities):
					d = self._densities[i][j]
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



	def run_csdp(self, show_output=False):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)

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
					obj_val = RDF(line.split()[-1]) * self._obj_value_factor
				if show_output:
					sys.stdout.write(line)
			except pexpect.EOF:
				break
		p.close()
		returncode = p.exitstatus
		
		sys.stdout.write("Returncode is %d. Objective value is %s.\n" % (returncode,
			obj_val))
			
		# TODO: if program is infeasible, a returncode of 1 is given,
		# and output contains "infeasible"

		with open(self._sdp_output_filename, "r") as f:
			self.process_sdp_output_file(f)


	def process_sdp_output_file(self, f):

		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)
		
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

		self._sdp_Q_matrices = [matrix(RDF, self._flag_bases[i].nrows(), self._flag_bases[i].nrows())
			for i in range(num_types)]
		
		for ti in range(num_types):
			B = self._flag_bases[ti]
			row_div = B.subdivisions()[0]
			self._sdp_Q_matrices[ti].subdivide(row_div, row_div)
		
		self._sdp_density_coeffs = [0.0 for i in range(num_densities)]
		
		for line in f:
			numbers = line.split()
			if numbers[0] != "2":
				continue
			ti = int(numbers[1]) - 2
			if ti == 2 * num_types + 1:
				j = int(numbers[2]) - 1
				self._sdp_density_coeffs[j] = RDF(numbers[4])
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

		for ti in range(num_types):
			self._sdp_Q_matrices[ti].set_immutable()


	def show_zero_eigenvalues(self, tolerance = 0.00001):
	
		num_types = len(self._types)

		for ti in range(num_types):
			eigvals = self._sdp_Q_matrices[ti].eigenvalues()
			zero_eigvals = sorted([e for e in eigvals if e < tolerance])
			if len(zero_eigvals) == 0:
				sys.stdout.write("Type %d. None.\n" % ti)
			else:
				sys.stdout.write("Type %d. %d possible: %s.\n" % (ti, len(zero_eigvals),
					" ".join("%s" % e for e in zero_eigvals)))


	def get_zero_eigenvectors(self, ti, tolerance = 0.00001):
	
		if not ti in range(len(self._types)):
			raise ValueError

		ns = len(self._sdp_Q_matrices[ti].subdivisions()[0]) + 1
	
		B = []
		for i in range(ns):
			QB = self._sdp_Q_matrices[ti].subdivision(i, i)
			eigvals, T = numpy.linalg.eigh(QB)
			M = matrix(RDF, 0, QB.ncols())
			for ei in range(len(eigvals)):
				if eigvals[ei] < tolerance:
					M = M.stack(matrix(numpy.matrix(T[:, ei])))
			B.append(M)								
		return block_diagonal_matrix(B)


	def check_construction(self, C, tolerance = 0.00001):
	
		num_types = len(self._types)

		for ti in range(num_types):
			M = C.zero_eigenvectors(self._types[ti], self._flags[ti], self._flag_bases[ti])
			if M.nrows() == 0:
				sys.stdout.write("Type %d. None.\n" % ti)
			else:
				R = M * self._sdp_Q_matrices[ti]
				norms = [R[i,:].norm() for i in range(R.nrows())]
				sys.stdout.write("Type %d. %d possible: %s.\n" % (ti, len(norms),
					" ".join("%s" % e for e in norms)))


	def check_floating_point_bound(self, tolerance = 0.00001):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)
		
		fbounds = [sum([self._densities[j][i] * self._sdp_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]
		
		for ti in range(num_types):
			num_flags = len(self._flags[ti])
			for gi in range(num_graphs):
				D = loads(self._product_densities_dumps[ti][gi])
				for j in range(D.nrows()):
					for k in range(j, D.nrows()):
						value = self._sdp_Q_matrices[ti][j, k]
						if j != k:
							value *= 2
						d = D[j, k]
						if d != 0:
							if not self._minimize:
								fbounds[gi] += d * value
							else:
								fbounds[gi] -= d * value

		if not self._minimize:
			bound = max(fbounds)
		else:
			bound = min(fbounds)
		
		if not self._target_bound is None:
			if abs(bound - RDF(self._target_bound)) < tolerance:
				sys.stdout.write("Bound of %s appears to have been met.\n" % self._target_bound)
			else:
				sys.stdout.write("Warning: bound of %s appears to have not been met.\n" % self._target_bound)
				
		apparently_sharp_graphs = [gi for gi in range(num_graphs) if abs(fbounds[gi] - bound) < tolerance]

		sys.stdout.write("The following %d graphs appear to be sharp:\n" % len(apparently_sharp_graphs))
		for gi in apparently_sharp_graphs:
			sys.stdout.write("%s : graph %d (%s)\n" % (fbounds[gi], gi, self._graphs[gi]))
		
		extra_sharp_graphs = [gi for gi in apparently_sharp_graphs if not gi in self._sharp_graphs]
		missing_sharp_graphs = [gi for gi in self._sharp_graphs if not gi in apparently_sharp_graphs]
		
		for gi in extra_sharp_graphs:
			sys.stdout.write("Warning: graph %d (%s) appears to be sharp.\n" % (gi, self._graphs[gi]))

		for gi in missing_sharp_graphs:
			sys.stdout.write("Warning: graph %d (%s) does not appear to be sharp.\n" % (gi, self._graphs[gi]))


	def make_exact(self, denominator=1024, cholesky=[], protect=[]):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_sharps = len(self._sharp_graphs)
		num_densities = len(self._densities)

		q_sizes = [self._sdp_Q_matrices[ti].nrows() for ti in range(num_types)]

		def rationalize (f):
			return Integer(round(f * denominator)) / denominator

		self._exact_Q_matrices = []
		for ti in range(num_types):
			
			if ti in cholesky:
				try:
					LF = numpy.linalg.cholesky(self._sdp_Q_matrices[ti])
				except numpy.linalg.linalg.LinAlgError:
					sys.stdout.write("Could not compute Cholesky decomposition for type %d.\n" % ti)
					return
				L = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
				for j in range(q_sizes[ti]):
					for k in range(j + 1): # only lower triangle
						L[j, k] = rationalize(LF[j, k])
				M = L * L.T
			
			else:
				M = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
				row_div = self._flag_bases[ti].subdivisions()[0]
				M.subdivide(row_div, row_div)
				for j in range(q_sizes[ti]):
					for k in range(j, q_sizes[ti]):
						value = rationalize(self._sdp_Q_matrices[ti][j, k])
						if value != 0:
							M[j, k] = value
							M[k, j] = value
			
			self._exact_Q_matrices.append(M)

		if num_densities == 1:
			self._exact_density_coeffs = [Integer(1)]
		else:
			self._exact_density_coeffs = [rationalize(self._sdp_density_coeffs[di]) for di in range(num_densities)]
			
		triples = [(ti, j, k) for ti in range(num_types) for j in range(q_sizes[ti])
			for k in range(j, q_sizes[ti])]
	
		num_triples = len(triples)
		triples.sort()
		triple_to_index = dict((triples[i], i) for i in range(num_triples))

		R = matrix(QQ, num_sharps, num_triples, sparse=True)

		for si in range(num_sharps):
			gi = self._sharp_graphs[si]
			for ti in range(num_types):
				D = loads(self._product_densities_dumps[ti][gi])
				for j in range(q_sizes[ti]):
					for k in range(j, q_sizes[ti]):
						trip = (ti, j, k)
						value = D[j, k]
						if j != k:
							value *= 2
						if self._minimize:
							value *= -1
						R[si, triple_to_index[trip]] = value
		
		density_cols_to_use = []
		DR = matrix(QQ, num_sharps, 0, sparse=True)
		rankDR = 0
		
		# Only if there is more than one density
		if num_densities > 1:
			
			for j in range(num_densities):
				newDR = DR.augment(matrix(QQ, [[self._densities[j][gi] for gi in self._sharp_graphs]]).T)
				if newDR.rank() > rankDR:
					rankDR += 1
					DR = newDR
					density_cols_to_use.append(j)
			
			sys.stdout.write("Chosen DR matrix of rank %d.\n" % rankDR)
				
		col_norms = {}
		for i in range(num_triples):
			n = sum(x**2 for x in R.column(i))
			if n != 0:
				col_norms[i] = n

		# Use columns with greatest non-zero norm - change minus to plus to
		# use smallest columns (not sure which is best, or maybe middle?)
		
		cols_in_order = sorted(col_norms.keys(), key = lambda i : -col_norms[i])
		cols_to_use = []
	
		for i in cols_in_order:
			ti, j, k = triples[i]
			if ti in protect: # don't use protected types
				continue
			newDR = DR.augment(R[:, i : i + 1])
			if newDR.rank() > rankDR:
				rankDR += 1
				DR = newDR
				cols_to_use.append(i)
		
		sys.stdout.write("Chosen DR matrix of rank %d.\n" % rankDR)
		
		T = matrix(QQ, num_sharps, 1, sparse=True)
		
		for si in range(num_sharps):
		
			gi = self._sharp_graphs[si]
			T[si, 0] = self._target_bound
			
			for j in range(num_densities):
				if not j in density_cols_to_use:
					T[si, 0] -= self._exact_density_coeffs[j] * self._densities[j][gi]
		
			for i in range(num_triples):
				if not i in cols_to_use:
					ti, j, k = triples[i]
					T[si, 0] -= self._exact_Q_matrices[ti][j, k] * R[si, i]

		X = DR.solve_right(T)
		RX = matrix(RDF, X.nrows(), 1)
	
		for i in range(len(density_cols_to_use)):
			di = density_cols_to_use[i]
			RX[i, 0] = self._exact_density_coeffs[di]
			self._exact_density_coeffs[di] = X[i, 0]

		for i in range(len(density_cols_to_use), X.nrows()):
			ti, j, k = triples[cols_to_use[i - len(density_cols_to_use)]]
			RX[i, 0] = self._sdp_Q_matrices[ti][j, k]
			self._exact_Q_matrices[ti][j, k] = X[i, 0]
			self._exact_Q_matrices[ti][k, j] = X[i, 0]
		
		for i in range(X.nrows()):
			print "%s : %s" % (RX[i,0], RDF(X[i,0]))

		for ti in range(num_types):
			self._exact_Q_matrices[ti].set_immutable()
	
	

	def compare_eigenvalues(self):
	
		for ti in range(len(self._types)):
			sys.stdout.write("Type %d:\n" % ti)
			original_eigvals = sorted(numpy.linalg.eigvalsh(self._sdp_Q_matrices[ti]))
			new_eigvals = sorted(numpy.linalg.eigvalsh(self._exact_Q_matrices[ti]))
			for i in range(len(original_eigvals)):
				sys.stdout.write("%s : %s\n" % (original_eigvals[i], RDF(new_eigvals[i])))


	def check_exact_bound(self):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)
		
		bounds = [sum([self._densities[j][i] * self._exact_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]
		
		for ti in range(num_types):
			num_flags = len(self._flags[ti])
			for gi in range(num_graphs):
				D = loads(self._product_densities_dumps[ti][gi])
				for j in range(D.nrows()):
					for k in range(j, D.nrows()):
						value = self._exact_Q_matrices[ti][j, k]
						if j != k:
							value *= 2
						d = D[j, k]
						if d != 0:
							if not self._minimize:
								bounds[gi] += d * value
							else:
								bounds[gi] -= d * value

		if not self._minimize:
			bound = max(bounds)
			violators = [gi for gi in range(num_graphs) if bounds[gi] > self._target_bound]
		else:
			bound = min(bounds)
			violators = [gi for gi in range(num_graphs) if bounds[gi] < self._target_bound]
		
		sys.stdout.write("Bound of %s attained by:\n" % self._target_bound)
		for gi in range(num_graphs):
			if bounds[gi] == self._target_bound:
				sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))

		if len(violators) > 0:
			sys.stdout.write("Bound violated by:")
			for gi in violators:
				sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))

		self._bounds = bounds


	def combine_densities(self, denominator=32, larger_than=0.0):
	
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)
	
		def rationalize (f):
			return Integer(round(f * denominator)) / denominator
	
		new_density_weights = [rationalize(n) for n in self._sdp_density_coeffs]
	
		for j in range(num_densities):
			if self._sdp_density_coeffs[j] < larger_than:
				new_density_weights[j] = Integer(0)
	
		print new_density_weights
	
		new_density = [Integer(0) for i in range(num_graphs)]
		for i in range(num_graphs):
			for j in range(num_densities):
				new_density[i] += self._densities[j][i] * new_density_weights[j]
				
		new_problem = copy(self)
		new_problem._sdp_Q_matrices = []
		new_problem._exact_Q_matrices = []
		new_problem._sdp_density_coeffs = []
		new_problem._exact_density_coeffs = []
		new_problem._bounds = []
		new_problem._densities = [new_density]

		return new_problem


	def lose_small_densities(self, larger_than):
	
		num_densities = len(self._densities)

		new_densities = []
		for j in range(num_densities):
			if self._sdp_density_coeffs[j] > larger_than:
				new_densities.append(self._densities[j])
		
		new_problem = copy(self)
		new_problem._sdp_Q_matrices = []
		new_problem._exact_Q_matrices = []
		new_problem._sdp_density_coeffs = []
		new_problem._exact_density_coeffs = []
		new_problem._bounds = []
		new_problem._densities = new_densities

		return new_problem

