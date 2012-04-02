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

import gzip
import json
import numpy
import os
import pexpect
import sys

from sage.structure.sage_object import SageObject
from sage.rings.all import Integer, QQ, RationalField, RDF
from sage.matrix.all import matrix, identity_matrix, block_diagonal_matrix
from sage.modules.misc import gram_schmidt
from sage.misc.misc import SAGE_TMP 

from flag import *
from flag_misc import *

# pexpect in Sage has a bug, which prevents it using commands with full paths.
# So for now, CSDP has to be in a directory in $PATH.

cdsp_cmd = "csdp"
sdpa_cmd = "sdpa"
sdpa_dd_cmd = "sdpa_dd"
sdpa_qd_cmd = "sdpa_qd"

class Problem(SageObject):

	def __init__(self, r=3, oriented=False):

		# Set some defaults...
		
		self._n = 0
		self._r = r
		self._oriented = oriented
		self._field = RationalField()
		self._approximate_field = RDF

		self._forbidden_edge_numbers = []
		self._forbidden_graphs = []
		self._forbidden_induced_graphs = []
		
		edge_graph = Flag("%d:" % r, r, oriented)
		edge_graph.add_edge(range(1, r + 1))
		self._density_graphs = [edge_graph]

		self._obj_value_factor = -1
		self._minimize = False
		self._force_sharps = False

		# Set these to be empty, as length is tested
		# TODO: get rid of them if not used
		self._block_bases = []
		self._flag_bases = []
		self._solution_bases = []


	def save(self, filename):

		if filename[-5:] != ".sobj":
			filename += ".sobj"

		with open(filename, "wb") as f:
			f.write(dumps(self.__dict__))


	@classmethod
	def load(cls, filename):

		if filename[-5:] != ".sobj":
			filename += ".sobj"
		
		obj = cls()
		
		with open(filename, "rb") as f:
			d = loads(f.read())
			obj.__dict__.update(d)
		
		return obj	

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
			forbidden_edge_numbers=self._forbidden_edge_numbers,
			forbidden_graphs=self._forbidden_graphs,
			forbidden_induced_graphs=self._forbidden_induced_graphs)
		sys.stdout.write("Generated %d graphs.\n" % len(self._graphs))
	
		self.calculate_densities()
	
		sys.stdout.write("Generating types and flags...\n")
		self._types = []
		self._flags = []
	
		for s in range(n % 2, n - 1, 2):
			
			these_types = generate_graphs(s, self._r, self._oriented,
				forbidden_edge_numbers=self._forbidden_edge_numbers,
				forbidden_graphs=self._forbidden_graphs,
				forbidden_induced_graphs=self._forbidden_induced_graphs)
			sys.stdout.write("Generated %d types of order %d, " % (
				len(these_types), s))

			m = (n + s) / 2
			these_flags = []
			for tg in these_types:
				these_flags.append(generate_flags(m, tg, self._r, self._oriented,
					forbidden_edge_numbers=self._forbidden_edge_numbers,
					forbidden_graphs=self._forbidden_graphs,
					forbidden_induced_graphs=self._forbidden_induced_graphs))
			sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))
						
			self._types.extend(these_types)
			self._flags.extend(these_flags)

		#self._flag_bases = [identity_matrix(QQ, len(self._flags[ti]), sparse=True) for ti in range(len(self._types))]
		#for ti in range(len(self._types)):
		#	self._flag_bases[ti].set_immutable()


	@property
	def graphs(self):
		return self._graphs

	@property
	def types(self):
		return self._types

	@property
	def flags(self):
		return self._flags

	@property
	def density_graphs(self):
		return self._density_graphs


	def calculate_densities(self):
	
		self._densities = [[sum(g.subgraph_density(dg) for dg in self._density_graphs)
			for g in self._graphs]]
	
	
	def set_density_graph(self, dg):

		self._density_graphs = [dg]
		self.calculate_densities()


	def set_density_graphs(self, dgs):

		self._density_graphs = dgs
		self.calculate_densities()


	def set_density_edge_number(self, k, ne):

		if k < self._r:
			raise ValueError

		max_e = binomial(k, self._r)
		if not ne in range(max_e + 1):
			raise ValueError

		graphs = generate_graphs(k, self._r, self._oriented,
			forbidden_edge_numbers=self._forbidden_edge_numbers,
			forbidden_graphs=self._forbidden_graphs,
			forbidden_induced_graphs=self._forbidden_induced_graphs)
		
		dgs = [g for g in graphs if g.ne == ne]

		self.set_density_graphs(dgs)


	def forbid_edge_number(self, k, ne):

		if k < self._r: 
			raise ValueError

		max_e = binomial(k, self._r)
		if not ne in range(max_e + 1):
			raise ValueError

		for i in range(ne, max_e + 1):
			self._forbidden_edge_numbers.append((k, i))


	def forbid_induced_edge_number(self, k, ne):

		if k < self._r:
			raise ValueError

		max_e = binomial(k, self._r)
		if not ne in range(max_e + 1):
			raise ValueError
	
		self._forbidden_edge_numbers.append((k, ne))


	def forbid_subgraph(self, h):

		self._forbidden_graphs.append(h)


	def forbid_induced_subgraph(self, h):

		self._forbidden_induced_graphs.append(h)


	# TODO: warn if already solved

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


	def use_construction(self, c):

		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)

		self._field = c.field

		# Ad Hoc constructions use target_bound() instead of subgraph_densities()
		pair = c.target_bound()

		if not pair is None:
		
			self._target_bound = pair[0]
			self._sharp_graphs = pair[1]

		else:

			sys.stdout.write("Determining which graphs appear in construction...\n")
			
			sharp_graphs = c.subgraph_densities(self._n)
			target_densities = [0 for j in range(num_densities)]
			self._sharp_graphs = []
			
			for pair in sharp_graphs:
				g, den = pair
				for gi in range(num_graphs):
					if g.is_equal(self._graphs[gi]):
						self._sharp_graphs.append(gi)
						for j in range(num_densities):
							target_densities[j] += self._densities[j][gi] * den
						break
				else:
					sys.stdout.write("Warning: non-admissible graph %s appears in construction!\n" % g)
	
			# set target_bound to equal the maximum - probably this will always be what is wanted...
			self._target_bound = max(target_densities)
			
		sys.stdout.write("Density of construction is %s.\n" % self._target_bound)
		
		self._zero_eigenvectors = []
		
		for ti in range(len(self._types)):

			self._zero_eigenvectors.append(c.zero_eigenvectors(self._types[ti], self._flags[ti]))
		
			sys.stdout.write("Found %d zero eigenvectors for type %d.\n" % (
				self._zero_eigenvectors[ti].nrows(), ti))
		
		for ti in range(len(self._types)):
			self._zero_eigenvectors[ti].set_immutable()

	

	# TODO: reinstate the per-block option?

	def add_zero_eigenvectors(self, ti, M):

		if len(self._flag_bases) > 0:
			NZ = (self._flag_bases[ti].T).solve_left(M)
		else:
			NZ = M

		self._zero_eigenvectors[ti] = self._zero_eigenvectors[ti].stack(NZ)
		self._zero_eigenvectors[ti].set_immutable()


	def change_problem_bases(self, use_blocks=True):

		self._flag_bases = self._create_new_bases(use_blocks)

	
	def change_solution_bases(self, use_blocks=True, use_smaller=False):

		self._solution_bases = self._create_new_bases(use_blocks)
		self._inverse_solution_bases = []

		for ti in range(len(self._types)):
			M = copy(self._solution_bases[ti])
			if use_smaller:
				self._inverse_solution_bases.append(M.T)
				for j in range(M.nrows()):
					M[j, :] /= sum([x**2 for x in M.row(j)])
				self._solution_bases[ti] = M
			else:
				for j in range(M.nrows()):
					M[j, :] /= sum([x**2 for x in M.row(j)])
				self._inverse_solution_bases.append(M.T)


	def create_block_bases(self):

		self._block_bases = []
		for ti in range(len(self._types)):
			B = flag_basis(self._types[ti], self._flags[ti])
			row_div = B.subdivisions()[0]
			div_sizes = row_div + [len(self._flags[ti])]
			for bi in range(1, len(div_sizes)):
				div_sizes[bi] -= div_sizes[bi - 1]
			sys.stdout.write("Type %d (%d flags) blocks: %s \n" % (ti, len(self._flags[ti]), div_sizes))
			self._block_bases.append(B)


	def _create_new_bases(self, use_blocks=True, keep_rows=False):

		num_types = len(self._types)

		if len(self._zero_eigenvectors) == 0:

			sys.stdout.write("No zero eigenvectors found.\n")
			
			if use_blocks:
				return [self._block_bases[ti] for ti in range(num_types)]
			else:
				return [identity_matrix(QQ, len(self._flags[ti]), sparse=True) for ti in range(num_types)]

		if use_blocks and len(self._block_bases) == 0:
			self.create_block_bases()

		new_bases = []

		for ti in range(num_types):
	
			if use_blocks:
			
				row_div = self._block_bases[ti].subdivisions()[0]
				div_sizes = row_div + [len(self._flags[ti])]
				for bi in range(1, len(div_sizes)):
					div_sizes[bi] -= div_sizes[bi - 1]
			
			else:
		
				div_sizes = [len(self._flags[ti])]
				
			BS = []
			
			for bi in range(len(div_sizes)):	
			
				if use_blocks:
					B = (self._block_bases[ti].subdivision(bi, 0) * self._zero_eigenvectors[ti].T).T
				else:
					B = self._zero_eigenvectors[ti]
				
				B = B.echelon_form()

				nzev = B.rank()
				B = B[:nzev, :]

				if nzev == 0:
					B = identity_matrix(QQ, div_sizes[bi], sparse=True)
				
				elif nzev == div_sizes[bi]:
					pass
					
				else:
					B = B.stack(B.right_kernel().basis_matrix())

				if use_blocks:
					B = B * self._block_bases[ti].subdivision(bi, 0)

				if not keep_rows:
					B = B[nzev:, :] # delete rows corresponding to zero eigenvectors

				if B.nrows() > 0:
					BS.append(B)

			M = block_matrix([[B] for B in BS], subdivide=True)
			div = M.subdivisions()

			if self._field == RationalField():
				
				M, mu = M.gram_schmidt()
				# .gram_schmidt doesn't appear to preserve sparsity
				M = matrix(self._field, M.rows(), sparse=True)

			else: # .gram_schmidt is broken for number fields in 4.8 !
				rows, mu = gram_schmidt(M.rows())
				M = matrix(self._field, rows, sparse=True)
			
			M.subdivide(div)
			
			M.set_immutable()
			new_bases.append(M)
	
		return new_bases


	def calculate_product_densities(self, compress=False):
 	
		graph_block = make_graph_block(self._graphs, self.n)
		
		self._product_densities_arrays = []
 		
		sys.stdout.write("Calculating product densities...\n")

		for ti in range(len(self._types)):

			tg = self._types[ti]
			s = tg.n
			m = (self._n + s) / 2
			nf = len(self._flags[ti])

			sys.stdout.write("Doing type %d (order %d; flags %d)...\n" % (ti, s, m))

			flags_block = make_graph_block(self._flags[ti], m)
			rarray = flag_products(graph_block, tg, flags_block, None)
			row_num = 0

			if len(self._flag_bases) > 0:

				sys.stdout.write("Transforming...\n")

				B = self._flag_bases[ti]
				new_rarray = numpy.zeros([0, 5], dtype=numpy.int)

				for gi in range(len(self._graphs)):
				
					D = matrix(QQ, nf, nf, sparse=True)
					for row in rarray:
						if row[0] == gi:
							j = row[1]
							k = row[2]
							value = Integer(row[3]) / Integer(row[4])
							D[j, k] = value
							D[k, j] = value
					ND = B * D * B.T
					# Only want upper triangle
					NDD = dict((pair, value) for pair, value in ND.dict().iteritems()
						if pair[0] <= pair[1])
					new_rarray.resize([row_num + len(NDD), 5], refcheck=False)
					for pair, value in NDD.iteritems():
						new_rarray[row_num, 0] = gi
						new_rarray[row_num, 1] = pair[0]
						new_rarray[row_num, 2] = pair[1]
						new_rarray[row_num, 3] = value.numerator()
						new_rarray[row_num, 4] = value.denominator()
						row_num += 1
				
				self._product_densities_arrays.append(new_rarray)

			else:
				self._product_densities_arrays.append(rarray)
	

	def write_blocks(self, f, write_sizes=False):

		num_graphs = len(self._graphs)
		num_types = len(self._types)
		
		inv_block_sizes = []
		anti_inv_block_sizes = []
		
		if len(self._flag_bases) > 0:
			for ti in range(num_types):
				row_div = self._flag_bases[ti].subdivisions()[0]
				if len(row_div) > 0:
					inv_block_sizes.append(row_div[0])
					anti_inv_block_sizes.append(self._flag_bases[ti].nrows() - row_div[0])
				else:
					inv_block_sizes.append(self._flag_bases[ti].nrows())
					anti_inv_block_sizes.append(1)
		else:
			for ti in range(num_types):
				inv_block_sizes.append(len(self._flags[ti]))
				anti_inv_block_sizes.append(1)

		if write_sizes:

			for i in range(num_types):
				f.write("%d " % inv_block_sizes[i])
				f.write("%d " % anti_inv_block_sizes[i])

			return

		for ti in range(num_types):
			for row in self._product_densities_arrays[ti]:
				gi = row[0]
				j = row[1]
				k = row[2]
				if j < inv_block_sizes[ti]:
					block = 2 * ti + 2
				else:
					block = 2 * ti + 3
					j -= inv_block_sizes[ti]
					k -= inv_block_sizes[ti]
				value = Integer(row[3]) / Integer(row[4])
				f.write("%d %d %d %d %s\n" % (gi + 1, block, j + 1, k + 1,
						value.n(digits=64)))

	# TODO: helpful error message if product densities have not been computed.
	
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
		
		sys.stdout.write("Writing SDP input file...\n")
		
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



	def run_sdp_solver(self, show_output=False, sdpa=False):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)

		self._sdp_output_filename = os.path.join(SAGE_TMP, "sdp.out")

		if not sdpa:

			cmd = "%s %s %s" % (cdsp_cmd, self._sdp_input_filename, self._sdp_output_filename)

		else:
		
			if sdpa == "dd":
				solver_cmd = sdpa_dd_cmd
			elif sdpa == "qd":
				solver_cmd = sdpa_qd_cmd
			else:
				solver_cmd = sdpa_cmd
		
			sdpa_output_filename = os.path.join(SAGE_TMP, "sdpa.out")
			cmd = "%s -ds %s -o %s" % (solver_cmd, self._sdp_input_filename, sdpa_output_filename)

		if not sdpa:
			sys.stdout.write("Running csdp...\n")
		else:
			sys.stdout.write("Running sdpa...\n")

		p = pexpect.spawn(cmd, timeout=60*60*24*7)
		obj_val = None
		while True:
			if p.eof():
				break
			try:
				p.expect("\r\n")
				line = p.before.strip() + "\n"
				if "Primal objective value:" in line: # CSDP
					obj_val = self._approximate_field(line.split()[-1]) * self._obj_value_factor
				elif "objValPrimal" in line: # SDPA
					obj_val = self._approximate_field(line.split()[-1]) * self._obj_value_factor
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

		if sdpa:
		
			with open(sdpa_output_filename, "r") as inf:
				with open(self._sdp_output_filename, "w") as f:

					found, diagonal = False, False
					t, row, col = 0, 1, 1
					
					for line in inf:
						if line[:6] == "yMat =":
							break
					else:
						raise ValueError
	
					for line in inf:
							
						if line == "}":
							break						
						elif line[:3] == "{ {":
							t += 1
							row = 1
							diagonal = False
						elif line[:2] == "{+" or line[:2] == "{-":
							t += 1
							row = 1
							diagonal = True	
						
						line = line.replace("{", "")
						line = line.replace("}", "")
						col = 1
						for a in line.split(","):
							try:
								v = a.strip()
								vf = float(v) # only done to see if we get ValueError
								if diagonal:
									f.write("2 %d %d %d %s\n" % (t, row, col, v))
									row += 1
								elif row <= col:
									f.write("2 %d %d %d %s\n" % (t, row, col, v))
								col += 1
							except ValueError:
								pass
						
						if col > 1: # at least one number found...
							row += 1

		with open(self._sdp_output_filename, "r") as f:
		
			inv_block_sizes = []
			anti_inv_block_sizes = []

			if len(self._flag_bases) > 0:

				for ti in range(num_types):
					row_div = self._flag_bases[ti].subdivisions()[0]
					if len(row_div) > 0:
						inv_block_sizes.append(row_div[0])
						anti_inv_block_sizes.append(self._flag_bases[ti].nrows() - row_div[0])
					else:
						inv_block_sizes.append(self._flag_bases[ti].nrows())
						anti_inv_block_sizes.append(1)

				self._sdp_Q_matrices = [matrix(self._approximate_field,
					self._flag_bases[ti].nrows(), self._flag_bases[ti].nrows())
					for ti in range(num_types)]
				
				for ti in range(num_types):
					B = self._flag_bases[ti]
					row_div = B.subdivisions()[0]
					self._sdp_Q_matrices[ti].subdivide(row_div, row_div)

			else:

				for ti in range(num_types):
					inv_block_sizes.append(len(self._flags[ti]))
					anti_inv_block_sizes.append(1)
	
				self._sdp_Q_matrices = [matrix(self._approximate_field,
					len(self._flags[ti]), len(self._flags[ti])) for ti in range(num_types)]
			
			self._sdp_density_coeffs = [self._approximate_field(0) for i in range(num_densities)]
			
			for line in f:
				numbers = line.split()
				if numbers[0] != "2":
					continue
				ti = int(numbers[1]) - 2
				if ti == 2 * num_types + 1:
					j = int(numbers[2]) - 1
					self._sdp_density_coeffs[j] = self._approximate_field(numbers[4])
					continue
				if ti < 0 or ti >= 2 * num_types:
					continue
				j = int(numbers[2]) - 1
				k = int(numbers[3]) - 1
				if ti % 2:
					ti = (ti - 1) / 2
					j += inv_block_sizes[ti]
					k += inv_block_sizes[ti]
					if j >= self._sdp_Q_matrices[ti].nrows():
						continue # Might be 'dummy' block for types with no anti-inv space
				else:
					ti /= 2
				self._sdp_Q_matrices[ti][j, k] = numbers[4]
				self._sdp_Q_matrices[ti][k, j] = self._sdp_Q_matrices[ti][j, k]

		for ti in range(num_types):
			self._sdp_Q_matrices[ti].set_immutable()


	def import_solution(self, directory):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)

		if num_densities > 1:
			raise NotImplementedError

		sys.path.insert(0, directory)
		dont_write_bytecode = sys.dont_write_bytecode 
		sys.dont_write_bytecode = True
		
		try:
			import flags
		except ImportError:
			sys.stdout.write("Cannot find flags.py in directory provided.\n")
			return

		# TODO: check admissible graphs, target bound, (sharps?).

		if not ("%d-graph" % self.r) in flags.description:
			raise ValueError
	
		if flags.n != self.n:
			raise ValueError

		if flags.num_H != num_graphs:
			raise ValueError

		if flags.num_types != num_types:
			raise ValueError

		ttr = []
		ftrs = []
	
		for ti in range(num_types):
			tg = Flag(flags.types[ti], self.r, self.oriented)
			for tj in range(num_types):
				if tg.is_equal(self._types[tj]):
					ttr.append(tj)
					break
			else:
				raise ValueError
		
			num_flags = len(self._flags[tj])
			ftr = []
			
			for fi in range(num_flags):
				fg = Flag(flags.flags[ti][fi], self.r, self.oriented)
				fg.t = tg.n
				for fj in range(num_flags):
					if fg.is_equal(self._flags[tj][fj]):
						ftr.append(fj)
						break
				else:
					print ti
					print ftr
					raise ValueError
					
			ftrs.append(ftr)
		
		print ttr
		print ftrs
		
		self._sdp_Q_matrices = [matrix(self._approximate_field, len(self._flags[ti]),
			len(self._flags[ti])) for ti in range(num_types)]

		try:
			f = open(directory + "/" + flags.out_filename, "r")
		except IOError:
			try:
				f = gzip.open(directory + "/" + flags.out_filename + ".gz", "rb")
			except IOError:
				print "Could not open %s or %s.gz" % (flags.out_filename,
					flags.out_filename)
				return

		for line in f:
			numbers = line.split()
			if numbers[0] != "2":
				continue
			ti = int(numbers[1]) - 2
			if ti >= 0 and ti < num_types:
				tj = ttr[ti]
				j = ftrs[ti][int(numbers[2]) - 1]
				k = ftrs[ti][int(numbers[3]) - 1]
				self._sdp_Q_matrices[tj][j, k] = numbers[4]
				self._sdp_Q_matrices[tj][k, j] = self._sdp_Q_matrices[tj][j, k]

		f.close()

		self._sdp_density_coeffs = [1.0]

		for ti in range(num_types):
			self._sdp_Q_matrices[ti].set_immutable()

		sys.path.remove(directory)
		sys.dont_write_bytecode = dont_write_bytecode


	def show_zero_eigenvalues(self, tolerance = 0.00001):
	
		num_types = len(self._types)

		for ti in range(num_types):
			eigvals = numpy.linalg.eigvalsh(self._sdp_Q_matrices[ti])
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
			M = matrix(self._approximate_field, 0, QB.ncols())
			for ei in range(len(eigvals)):
				if eigvals[ei] < tolerance:
					M = M.stack(matrix(numpy.matrix(T[:, ei])))
			B.append(M)
		return block_diagonal_matrix(B)


	def check_construction(self, C, tolerance = 0.00001):
	
		# TODO: transform with _flag_bases if present.
		
		num_types = len(self._types)

		for ti in range(num_types):
			M = C.zero_eigenvectors(self._types[ti], self._flags[ti]) 
			if M.nrows() == 0:
				sys.stdout.write("Type %d. None.\n" % ti)
			else:
				R = M * self._sdp_Q_matrices[ti]
				norms = [R[i,:].norm() for i in range(R.nrows())]
				sys.stdout.write("Type %d. %d possible: %s.\n" % (ti, len(norms),
					" ".join("%s" % e for e in norms)))



	def check_floating_point_bound(self, tolerance = 0.00001, show_all=False):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)
		
		sys.stdout.write("Checking numerical bound...\n")
		
		fbounds = [sum([self._densities[j][i] * self._sdp_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]

		for ti in range(num_types):
			for row in self._product_densities_arrays[ti]:
				gi, j, k, numer, denom = row
				d = Integer(numer) / Integer(denom)
				value = self._sdp_Q_matrices[ti][j, k]
				if j != k:
					d *= 2
				if not self._minimize:
					fbounds[gi] += d * value
				else:
					fbounds[gi] -= d * value

		if not self._minimize:
			bound = max(fbounds)
		else:
			bound = min(fbounds)
		
		if not self._target_bound is None:
			if abs(bound - self._approximate_field(self._target_bound)) < tolerance:
				sys.stdout.write("Bound of %s appears to have been met.\n" % self._target_bound)
			else:
				sys.stdout.write("Warning: bound of %s appears to have not been met.\n" % self._target_bound)
				
		apparently_sharp_graphs = [gi for gi in range(num_graphs) if abs(fbounds[gi] - bound) < tolerance]

		if show_all:
		
			sorted_indices = sorted(range(num_graphs), key = lambda i : fbounds[i])

			for gi in sorted_indices:
				sys.stdout.write("%s : graph %d (%s) " % (fbounds[gi], gi, self._graphs[gi]))
				if gi in self._sharp_graphs:
					sys.stdout.write("S")
				if gi in apparently_sharp_graphs:
					sys.stdout.write("*")
				sys.stdout.write("\n")	

		else:
	
			sys.stdout.write("The following %d graphs appear to be sharp:\n" % len(apparently_sharp_graphs))
			for gi in apparently_sharp_graphs:
				sys.stdout.write("%s : graph %d (%s)\n" % (fbounds[gi], gi, self._graphs[gi]))
			
		extra_sharp_graphs = [gi for gi in apparently_sharp_graphs if not gi in self._sharp_graphs]
		missing_sharp_graphs = [gi for gi in self._sharp_graphs if not gi in apparently_sharp_graphs]
			
		#for gi in extra_sharp_graphs:
		#	sys.stdout.write("Warning: graph %d (%s) appears to be sharp.\n" % (gi, self._graphs[gi]))
	
		if len(extra_sharp_graphs) > 0:
			sys.stdout.write("Warning: additional sharp graphs: %s\n" % (extra_sharp_graphs,))	
	
		for gi in missing_sharp_graphs:
			sys.stdout.write("Warning: graph %d (%s) does not appear to be sharp.\n" % (gi, self._graphs[gi]))


	def make_exact(self, denominator=1024, cholesky=[], protect=[], show_changes=False,
		use_densities=True):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_sharps = len(self._sharp_graphs)
		num_densities = len(self._densities)

		if len(self._solution_bases) > 0:

			sys.stdout.write("Transforming matrices")
		
			self._sdp_Qdash_matrices = []

			for ti in range(num_types):
				
				B = self._solution_bases[ti]
				row_div = B.subdivisions()[0]
				M = B * self._sdp_Q_matrices[ti] * B.T
				M.subdivide(row_div, row_div)
				# zero out bits that should be zero. Note the copy() seems to be needed.
				M = block_diagonal_matrix([copy(M.subdivision(i,i)) for i in range(len(row_div) + 1)])
				M.set_immutable()
				self._sdp_Qdash_matrices.append(M)
				sys.stdout.write(".")
				sys.stdout.flush()
		
			sys.stdout.write("\n")
		
		else:
		
			self._sdp_Qdash_matrices = self._sdp_Q_matrices

		q_sizes = [self._sdp_Qdash_matrices[ti].nrows() for ti in range(num_types)]

		def rationalize (f):
			return Integer(round(f * denominator)) / denominator

		sys.stdout.write("Rounding matrices")

		self._exact_Qdash_matrices = []
		for ti in range(num_types):
			
			if ti in cholesky:
				try:
					LF = numpy.linalg.cholesky(self._sdp_Qdash_matrices[ti])
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
				for j in range(q_sizes[ti]):
					for k in range(j, q_sizes[ti]):
						value = rationalize(self._sdp_Qdash_matrices[ti][j, k])
						if value != 0:
							M[j, k] = value
							M[k, j] = value
			
			row_div = self._sdp_Qdash_matrices[ti].subdivisions()[0]
			M.subdivide(row_div, row_div)
			self._exact_Qdash_matrices.append(matrix(self._field, M))
			sys.stdout.write(".")
			sys.stdout.flush()

		sys.stdout.write("\n")

		if num_densities == 1:
			self._exact_density_coeffs = [Integer(1)]
		else:
			self._exact_density_coeffs = [rationalize(self._sdp_density_coeffs[di]) for di in range(num_densities)]
			
		triples = [(ti, j, k) for ti in range(num_types) for j in range(q_sizes[ti])
			for k in range(j, q_sizes[ti])]
	
		num_triples = len(triples)
		triples.sort()
		triple_to_index = dict((triples[i], i) for i in range(num_triples))

		R = matrix(self._field, num_sharps, num_triples, sparse=True)

		sys.stdout.write("Constructing R matrix")

		# TODO: only use triples that correspond to middle blocks.

		for ti in range(num_types):

			Ds = [matrix(QQ, len(self._flags[ti]), len(self._flags[ti]))
				for si in range(num_sharps)]

			for row in self._product_densities_arrays[ti]:
				gi = row[0]
				if not gi in self._sharp_graphs:
					continue
				si = self._sharp_graphs.index(gi)
				j = row[1]
				k = row[2]
				value = Integer(row[3]) / Integer(row[4])
				Ds[si][j, k] = value
				Ds[si][k, j] = value

			if len(self._solution_bases) > 0:
				B = self._inverse_solution_bases[ti]
				for si in range(num_sharps):
					Ds[si] = B.T * Ds[si] * B

			for si in range(num_sharps):
				for j in range(q_sizes[ti]):
					for k in range(j, q_sizes[ti]):
						trip = (ti, j, k)
						value = Ds[si][j, k]
						if j != k:
							value *= 2
						if self._minimize:
							value *= -1
						R[si, triple_to_index[trip]] = value

			sys.stdout.write(".")
			sys.stdout.flush()
		sys.stdout.write("\n")
		
		density_cols_to_use = []
		DR = matrix(self._field, 0, num_sharps, sparse=True)
		EDR = matrix(self._field, 0, num_sharps, sparse=True)
		
		sys.stdout.write("Constructing DR matrix")
		
		# Only if there is more than one density
		if num_densities > 1 and use_densities:
			
			for j in range(num_densities):
				new_row = matrix(QQ, [[self._densities[j][gi] for gi in self._sharp_graphs]], sparse=True)
				if new_row.is_zero():
					continue
				try:
					X = EDR.solve_left(new_row)
					continue
				except ValueError:
					DR = DR.stack(new_row)
					EDR = EDR.stack(new_row)
					EDR.echelonize()
					density_cols_to_use.append(j)
					sys.stdout.write(".")
					sys.stdout.flush()
			
			sys.stdout.write("\n")
			sys.stdout.write("DR matrix (density part) has rank %d.\n" % DR.nrows())
				
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
			new_row = R[:, i : i + 1].T	
			if new_row.is_zero():
				continue
			try:
				X = EDR.solve_left(new_row)
				continue
			except ValueError:
				DR = DR.stack(new_row)
				EDR = EDR.stack(new_row)
				EDR.echelonize()
				cols_to_use.append(i)
				sys.stdout.write(".")
				sys.stdout.flush()
				# TODO: add this check to density cols loop
				if DR.nrows() == num_sharps:
					sys.stdout.write(" got enough.")
					break
		
		sys.stdout.write("\n")
		sys.stdout.write("DR matrix has rank %d.\n" % DR.nrows())
		
		T = matrix(self._field, num_sharps, 1)
				
		for si in range(num_sharps):
		
			gi = self._sharp_graphs[si]
			T[si, 0] = self._target_bound
			
			for j in range(num_densities):
				if not j in density_cols_to_use:
					T[si, 0] -= self._exact_density_coeffs[j] * self._densities[j][gi]
		
			for i in range(num_triples):
				if not i in cols_to_use:
					ti, j, k = triples[i]
					T[si, 0] -= self._exact_Qdash_matrices[ti][j, k] * R[si, i]

		FDR = matrix(self._field, DR.T)
		X = FDR.solve_right(T)
		RX = matrix(self._approximate_field, X.nrows(), 1)
	
		for i in range(len(density_cols_to_use)):
			di = density_cols_to_use[i]
			RX[i, 0] = self._exact_density_coeffs[di]
			self._exact_density_coeffs[di] = X[i, 0]

		for i in range(len(density_cols_to_use), X.nrows()):
			ti, j, k = triples[cols_to_use[i - len(density_cols_to_use)]]
			RX[i, 0] = self._sdp_Qdash_matrices[ti][j, k]
			self._exact_Qdash_matrices[ti][j, k] = X[i, 0]
			self._exact_Qdash_matrices[ti][k, j] = X[i, 0]
		
		if show_changes:
			for i in range(X.nrows()):
				sys.stdout.write("%.11s -> %.11s " % (RX[i,0], RDF(X[i,0])))
				if i < len(density_cols_to_use):
					sys.stdout.write("(density %d)\n" % i)
				else:
					sys.stdout.write("(matrix %d, entry [%d, %d])\n" % triples[cols_to_use[i - len(density_cols_to_use)]])

		for ti in range(num_types):
			self._exact_Qdash_matrices[ti].set_immutable()
	

	def compare_eigenvalues(self, only_smallest=True):
	
		for ti in range(len(self._types)):
			sys.stdout.write("Type %d:\n" % ti)
			original_eigvals = sorted(numpy.linalg.eigvalsh(self._sdp_Qdash_matrices[ti]))
			new_eigvals = sorted(numpy.linalg.eigvalsh(self._exact_Qdash_matrices[ti]))
			if only_smallest:
				num = 1
			else:
				num = len(original_eigvals)
			for i in range(num):
				sys.stdout.write("%.11f : %.11f\n" % (original_eigvals[i], new_eigvals[i]))


	def check_exact_bound(self):
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)
		
		self._exact_Q_matrices = []
		
		if len(self._solution_bases) > 0:
			for ti in range(num_types):
				B = self._inverse_solution_bases[ti]
				M = B * self._exact_Qdash_matrices[ti] * B.T
				self._exact_Q_matrices.append(M)
		else:
			self._exact_Q_matrices = self._exact_Qdash_matrices
		
		bounds = [sum([self._densities[j][i] * self._exact_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]

		for ti in range(num_types):
			for row in self._product_densities_arrays[ti]:
				gi, j, k, numer, denom = row
				d = Integer(numer) / Integer(denom)
				value = self._exact_Q_matrices[ti][j, k]
				if j != k:
					value *= 2
				if not self._minimize:
					bounds[gi] += d * value
				else:
					bounds[gi] -= d * value

		if self._field == RationalField():

			if not self._minimize:
				bound = max(bounds)
				violators = [gi for gi in range(num_graphs) if bounds[gi] > self._target_bound]
			else:
				bound = min(bounds)
				violators = [gi for gi in range(num_graphs) if bounds[gi] < self._target_bound]

		else:

			# Sorting doesn't currently work for number fields with embeddings, so use float approximation.

			if not self._minimize:
				bound = max(bounds, key = lambda x : float(x))
				violators = [gi for gi in range(num_graphs) if float(bounds[gi]) > float(self._target_bound)]
			else:
				bound = min(bounds, key = lambda x : float(x))
				violators = [gi for gi in range(num_graphs) if float(bounds[gi]) < float(self._target_bound)]
			
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




	def get_large_densities(self, larger_than=0.0):

		num_densities = len(self._densities)

		densities_to_use = []
		for j in range(num_densities):
			if self._sdp_density_coeffs[j] > larger_than:
				densities_to_use.append(j)

		sys.stdout.write("Densities: %s\n" % (densities_to_use,))

		sys.stdout.write("Coefficients: %s\n" % ([self._sdp_density_coeffs[j] for j in densities_to_use],))


	def get_independent_densities(self):
	
		num_sharps = len(self._sharp_graphs)
		num_densities = len(self._densities)
	
		densities_to_use = []
		
		if len(self._sdp_density_coeffs) > 0:
			density_indices = sorted(range(num_densities), key = lambda i : -self._sdp_density_coeffs[i])
		else:
			density_indices = range(num_densities)
		
		DR = matrix(self._field, 0, num_sharps, sparse=True)
		EDR = matrix(self._field, 0, num_sharps, sparse=True)
				
		sys.stdout.write("Constructing DR matrix")
		
		for j in density_indices:
			new_row = matrix(QQ, [[self._densities[j][gi] for gi in self._sharp_graphs]], sparse=True)
			if new_row.is_zero():
				continue
			try:
				X = EDR.solve_left(new_row)
				continue
			except ValueError:
				DR = DR.stack(new_row)
				EDR = EDR.stack(new_row)
				EDR.echelonize()
				densities_to_use.append(j)
				sys.stdout.write(".")
				sys.stdout.flush()
			
		sys.stdout.write("\n")
		sys.stdout.write("Rank is %d.\n" % DR.nrows())

		sys.stdout.write("Densities: %s\n" % (densities_to_use,))

		sys.stdout.write("Coefficients: %s\n" % ([self._sdp_density_coeffs[j] for j in densities_to_use],))


	def problem_with_densities(self, densities_to_use):
	
		if len(densities_to_use) == 0:
			raise ValueError
	
		num_densities = len(self._densities)

		new_densities = []
		for j in densities_to_use:
			new_densities.append(self._densities[j])
		
		new_problem = copy(self)
		new_problem._sdp_Q_matrices = []
		new_problem._sdp_Qdash_matrices = []
		new_problem._exact_Q_matrices = []
		new_problem._exact_Qdash_matrices = []
		new_problem._sdp_density_coeffs = []
		new_problem._exact_density_coeffs = []
		new_problem._bounds = []
		new_problem._densities = new_densities

		return new_problem


	def diagonalize(self):
		
		def LDLdecomposition(M):
	
			MS = M.parent()
			D = MS.matrix()
			L = copy(MS.identity_matrix())
			for i in xrange(M.nrows()):
				for j in xrange(i):
					L[i, j] = (Integer(1) / D[j, j]) * (M[i, j] - sum(L[i, k] * L[j, k] * D[k, k] for k in xrange(j)))
				D[i, i] = M[i, i] - sum(L[i, k]**2 * D[k, k]
					for k in xrange(i))
			return L, D

		self._exact_diagonal_matrices = []
		self._exact_r_matrices = []

		sys.stdout.write("Diagonalizing")

		for ti in range(len(self._types)):
			R, M = LDLdecomposition(self._exact_Qdash_matrices[ti])
			self._exact_diagonal_matrices.append(M)
			self._exact_r_matrices.append(R)			
			sys.stdout.write(".")
			sys.stdout.flush()
		
		sys.stdout.write("\n")
		
