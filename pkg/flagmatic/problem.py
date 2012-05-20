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
from sage.functions.other import floor
from sage.matrix.all import matrix, identity_matrix, block_matrix, block_diagonal_matrix
from sage.modules.misc import gram_schmidt
from sage.misc.misc import SAGE_TMP 
from copy import copy

from hypergraph_flag import make_graph_block
from three_graph_flag import *
from graph_flag import *
from oriented_graph_flag import *
from multigraph_flag import *

# pexpect in Sage has a bug, which prevents it using commands with full paths.
# So for now, CSDP has to be in a directory in $PATH.

cdsp_cmd = "csdp"
sdpa_cmd = "sdpa"
sdpa_dd_cmd = "sdpa_dd"
sdpa_qd_cmd = "sdpa_qd"


def block_structure(M):
	"""
	Returns a tuple. The first entry is the number of blocks, the second is a list of
	the block sizes, and the third is a list containing the row in which each block begins.
	Note that only the row subdivisions are looked at.
	"""
	row_div = M.subdivisions()[0]
	div_offsets = [0] + row_div
	div_sizes = row_div + [M.nrows()]
	num_blocks = len(div_sizes)
	
	for i in range(1, num_blocks):
		div_sizes[i] -= div_sizes[i - 1]

	return num_blocks, div_sizes, div_offsets


def safe_gram_schmidt(M):
	"""
	In Sage 4.8, .gram_schmidt() is broken for matrices over number fields. This function
	will use the new implementation for matrices over the rationals, and the old
	implementation in all other cases. In addition, the sparsity and subdivisions are
	preserved.
	"""

	div = M.subdivisions()
	BR = M.base_ring()

	if BR == QQ:
		M, mu = M.gram_schmidt()
		# .gram_schmidt doesn't appear to preserve sparsity
		M = matrix(QQ, M.rows(), sparse=True)

	else: # .gram_schmidt is broken for number fields in 4.8 !
		rows, mu = gram_schmidt(M.rows())
		M = matrix(BR, rows, sparse=True)

	M.subdivide(div)
	return M


class Problem(SageObject):


	def __init__(self, flag_cls):
	
		self._flagmatic_version = "2.0"
	
		if flag_cls in [ThreeGraphFlag, GraphFlag, OrientedGraphFlag, TwoMultigraphFlag, ThreeMultigraphFlag]:
			self._flag_cls = flag_cls
		else:
			raise ValueError
		
		self._density_graphs = flag_cls.default_density_graphs()
		self._n = 0

		self._field = QQ
		self._approximate_field = RDF

		self._forbidden_edge_numbers = []
		self._forbidden_graphs = []
		self._forbidden_induced_graphs = []
		
		self._obj_value_factor = -1
		self._minimize = False
		self._force_sharps = False

		self._register_progression("specify", "set")
		

	def _register_progression(self, cur_sn, action):

		states = {
			"specify" : {
				"requires" : [],
				"depends" : []
			},
			"set_objective" : {
				"requires" : [],
				"depends" : []
			},
			"compute_flags" : {
				"requires" : [],
				"depends" : ["specify"]
			},
			"set_construction" : {
				"requires" : ["compute_flags"],
				"depends" : []
			},
			"add_zero_eigenvectors" : {
				"requires" : ["set_construction"],
				"depends" : []
			},
			"compute_block_bases" : {
				"requires" : ["compute_flags"],
				"depends" : []
			},
			"compute_flag_bases" : {
				"requires" : ["set_construction"],
				"depends" : ["add_zero_eigenvectors"]
			},
			"transform_problem" : {
				"requires" : ["compute_flag_bases"],
				"depends" : []
			},
			"compute_products" : {
				"requires" : ["compute_flags"],
				"depends" : []
			},			
			"set_active_types" : {
				"requires" : ["compute_flags"],
				"depends" : []
			},
			"write_sdp_input_file" : {
				"requires" : ["compute_flags"],
				"depends" : ["set_objective", "transform_problem", "set_active_types"]
			},
			"run_sdp_solver" : {
				"requires" : ["write_sdp_input_file"],
				"depends" : []
			},
			"transform_solution" : {
				"requires" : ["compute_flag_bases", "run_sdp_solver"],
				"depends" : []
			},
			"add_sharp_graphs" : {
				"requires" : ["set_construction"],
				"depends" : []
			},
			"make_exact" : {
				"requires" : ["run_sdp_solver"],
				"depends" : ["transform_solution", "add_sharp_graphs"]
			},
			"check_exact" : {
				"requires" : ["make_exact"],
				"depends" : []
			}
		}

		state_names = states.keys()

		if not hasattr(self, "_states"):
			self._states = dict((sn, "no") for sn in state_names)

		if not cur_sn in state_names:
			raise ValueError("unknown state.")

		if action == "set":
			if not all(self._states[sn] == "set" for sn in states[cur_sn]["requires"]):
				raise NotImplementedError("not ready for this yet!")
			self._states[cur_sn] = "set"	
			for sn in states:
				if cur_sn in states[sn]["depends"] and self._states[sn] == "set":
					self._register_progression(sn, "stale")
				
		elif action == "stale":
			self._states[cur_sn] = "stale"
			for sn in states:
				if cur_sn in states[sn]["depends"]:
					self._register_progression(sn, "stale")

		elif action == "ensure":
			if self._states[cur_sn] != "set":
				raise NotImplementedError("not ready for this yet!")

		elif action == "ensure-not-no":
			if self._states[cur_sn] == "no":
				raise NotImplementedError("not ready for this yet!")

		elif action == "query":
			return self._states[cur_sn]


	@property
	def flag_cls(self):
		return self._flag_cls


	@property
	def n(self):
		return self._n


	@n.setter
	def n(self, n):
		self.generate_flags(n)

	
	def generate_flags(self, n, type_orders="all", max_flags=None, compute_products=True):

		if type_orders == "all":
			type_orders = range(n % 2, n - 1, 2)
		else:
			if not all(s in range(n + 1) for s in type_orders):
				raise ValueError
			type_orders = sorted(list(set(type_orders)))

		self._register_progression("compute_flags", "set")

		self._n = n

		sys.stdout.write("Generating graphs...\n")
		self._graphs = self._flag_cls.generate_graphs(n,
			forbidden_edge_numbers=self._forbidden_edge_numbers,
			forbidden_graphs=self._forbidden_graphs,
			forbidden_induced_graphs=self._forbidden_induced_graphs)
		sys.stdout.write("Generated %d graphs.\n" % len(self._graphs))
	
		self._compute_densities()
	
		sys.stdout.write("Generating types and flags...\n")
		self._types = []
		self._flags = []
		
		for s in type_orders:
			
			these_types = self._flag_cls.generate_graphs(s,
				forbidden_edge_numbers=self._forbidden_edge_numbers,
				forbidden_graphs=self._forbidden_graphs,
				forbidden_induced_graphs=self._forbidden_induced_graphs)
			sys.stdout.write("Generated %d types of order %d, " % (
				len(these_types), s))

			m = floor((n + s) / 2)
			these_flags = []
			for tg in these_types:
				these_flags.append(self._flag_cls.generate_flags(m, tg,
					forbidden_edge_numbers=self._forbidden_edge_numbers,
					forbidden_graphs=self._forbidden_graphs,
					forbidden_induced_graphs=self._forbidden_induced_graphs))
			sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))
			
			self._types.extend(these_types)
			self._flags.extend(these_flags)

		num_types = len(self._types)

		if not max_flags is None:
			bad_indices = [i for i in range(num_types) if len(self._flags[i]) > max_flags]
			if len(bad_indices) > 0:
				good_indices = [i for i in range(num_types) if not i in bad_indices]
				self._types = [self._types[i] for i in good_indices]
				self._flags = [self._flags[i] for i in good_indices]
				sys.stdout.write("Removed types %s as they have too many flags.\n" % bad_indices)
		
		num_types = len(self._types) # may have changed!
		
		self._zero_eigenvectors = [matrix(self._field, 0, len(self._flags[ti]), sparse=True)
			for ti in range(num_types)]

		self._active_types = range(num_types)
		
		if compute_products:
			self.compute_products()


	@property
	def graphs(self):
		"""
		Read-only. A (copy) of the list of admissible graphs. Modifying the list will have
		no effect on the Problem.
		 
		"""
		return copy(self._graphs)

	@property
	def types(self):
		"""
		Read-only. A (copy) of the list of types. Modifying the list will have
		no effect on the Problem.
		 
		"""
		return copy(self._types)

	@property
	def flags(self):
		"""
		Read-only. A (copy) of the flags, as a list of lists. The flags are listed by
		type, in the same order as the types appear in the types list. Modifying the lists
		will have no effect on the Problem.
		 
		"""
		return copy(self._flags)

	@property
	def density_graphs(self):
		"""
		Read-only. A (copy) of the list of density graphs. Modifying the lists
		will have no effect on the Problem.
		 
		"""
		return copy(self._density_graphs)


	@property
	def minimize(self):
		"""
		If minimize is set to True, then the problem is a "minimization" problem;
		this means that the SDP should be set up so that the largest possible lower bound
		on the density is found.
		
		The minimize property is always the negation of the maximize property.
		
		By default, minimize is False.
		
		"""
		return self._minimize

	@minimize.setter
	def minimize(self, value):
		if not type(value) == bool:
			raise ValueError
		self._minimize = value

	@property
	def maximize(self):
		"""
		If maximize is set to True, then the problem is a "maximization" problem;
		this means that the SDP should be set up so that the smallest possible upper bound
		on the density is found.
		
		The maximize property is always the negation of the minimize property.
		
		By default, maximize is True.
		
		"""

		return not self._minimize

	@maximize.setter
	def maximize(self, value):

		self._register_progression("set_objective", "set")
	
		if not type(value) == bool:
			raise ValueError
		self._minimize = not value


	def _compute_densities(self):
	
		self._densities = [[sum(g.subgraph_density(dg) for dg in self._density_graphs)
			for g in self._graphs]]


	def set_density(self, *args):

		self._register_progression("set_objective", "set")

		density_graphs = []
		orders = []

		for h in args:
		
			if type(h) == str and "." in h:
				h = tuple(map(int, h.split(".")))
		
			if type(h) == tuple:
				k, ne = h
				if k < self._flag_cls().r:
					raise ValueError
				max_e = self._flag_cls.max_number_edges(k)
				if not ne in range(max_e + 1):
					raise ValueError
		
				# Don't forbid anything - if we do, we'll have to keep list
				# updated whenever forbidden things change. Also it seems the most
				# appropriate thing to do.
				graphs = self._flag_cls.generate_graphs(k)
				for g in graphs:
					if g.ne == ne:
						density_graphs.append(g)
				orders.append(k)
				continue
			
			if type(h) == str:
				h = self._flag_cls(h)
				
			if type(h) != self._flag_cls:
				raise ValueError

			density_graphs.append(copy(h))
			orders.append(h.n)
		
		if len(density_graphs) == 0:
			raise ValueError
		
		if len(set(orders)) > 1:
			raise ValueError("Density graphs must all contain the same number of vertices.")
		
		self._density_graphs = density_graphs
		self._compute_densities()


	def _forbid_graph(self, h, induced):

		self._register_progression("specify", "set")

		if type(h) == str and "." in h:
			h = tuple(map(int, h.split(".")))
		
		if type(h) == tuple:
			k, ne = h
			if k < self._flag_cls().r:
				raise ValueError
			max_e = self._flag_cls.max_number_edges(k)
			if not ne in range(max_e + 1):
				raise ValueError
			if induced:
				self._forbidden_edge_numbers.append((k, ne))
			else:
				for i in range(ne, max_e + 1):
					self._forbidden_edge_numbers.append((k, i))
			return
		
		if type(h) == str:
			h = self._flag_cls(h)
		
		if type(h) != self._flag_cls:
			raise ValueError

		if induced:
			self._forbidden_induced_graphs.append(copy(h))
			self._forbidden_induced_graphs.sort(key = lambda g : (g.n, g.ne))
		else:
			self._forbidden_graphs.append(copy(h))
			self._forbidden_graphs.sort(key = lambda g : (g.n, g.ne))


	# TODO: allow lists

	def forbid_subgraph(self, *args):
		for h in args:
			self._forbid_graph(h, False)


	def forbid_induced_subgraph(self, *args):
		for h in args:
			self._forbid_graph(h, True)
	

	def forbid_homomorphic_images(self):
	
		L = sum([g.homomorphic_images() for g in self._forbidden_graphs], [])
		LM = self._flag_cls.minimal_by_inclusion(L)
		if len(LM) == 0:
			return
		sys.stdout.write("Forbidding")
		for g in LM:
			sys.stdout.write(" %s" % repr(g))
			self._forbid_graph(g, False)
		sys.stdout.write("\n")


	# TODO: warn if already solved

	def set_inactive_types(self, *args):
	
		self._register_progression("set_active_types", "set")
	
		num_types = len(self._types)
		
		for arg in args:
			ti = int(arg)
			if not ti in range(num_types):
				raise ValueError
			if ti in self._active_types:
				self._active_types.remove(ti)
			else:
				sys.stdout.write("Warning: type %d is already inactive.\n" % ti)


	def set_extremal_construction(self, c):

		self._register_progression("set_construction", "set")

		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)

		self._construction = c
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
			self._sharp_graph_densities = []
			
			for pair in sharp_graphs:
				g, den = pair
				for gi in range(num_graphs):
					if g.is_labelled_isomorphic(self._graphs[gi]):
						self._sharp_graphs.append(gi)
						self._sharp_graph_densities.append(den)
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

	def add_zero_eigenvectors(self, ti, M, use_bases=False):

		self._register_progression("add_zero_eigenvectors", "set")

		if use_bases:
			self._register_progression("compute_flag_bases", "ensure-not-no")
			NZ = (self._flag_bases[ti].T).solve_left(M)
		else:
			NZ = M

		self._zero_eigenvectors[ti] = self._zero_eigenvectors[ti].stack(NZ)
		self._zero_eigenvectors[ti].set_immutable()

	
	# TODO: is good idea to assume zero densities?
	
	def add_sharp_graphs(self, *args):
	
		self._register_progression("add_sharp_graphs", "set")
	
		num_graphs = len(self._graphs)
	
		for arg in args:
			si = int(arg)
			if not si in range(num_graphs):
				raise ValueError
			if not si in self._sharp_graphs:
				self._sharp_graphs.append(si)
				self._sharp_graph_densities.append(Integer(0))
			else:
				sys.stdout.write("Warning: graph %d is already marked as sharp.\n" % si)


	# TODO: handle non-rational flag bases?

	def change_problem_bases(self, use_blocks=True):

		if self._register_progression("compute_flag_bases", "query") != "set":
			self.compute_flag_bases(use_blocks)

		self._register_progression("transform_problem", "set")

		num_graphs = len(self._graphs)
		num_types = len(self._types)

		for ti in range(num_types):

			sys.stdout.write("Transforming type %d products...\n" % ti)

			nf = len(self._flags[ti])
			B = self._flag_bases[ti]

			rarray = self._product_densities_arrays[ti]
			new_rarray = numpy.zeros([0, 5], dtype=numpy.int)
			row_num = 0

			for gi in range(num_graphs):
			
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
			
			self._product_densities_arrays[ti] = new_rarray

	
	def change_solution_bases(self, use_blocks=True):

		if self._register_progression("compute_flag_bases", "query") != "set":
			self.compute_flag_bases(use_blocks)

		self._register_progression("transform_solution", "set")

		num_types = len(self._types)
	
		sys.stdout.write("Transforming matrices")
		
		self._sdp_Qdash_matrices = []

		for ti in range(num_types):
			
			B = self._flag_bases[ti]
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


	def compute_block_bases(self):

		self._register_progression("compute_block_bases", "set")

		self._block_bases = []
		
		for ti in range(len(self._types)):
		
			B = self._flag_cls.flag_basis(self._types[ti], self._flags[ti])
			num_blocks, block_sizes, block_offsets = block_structure(B)
			sys.stdout.write("Type %d (%d flags) blocks: %s \n" % (ti, len(self._flags[ti]), block_sizes))
			self._block_bases.append(B)


	def compute_flag_bases(self, use_blocks=True, keep_rows=False, use_smaller=False):

		self._register_progression("compute_flag_bases", "set")

		num_types = len(self._types)

		if use_blocks and self._register_progression("compute_block_bases", "query") != "set":
			self.compute_block_bases()

		self._flag_bases = []
		
		sys.stdout.write("Creating bases")
		sys.stdout.flush()

		for ti in range(num_types):
	
			if use_blocks:
				num_blocks, block_sizes, block_offsets = block_structure(self._block_bases[ti])
			else:
				num_blocks, block_sizes, block_offsets = 1, [len(self._flags[ti])], [0]
				
			BS = []
			
			for bi in range(num_blocks):	
			
				Z = self._zero_eigenvectors[ti]
				
				if use_blocks:
					B = (self._block_bases[ti].subdivision(bi, 0) * Z.T).T
				else:
					B = Z
				
				B = B.echelon_form()

				nzev = B.rank()
				B = B[:nzev, :]

				if nzev == 0:
					B = identity_matrix(QQ, block_sizes[bi], sparse=True)
				
				elif nzev == block_sizes[bi]:
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
			
			if M.nrows() == 0:
				M = matrix(self._field, 0, len(self._flags[ti]), sparse=True)
			else:
				M = safe_gram_schmidt(M)
			
			M.set_immutable()
			self._flag_bases.append(M)
	
			sys.stdout.write(".")
			sys.stdout.flush()

		sys.stdout.write("\n")
	
		self._inverse_flag_bases = []

		for ti in range(num_types):
		
			M = copy(self._flag_bases[ti])
			
			if use_smaller:
				MT = M.T
				MT.set_immutable()
				self._inverse_flag_bases.append(MT)
				for j in range(M.nrows()):
					M[j, :] /= sum([x**2 for x in M.row(j)])
				M.set_immutable()	
				self._flag_bases[ti] = M
			
			else:
				for j in range(M.nrows()):
					M[j, :] /= sum([x**2 for x in M.row(j)])
				MT = M.T
				MT.set_immutable()	
				self._inverse_flag_bases.append(MT)


	def compute_products(self):
 	
	 	self._register_progression("compute_products", "set")
 	
	 	num_types = len(self._types)
 		graph_block = make_graph_block(self._graphs, self.n)
		self._product_densities_arrays = []
 		
		sys.stdout.write("Computing products")

		for ti in range(num_types):

			tg = self._types[ti]
			s = tg.n
			m = (self._n + s) / 2
			
			flags_block = make_graph_block(self._flags[ti], m)
			rarray = self._flag_cls.flag_products(graph_block, tg, flags_block, None)
			self._product_densities_arrays.append(rarray)
	
			sys.stdout.write(".")
			sys.stdout.flush()
			
		sys.stdout.write("\n")


	def _set_block_matrix_structure(self):
	
		self._block_matrix_structure = []

		for ti in self._active_types:

			if self._register_progression("transform_problem", "query") == "set":
				num_blocks, block_sizes, block_offsets = block_structure(self._flag_bases[ti])
			else:
				num_blocks, block_sizes, block_offsets = 1, [len(self._flags[ti])], [0]
			
			# Remove zero-sized blocks
			bi = 0
			while bi < num_blocks:
				if block_sizes[bi] == 0:
					num_blocks -= 1
					del block_sizes[bi]
					del block_offsets[bi]
				else:
					bi += 1
			
			for bi in range(num_blocks):
				self._block_matrix_structure.append((ti, block_sizes[bi], block_offsets[bi]))


	def _get_block_matrix_structure(self, ti):
		
		num_blocks = 0
		block_indices = []
		block_offsets = []
		block_sizes = []

		for bi in range(len(self._block_matrix_structure)):
			b = self._block_matrix_structure[bi]
			if b[0] == ti:
				num_blocks += 1
				block_indices.append(bi)
				block_sizes.append(b[1])
				block_offsets.append(b[2])
		
		return (num_blocks, block_sizes, block_offsets, block_indices)


	# TODO: helpful error message if product densities have not been computed.
	# TODO: add option for forcing sharps
	# TODO: should block matrix structure be set elsewhere?
	
	def write_sdp_input_file(self, no_output=False):
	
		self._register_progression("write_sdp_input_file", "set")
		
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)
		
		# For maximization problems, the objective value returned by the SDP solver
		# must be negated.
		self._obj_value_factor = 1.0 if self._minimize else -1.0
		
		if num_densities < 1:
			raise NotImplementedError("there must be at least one density.")
		
		self._sdp_input_filename = os.path.join(SAGE_TMP, "sdp.dat-s")
	
		self._set_block_matrix_structure()
		total_num_blocks = len(self._block_matrix_structure)
		
		if no_output:
			return
		
		sys.stdout.write("Writing SDP input file...\n")
				
		with open(self._sdp_input_filename, "w") as f:
	
			f.write("%d\n" % (num_graphs + 1,))
			f.write("%d\n" % (total_num_blocks + 3,))
			
			f.write("1 ")
			for b in self._block_matrix_structure:
				f.write("%d " % b[1])
			
			f.write("-%d -%d\n" % (num_graphs, num_densities))
			f.write("0.0 " * num_graphs)
			f.write("1.0\n")
			
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
					f.write("%d %d %d %d 1.0\n" % (i + 1, total_num_blocks + 2, i + 1, i + 1))
	
			for i in range(num_graphs):
				for j in range(num_densities):
					d = self._densities[j][i]
					if d != 0:
						if self._minimize:
							d *= -1
						f.write("%d %d %d %d %s\n" % (i + 1, total_num_blocks + 3, j + 1, j + 1, d.n(digits=64)))

			non_free_densities = range(num_densities)
			if hasattr(self, "_free_densities"):
				for j in self._free_densities:
					non_free_densities.remove(j)
			
			for j in non_free_densities:
				f.write("%d %d %d %d 1.0\n" % (num_graphs + 1, total_num_blocks + 3, j + 1, j + 1))
		
			for ti in self._active_types:

				num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)
		
				for row in self._product_densities_arrays[ti]:
					gi = row[0]
					j = row[1]
					k = row[2]
					bi = num_blocks - 1
					if bi > 0:
						while block_offsets[bi] > j:
							bi -= 1
						j -= block_offsets[bi]
						k -= block_offsets[bi]
					value = Integer(row[3]) / Integer(row[4])
					f.write("%d %d %d %d %s\n" %
						(gi + 1, block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))


	def write_sdp_initial_point_file(self, small_change=1/Integer(10)):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)
		total_num_blocks = len(self._block_matrix_structure)
		
		self._sdp_initial_point_filename = os.path.join(SAGE_TMP, "sdp.ini-s")
	
		self._set_block_matrix_structure()
		num_blocks = len(self._block_matrix_structure)
		
		sys.stdout.write("Writing SDP initial point file...\n")
				
		with open(self._sdp_initial_point_filename, "w") as f:
	
			for gi in range(num_graphs):
				if gi in self._sharp_graphs:
					si = self._sharp_graphs.index(gi)
					f.write("%s " % self._sharp_graph_densities[si].n(digits=64))
				else:
					f.write("0.0 ")
			
			if not self._minimize:
				f.write("%s\n" % (-self._target_bound).n(digits=64))
			else:
				f.write("%s\n" % self._target_bound.n(digits=64))
			
			f.write("1 1 1 1 %s\n" % small_change.n(digits=64))
			
			for ti in range(num_types):
				
				nf = len(self._flags[ti])
				z_matrix = matrix(QQ, nf, nf)
				
				for row in self._product_densities_arrays[ti]:
					gi = row[0]
					if not gi in self._sharp_graphs:
						continue
					si = self._sharp_graphs.index(gi)
					j = row[1]
					k = row[2]
					value = Integer(row[3]) / Integer(row[4])
					z_matrix[j, k] += value * self._sharp_graph_densities[si]

				for j in range(nf):
					z_matrix[j, j] += small_change

				num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

				for bi in range(num_blocks):
					for j in range(block_sizes[bi]):
						for k in range(j, block_sizes[bi]):
							value = z_matrix[block_offsets[bi] + j, block_offsets[bi] + k]
							if value != 0:
								f.write("1 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

			for gi in range(num_graphs):
				if gi in self._sharp_graphs:
					si = self._sharp_graphs.index(gi)
					value = self._sharp_graph_densities[si]
				else:
					value = 0
				if value <= 0:
					value = small_change
				f.write("1 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

			for j in range(num_densities):
				f.write("1 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, small_change.n(digits=64)))

			if hasattr(self, "_exact_Q_matrices"):

				f.write("2 1 1 1 %s\n" % self._bound.n(digits=64))

				for ti in range(num_types):

					num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)
	
					for bi in range(num_blocks):
						for j in range(block_sizes[bi]):
							for k in range(j, block_sizes[bi]):
								value = self._exact_Q_matrices[ti][block_offsets[bi] + j, block_offsets[bi] + k]
								if j == k:
									value += small_change
								if value != 0:
									f.write("2 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

				for gi in range(num_graphs):
					value = self._bound - self._bounds[gi]
					if value <= 0:
						value = small_change
					f.write("2 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

				for j in range(num_densities):
					value = self._exact_density_coeffs[j]
					if value <= 0:
						value = small_change
					f.write("2 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, value.n(digits=64)))

			else:

				densities = [sum([self._densities[j][gi] for j in range(num_densities)])
						/ num_densities for gi in range(num_graphs)]
				
				bound = min(densities) if self._minimize else max(densities)
				
				value = bound
				if value <= 0:
					value = small_change
				f.write("2 1 1 1 %s\n" % value.n(digits=64))
				
				num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)
	
				for ti in range(num_types):
				
					num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)
					
					for bi in range(num_blocks):
						for j in range(block_sizes[bi]):
							f.write("2 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, j + 1, small_change.n(digits=64)))
				
				for gi in range(num_graphs):
					value = (bound - densities[gi]) * (-1 if self._minimize else 1)
					if value <= 0:
						value = small_change
					f.write("2 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))
				
				value = Integer(1) / num_densities
				for j in range(num_densities):
					f.write("2 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, value.n(digits=64)))


	# TODO: report error if problem infeasible
	
	def run_sdp_solver(self, show_output=False, sdpa=False, output_file=None, use_initial_point=True):
	
		self._register_progression("run_sdp_solver", "set")
	
		if not output_file is None:
			self._sdp_output_filename = output_file
			self._read_sdp_output_file()
			return
	
		self._sdp_output_filename = os.path.join(SAGE_TMP, "sdp.out")

		if not sdpa:

			cmd = "%s %s %s" % (cdsp_cmd, self._sdp_input_filename, self._sdp_output_filename)

			if use_initial_point and hasattr(self, "_sdp_initial_point_filename"):
				cmd += " %s" % self._sdp_initial_point_filename

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

		self._read_sdp_output_file()


	def _read_sdp_output_file(self):

		num_graphs = len(self._graphs)
		num_types = len(self._types)
		num_densities = len(self._densities)

		with open(self._sdp_output_filename, "r") as f:
		
			if self._register_progression("transform_problem", "query") == "set":

				self._sdp_Q_matrices = [matrix(self._approximate_field,
					self._flag_bases[ti].nrows(), self._flag_bases[ti].nrows())
					for ti in range(num_types)]
				
				for ti in range(num_types):
					B = self._flag_bases[ti]
					row_div = B.subdivisions()[0]
					self._sdp_Q_matrices[ti].subdivide(row_div, row_div)

			else:
	
				self._sdp_Q_matrices = [matrix(self._approximate_field,
					len(self._flags[ti]), len(self._flags[ti])) for ti in range(num_types)]
			
			self._sdp_density_coeffs = [self._approximate_field(0) for i in range(num_densities)]

			num_blocks = len(self._block_matrix_structure)
			
			for line in f:
				numbers = line.split()
				if numbers[0] != "2":
					continue
				bi = int(numbers[1]) - 2
				if bi == num_blocks + 1:
					j = int(numbers[2]) - 1
					self._sdp_density_coeffs[j] = self._approximate_field(numbers[4])
					continue
				if bi < 0 or bi >= num_blocks:
					continue
				j = int(numbers[2]) - 1
				k = int(numbers[3]) - 1
				ti, size, offset = self._block_matrix_structure[bi]
				j += offset
				k += offset 

				self._sdp_Q_matrices[ti][j, k] = numbers[4]
				self._sdp_Q_matrices[ti][k, j] = self._sdp_Q_matrices[ti][j, k]

		for ti in range(num_types):
			self._sdp_Q_matrices[ti].set_immutable()


	def check_floating_point_bound(self, tolerance = 0.00001, show_sorted=False, show_all=False):
	
		self._register_progression("run_sdp_solver", "ensure")
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)
		
		sys.stdout.write("Checking numerical bound...\n")
		
		fbounds = [sum([self._densities[j][i] * self._sdp_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]

		for ti in self._active_types:
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
		
		self._sdp_bounds = fbounds

		if self._register_progression("set_construction", "query") == "set":

			if abs(bound - self._approximate_field(self._target_bound)) < tolerance:
				sys.stdout.write("Bound of %s appears to have been met.\n" % self._target_bound)
			else:
				sys.stdout.write("Warning: bound of %s appears to have not been met.\n" % self._target_bound)
				return
			sharp_graphs = self._sharp_graphs
		else:
			sharp_graphs = [] # set dummy sharp_graphs
				
		apparently_sharp_graphs = [gi for gi in range(num_graphs) if abs(fbounds[gi] - bound) < tolerance]

		if show_sorted or show_all:
		
			if not self._minimize:
				sorted_indices = sorted(range(num_graphs), key = lambda i : -fbounds[i])
			else:
				sorted_indices = sorted(range(num_graphs), key = lambda i : fbounds[i])

			for gi in sorted_indices:
				if gi in apparently_sharp_graphs:
					sys.stdout.write("S")
				elif not show_all:
					break
				else:
					sys.stdout.write(" ")
				if gi in self._sharp_graphs:
					sys.stdout.write("C")
				else:
					sys.stdout.write(" ")
				sys.stdout.write(" %s : graph %d (%s) " % (fbounds[gi], gi, self._graphs[gi]))
				sys.stdout.write("\n")

		else:
	
			sys.stdout.write("The following %d graphs appear to be sharp:\n" % len(apparently_sharp_graphs))
			for gi in apparently_sharp_graphs:
				sys.stdout.write("%.12f : graph %d (%s)\n" % (fbounds[gi], gi, self._graphs[gi]))
		
		if self._register_progression("set_construction", "query") != "set":
			return
			
		extra_sharp_graphs = [gi for gi in apparently_sharp_graphs if not gi in self._sharp_graphs]
		missing_sharp_graphs = [gi for gi in self._sharp_graphs if not gi in apparently_sharp_graphs]
			
		if len(extra_sharp_graphs) > 0:
			sys.stdout.write("Warning: additional sharp graphs: %s\n" % (extra_sharp_graphs,))	
	
		for gi in missing_sharp_graphs:
			sys.stdout.write("Warning: graph %d (%s) does not appear to be sharp.\n" % (gi, self._graphs[gi]))



	def solve_sdp(self, show_output=False, sdpa=False, tolerance=0.00001, show_all=False, output_file=None):
		"""
		Creates and solves the semi-definite program corresponding to the problem.
		
		Equivalent to calling write_sdp_input_file, followed by run_sdp_solver,
		followed by check_floating_point_bound.
		
		"""
		
		if output_file is None:
			self.write_sdp_input_file()
			self.write_sdp_initial_point_file()
			self.run_sdp_solver(show_output=show_output, sdpa=sdpa)
		else:
			self.write_sdp_input_file(no_output=True)
			self.run_sdp_solver(show_output=show_output, sdpa=sdpa, output_file=output_file)

		self.check_floating_point_bound(tolerance=tolerance, show_all=show_all)


	def import_solution(self, directory, complement=False):
	
		self._register_progression("write_sdp_input_file", "set")
		self._register_progression("run_sdp_solver", "set")
	
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

		if flags.n != self.n:
			raise ValueError

		if self._flag_cls().oriented:
			if not ("oriented %d-graph" % self._flag_cls().r) in flags.description:
				raise ValueError
		else:
			if not ("%d-graph" % self._flag_cls().r) in flags.description:
				raise ValueError

		if flags.num_H != num_graphs:
			raise ValueError

		if flags.num_types != num_types:
			raise ValueError

		ttr = []
		ftrs = []
	
		for ti in range(num_types):
			tg = self._flag_cls(flags.types[ti])
			if complement:
				tg = tg.complement(True)
			for tj in range(num_types):
				if tg.is_labelled_isomorphic(self._types[tj]):
					ttr.append(tj)
					break
			else:
				raise ValueError
		
			num_flags = len(self._flags[tj])
			ftr = []
			
			for fi in range(num_flags):
				fg = self._flag_cls(flags.flags[ti][fi])
				fg.t = tg.n
				if complement:
					fg = fg.complement(True)
				for fj in range(num_flags):
					if fg.is_labelled_isomorphic(self._flags[tj][fj]):
						ftr.append(fj)
						break
				else:
					#print ti
					#print ftr
					raise ValueError
					
			ftrs.append(ftr)
		
		#print ttr
		#print ftrs
		
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
				if tj in self._active_types:
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


	def show_zero_eigenvalues(self, tolerance=0.00001, types=None):
	
		self._register_progression("run_sdp_solver", "ensure")

		if types is None:
			types = self._active_types
	
		for ti in types:
			eigvals = numpy.linalg.eigvalsh(self._sdp_Q_matrices[ti])
			zero_eigvals = sorted([e for e in eigvals if e < tolerance])
			if len(zero_eigvals) == 0:
				sys.stdout.write("Type %d. None.\n" % ti)
			else:
				sys.stdout.write("Type %d. %d possible: %s.\n" % (ti, len(zero_eigvals),
					" ".join("%s" % e for e in zero_eigvals)))


	# TODO: transform with _flag_bases if present.

	def check_construction(self, C, tolerance = 0.00001, types=None):
	
		self._register_progression("run_sdp_solver", "ensure")
		
		if types is None:
			types = self._active_types
		
		for ti in types:
			M = C.zero_eigenvectors(self._types[ti], self._flags[ti]) 
			if M.nrows() == 0:
				sys.stdout.write("Type %d. None.\n" % ti)
			else:
				R = M * self._sdp_Q_matrices[ti]
				norms = [R[i,:].norm() for i in range(R.nrows())]
				sys.stdout.write("Type %d. %d possible: %s.\n" % (ti, len(norms),
					" ".join("%s" % e for e in norms)))


	def find_extra_zero_eigenvectors(self, ti, tolerance=1e-5, threshold=1e-10, echelonize=True, denominator=None):
	
		self._register_progression("run_sdp_solver", "ensure")
	
		if not ti in self._active_types:
			raise ValueError("Type is not active.")
	
		M = self._zero_eigenvectors[ti].echelon_form()
		E = copy(self.get_zero_eigenvectors(ti, tolerance=tolerance, use_bases=False))
		
		pivots = M.pivots()
		for i in range(len(pivots)):
			c = pivots[i]
			for r in range(E.nrows()):
				E[r, :] -= E[r, c] * M[i, :]
	
		if not echelonize:
			return E
	
		for r in range(E.nrows()):
			for c in range(E.ncols()):
				if abs(E[r, c]) > threshold:
					E[r, :] /= E[r, c]
					for s in range(E.nrows()):
						if s != r:
							E[s, :] -= E[s, c] * E[r, :]
					break
	
		for r in range(E.nrows()):
			for c in range(E.ncols()):
				if abs(E[r, c]) < threshold:
					 E[r, c] = 0
	
		r = E.nrows() - 1
		while r >= 0 and E[r, :].is_zero():
			E = E[:r, :]
			r -= 1
	
		if denominator is None:
			return E
		
		ER = matrix(QQ, E.nrows(), E.ncols())
	
		def rationalize (f):
			return Integer(round(f * denominator)) / denominator
	
		for r in range(E.nrows()):
			for c in range(E.ncols()):
				ER[r, c] = rationalize(E[r, c])
	
		return ER


	def get_zero_eigenvectors(self, ti, tolerance=1e-5, use_bases=True):

		self._register_progression("run_sdp_solver", "ensure")
	
		if not ti in self._active_types:
			raise ValueError("Type is not active.")

		if use_bases and self._register_progression("transform_solution", "query"):
			QM = self._sdp_Qdash_matrices[ti]
		else:
			QM = self._sdp_Q_matrices[ti]

		ns = len(QM.subdivisions()[0]) + 1
	
		B = []
		for i in range(ns):
			QB = QM.subdivision(i, i)
			eigvals, T = numpy.linalg.eigh(QB)
			M = matrix(self._approximate_field, 0, QB.ncols())
			for ei in range(len(eigvals)):
				if eigvals[ei] < tolerance:
					M = M.stack(matrix(numpy.matrix(T[:, ei])))
			B.append(M)
		return block_diagonal_matrix(B)



	def make_exact(self, denominator=1024, cholesky=[], protect=[], meet_target_bound=True, show_changes=False,
		use_densities=True):
	
		self._register_progression("make_exact", "set")
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)

		if cholesky == "all":
			cholesky = self._active_types

		if meet_target_bound:
			self._register_progression("set_construction", "ensure")
			num_sharps = len(self._sharp_graphs)
		else:
			if not all(ti in cholesky for ti in self._active_types):
				raise NotImplementedError("If construction is not set, then cholesky must be used.")
			num_sharps = 0
			
		if self._register_progression("transform_solution", "query") != "set":
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
		
		if not meet_target_bound:
			for ti in range(num_types):
				self._exact_Qdash_matrices[ti].set_immutable()
			return
		
		
		triples = [(ti, j, k) for ti in self._active_types for j in range(q_sizes[ti])
			for k in range(j, q_sizes[ti])]
	
		num_triples = len(triples)
		triples.sort()
		triple_to_index = dict((triples[i], i) for i in range(num_triples))

		R = matrix(self._field, num_sharps, num_triples, sparse=True)

		sys.stdout.write("Constructing R matrix")

		# TODO: only use triples that correspond to middle blocks.

		for ti in self._active_types:

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

			if self._register_progression("transform_solution", "query") == "set":
				B = self._inverse_flag_bases[ti]
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
		DR = matrix(self._field, num_sharps, 0) # sparsity harms performance too much here
		EDR = DR.T
		
		sys.stdout.write("Constructing DR matrix")
		
		# Only if there is more than one density
		if num_densities > 1 and use_densities:
			
			for j in range(num_densities):
				new_col = matrix(QQ, [[self._densities[j][gi]] for gi in self._sharp_graphs])
				EDR = EDR.stack(new_col.T)
				EDR.echelonize()
				if EDR[-1, :].is_zero():
					EDR = EDR[:-1, :]
					sys.stdout.write("~")
					sys.stdout.flush()
					continue		
				
				DR = DR.augment(new_col)
				density_cols_to_use.append(j)
				sys.stdout.write(".")
				sys.stdout.flush()
			
			sys.stdout.write("\n")
			sys.stdout.write("DR matrix (density part) has rank %d.\n" % DR.ncols())
				
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
			new_col = R[:, i : i + 1]
			EDR = EDR.stack(new_col.T)
			EDR.echelonize()
			if EDR[-1, :].is_zero():
				EDR = EDR[:-1, :]
				sys.stdout.write("~")
				sys.stdout.flush()
				continue		
			
			DR = DR.augment(new_col)
			cols_to_use.append(i)
			sys.stdout.write(".")
			sys.stdout.flush()
		
		sys.stdout.write("\n")
		sys.stdout.write("DR matrix has rank %d.\n" % DR.ncols())
		
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

		FDR = matrix(self._field, DR)
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
	
		self._register_progression("make_exact", "ensure")
	
		for ti in self._active_types:
			sys.stdout.write("Type %d:\n" % ti)
			original_eigvals = sorted(numpy.linalg.eigvalsh(self._sdp_Qdash_matrices[ti]))
			new_eigvals = sorted(numpy.linalg.eigvalsh(self._exact_Qdash_matrices[ti]))
			if only_smallest:
				num = 1
			else:
				num = len(original_eigvals)
			for i in range(num):
				sys.stdout.write("%.11f : %.11f\n" % (original_eigvals[i], new_eigvals[i]))


	# TODO: check numerically for negative (and zero) eigenvalues
	# TODO: check for negative density coefficients 

	def check_exact_bound(self):
	
		self._register_progression("check_exact", "set")
	
		num_types = len(self._types)
		num_graphs = len(self._graphs)
		num_densities = len(self._densities)
		
		self._exact_Q_matrices = []
		
		if self._register_progression("transform_solution", "query") == "set":
			for ti in range(num_types):
				B = self._inverse_flag_bases[ti]
				M = B * self._exact_Qdash_matrices[ti] * B.T
				self._exact_Q_matrices.append(M)
		else:
			self._exact_Q_matrices = self._exact_Qdash_matrices
		
		bounds = [sum([self._densities[j][i] * self._exact_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]

		for ti in self._active_types:
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
			else:
				bound = min(bounds)
		else:
			# Sorting doesn't currently work for number fields with embeddings, so use float approximation.
			if not self._minimize:
				bound = max(bounds, key = lambda x : float(x))
			else:
				bound = min(bounds, key = lambda x : float(x))

		self._bounds = bounds
		self._bound = bound

		if self._register_progression("set_construction", "query") != "set":
			return

		if self._field == RationalField():
			if not self._minimize:
				violators = [gi for gi in range(num_graphs) if bounds[gi] > self._target_bound]
			else:
				violators = [gi for gi in range(num_graphs) if bounds[gi] < self._target_bound]
		else:
			if not self._minimize:
				violators = [gi for gi in range(num_graphs) if float(bounds[gi]) > float(self._target_bound)]
			else:
				violators = [gi for gi in range(num_graphs) if float(bounds[gi]) < float(self._target_bound)]

		sys.stdout.write("Bound of %s attained by:\n" % self._target_bound)
		for gi in range(num_graphs):
			if bounds[gi] == self._target_bound:
				sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))

		if len(violators) > 0:
			sys.stdout.write("Bound violated by:")
			for gi in violators:
				sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))


	# TODO: use self._register_progression

	def diagonalize(self):
		
		def LDLdecomposition(M):	# TODO: does this handle matrices with zero eigenvalues?
	
			MS = M.parent()
			D = MS.matrix()
			if M.is_zero():
				D.set_immutable()
				return D, D
			L = copy(MS.identity_matrix())
			for i in xrange(M.nrows()):
				for j in xrange(i):
					L[i, j] = (Integer(1) / D[j, j]) * (M[i, j] - sum(L[i, k] * L[j, k] * D[k, k] for k in xrange(j)))
				D[i, i] = M[i, i] - sum(L[i, k]**2 * D[k, k]
					for k in xrange(i))
			L.set_immutable()
			D.set_immutable()
			return L, D

		self._exact_diagonal_matrices = []
		self._exact_r_matrices = []

		sys.stdout.write("Diagonalizing")

		# TODO: what to do if problem has been transformed?

		for ti in range(len(self._types)):
			R, M = LDLdecomposition(self._exact_Qdash_matrices[ti])
			self._exact_diagonal_matrices.append(M)
			if self._register_progression("transform_solution", "query") == "set":
				R = self._inverse_flag_bases[ti] * R
			self._exact_r_matrices.append(R)
			sys.stdout.write(".")
			sys.stdout.flush()
		
		# Q can now be computed as Q = R * D * R.T
		
		sys.stdout.write("\nVerifying")

		for ti in range(len(self._types)):
			Q = self._exact_r_matrices[ti] * self._exact_diagonal_matrices[ti] * self._exact_r_matrices[ti].T
			if Q != self._exact_Q_matrices[ti]:
				raise ValueError # TODO: choose appropriate error
			sys.stdout.write(".")
			sys.stdout.flush()
		
		sys.stdout.write("\n")


# 	def write_certificate(self, filename):
# 
# 
# 		data = {
# 			"description" : flags.description,
# 			"bound" : self._bound,
# 			"order_of_admissible_graphs" : self._n,
# 			"number_of_admissible_graphs" : len(self._graphs),
# 			"admissible_graphs" : self._graphs,
# 			"admissible_graph_densities" : array_to_json(self._densities[0]), # TODO: multiple densities
# 			"number_of_types" : len(self._types),
# 			"types" : self._types,
# 			"numbers_of_flags" : [len(L) for L in self._flags],
# 			"flags" : self._flags,
# 			"qdash_matrices" : array_to_json(self._exact_diagonal_matrices),
# 			"r_matrices" : array_to_json(self._exact_r_matrices)
# 		}


def ThreeGraphProblem():
	return Problem(ThreeGraphFlag)

def GraphProblem():
	return Problem(GraphFlag)

def OrientedGraphProblem():
	return Problem(OrientedGraphFlag)

def TwoMultigraphProblem():
	return Problem(TwoMultigraphFlag)

def ThreeMultigraphProblem():
	return Problem(ThreeMultigraphFlag)

