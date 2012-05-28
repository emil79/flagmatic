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

from sage.rings.all import Integer

from hypergraph_flag import make_graph_block
from three_graph_flag import *
from graph_flag import *
from oriented_graph_flag import *
from problem import *


class AxiomsProblem(Problem):

	def __init__(self, flag_cls):
	
		Problem.__init__(self, flag_cls)
		self._axioms = []
		self._axiom_flags = []
		

	def clear_densities(self):
		self._densities = []


	def add_axiom(self, tg, terms, make_free=True):

		self.state("set_objective", "yes")

		m = self.n - max([t[0].n for t in terms]) + tg.n

		axiom_flags = self._flag_cls.generate_flags(m, tg,
			forbidden_edge_numbers=self._forbidden_edge_numbers,
			forbidden_graphs=self._forbidden_graphs,
			forbidden_induced_graphs=self._forbidden_induced_graphs)
		
		num_densities = len(axiom_flags)
		sys.stdout.write("Added %d quantum graphs.\n" % num_densities)
		
		num_graphs = len(self._graphs)
		quantum_graphs = [[Integer(0) for i in range(num_graphs)] for j in range(num_densities)]
		
		axiom_flags_block = make_graph_block(axiom_flags, m)
		graph_block = make_graph_block(self._graphs, self.n)

		for i in range(len(terms)):
			fg = terms[i][0]
			flags_block = make_graph_block([fg], fg.n)
			rarray = self._flag_cls.flag_products(graph_block, tg, flags_block, axiom_flags_block)
			for row in rarray:
				gi = row[0]
				j = row[1] # always 0
				k = row[2]
				value = Integer(row[3]) / Integer(row[4])
				quantum_graphs[k][gi] += value * terms[i][1]
		
		self._axioms.append((tg, terms))
		self._axiom_flags.append(axiom_flags)
		
		num_densities = len(self._densities)
		self._densities.extend(quantum_graphs)
		if make_free:
			if not hasattr(self, "_free_densities"):
				self._free_densities = []
			self._free_densities.extend(range(num_densities, len(self._densities)))
	
	
	def add_codegree_axiom(self, value, make_free=True):

		if not self._flag_cls().r == 3:
			raise NotImplementedError
	
		tg = ThreeGraphFlag("2:")
		f1 = ThreeGraphFlag("3:123(2)")
		f2 = ThreeGraphFlag("2:(2)")
		self.add_axiom(tg, [(f1, Integer(1)), (f2, -value)], make_free=make_free)


	def add_degree_axiom(self, value, make_free=True):
	
		if self._flag_cls().oriented:
			raise NotImplementedError
	
		if self._flag_cls().r == 3:
	
			tg = ThreeGraphFlag("1:")
			f1 = ThreeGraphFlag("3:123(1)")
			f2 = ThreeGraphFlag("1:(1)")
			self.add_axiom(tg, [(f1, Integer(1)), (f2, -value)], make_free=make_free)

		elif self._flag_cls().r == 2:

			tg = GraphFlag("1:")
			f1 = GraphFlag("2:12(1)")
			f2 = GraphFlag("1:(1)")
			self.add_axiom(tg, [(f1, Integer(1)), (f2, -value)], make_free=make_free)


	# TODO: fix this!

	def add_equal_degrees_axiom(self, make_free=True):
	
		if self._flag_cls().oriented:
			raise NotImplementedError
	
		if self._flag_cls().r == 3:
	
			tg = ThreeGraphFlag("2:")
			f1 = ThreeGraphFlag("4:134(2)")
			f2 = ThreeGraphFlag("4:234(2)")
			self.add_axiom(tg, [(f1, Integer(1)), (f2, -Integer(1))], make_free=make_free)

		elif self._flag_cls().r == 2:

			tg = GraphFlag("2:")
			f1 = GraphFlag("3:13(2)")
			f2 = GraphFlag("3:23(2)")
			self.add_axiom(tg, [(f1, Integer(1)), (f2, -Integer(1))], make_free=make_free)
			tg = GraphFlag("2:12")
			f1 = GraphFlag("3:1213(2)")
			f2 = GraphFlag("3:1223(2)")
			self.add_axiom(tg, [(f1, Integer(1)), (f2, -Integer(1))], make_free=make_free)
		
	
	def add_out_degree_axiom(self, value, make_free=True):
	
		if not (self._flag_cls().r == 2 and self._flag_cls().oriented):
			raise NotImplementedError
	
		tg = OrientedGraphFlag("1:")
		f1 = OrientedGraphFlag("2:12(1)")
		f2 = OrientedGraphFlag("1:(1)")
		self.add_axiom(tg, [(f1, Integer(1)), (f2, -value)], make_free=make_free)


	def add_in_degree_axiom(self, value, make_free=True):
	
		if not (self._flag_cls().r == 2 and self._flag_cls().oriented):
			raise NotImplementedError
	
		tg = OrientedGraphFlag("1:")
		f1 = OrientedGraphFlag("2:21(1)")
		f2 = OrientedGraphFlag("1:(1)")
		self.add_axiom(tg, [(f1, Integer(1)), (f2, -value)], make_free=make_free)	


	def show_large_densities(self, larger_than=0.0):

		self.state("run_sdp_solver", "ensure_yes")

		num_densities = len(self._densities)

		densities_to_use = []
		for j in range(num_densities):
			if self._sdp_density_coeffs[j] > larger_than:
				densities_to_use.append(j)

		sys.stdout.write("Densities: %s\n" % (densities_to_use,))

		sys.stdout.write("Coefficients: %s\n" % ([self._sdp_density_coeffs[j] for j in densities_to_use],))


	def show_independent_densities(self):

		self.state("run_sdp_solver", "ensure_yes")
	
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
	
		if len(self._axioms) != 1 or hasattr(self, "_free_densities"):
			raise NotImplementedError
	
		new_densities = [self._densities[j] for j in densities_to_use]
		new_axiom_flags =  [self._axiom_flags[0][j] for j in densities_to_use]
		
		new_problem = copy(self)
		new_problem._densities = new_densities
		new_problem._axiom_flags = [new_axiom_flags]
		
		if hasattr(new_problem, "_sdp_Q_matrices"):
			del new_problem._sdp_Q_matrices
		if hasattr(new_problem, "_sdp_Qdash_matrices"):
			del new_problem._sdp_Qdash_matrices
		if hasattr(new_problem, "_exact_Q_matrices"):
			del new_problem._exact_Q_matrices
		if hasattr(new_problem, "_exact_Qdash_matrices"):
			del new_problem._exact_Qdash_matrices
		if hasattr(new_problem, "_sdp_density_coeffs"):
			del new_problem._sdp_density_coeffs
		if hasattr(new_problem, "_exact_density_coeffs"):
			del new_problem._exact_density_coeffs
		if hasattr(new_problem, "_sdp_bounds"):
			del new_problem._sdp_bounds
		if hasattr(new_problem, "_bounds"):
			del new_problem._bounds
		
		new_problem.state("set_objective", "yes")

		return new_problem


def ThreeGraphAxiomsProblem():
	return AxiomsProblem(ThreeGraphFlag)

def GraphAxiomsProblem():
	return AxiomsProblem(GraphFlag)

def OrientedGraphAxiomsProblem():
	return AxiomsProblem(OrientedGraphFlag)
