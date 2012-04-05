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

from flagmatic.all import *

def test_graphs():

	gsl = ["6:", "6:123", "6:123124", "6:123124", "6:123124", "6:123124", "6:123124125",
		"6:123124135145", "6:123124125126", "6:123124125136", "6:123124125136146",
		"6:123124125136146156", "6:123124125136146156", "6:123124",
		"6:123124125136146156345", "6:123124135245345", "6:123124135146345346",
		"6:123124125136146356456", "6:123124125136146256356456"]
		
	for gs in gsl:
		print gs
		P = Problem()
		P.forbid_subgraph((4, 3))
		P.forbid_induced_subgraph(Flag(gs))
		P.n = 6
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()

def example(prob):
	
	if prob == "ff83":
	
		P = Problem()
		P.forbid_subgraph("5:123124345")
		P.n = 6
		C = BlowupConstruction(Flag("3:123"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P._approximate_field = RealField(113)
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
	
	elif prob == "ff83n":
	
		P = Problem()
		P.forbid_subgraph("5:123124345")
		P.n = 6
		C = BlowupConstruction(Flag("3:123"))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
	
	
	elif prob == "k4-":

		P = Problem()
		P.forbid_subgraph((4, 3))
		P.n = 7
		#
		C = None
		P.calculate_product_densities()
		#P.write_sdp_input_file()
		#P.run_sdp_solver()


	elif prob == "f32":

		P = Problem()
		P.forbid_subgraph("5:123124125345")
		P.n = 6
		C = UnbalancedBlowupConstruction(Flag("2:122"), weights=[1,2])
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "razb":
	
		P = Problem()
		P.forbid_subgraph((4, 4))
		P.forbid_induced_subgraph((4, 1))
		P.n = 6
		C = BlowupConstruction(Flag("3:112223331123"))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()

	
	elif prob == "k5weak":
	
		P = Problem()
		P.forbid_subgraph((5, 10))
		P.forbid_induced_subgraph((5, 8))
		P.n = 6
		C = BlowupConstruction(Flag("2:112122"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "turank4-":
	
		P = Problem()
		P.forbid_subgraph((4, 4))
		P.n = 6
		P.set_density((4,3))
		C = BlowupConstruction(Flag("3:112223331123"))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "max59":
	
		P = Problem()
		P.n = 6
		P.set_density((5, 9))
		C = BlowupConstruction(Flag("2:112122"))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "max56":
	
		P = Problem()
		P.n = 6
		P.set_density_edge_number(5, 6)
		C = BlowupConstruction(Flag("3:112223331122233311"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "38":

		P = Problem()
		P.forbid_subgraph(Flag("5:123124125345"))
		P.forbid_subgraph(Flag("5:123124125134135145"))
		P.n = 6
		C = BlowupConstruction(Flag("4:123124134234"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "38n":
	
		P = Problem()
		P.forbid_subgraph(Flag("5:123124125345"))
		P.forbid_subgraph(Flag("5:123124125134135145"))
		P.n = 6
		C = BlowupConstruction(Flag("4:123124134234"))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases(use_blocks=False)
		P.make_exact(1024)
		P.check_exact_bound()
		

	elif prob == "max42":

		P = Problem()
		P.n = 5
		P.remove_types([0])
		P.set_density_graph(Flag("4:123124"))
		C = AdHocConstruction("max42")
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
	
	elif prob == "k4-f32":

		P = Problem()
		P.forbid_subgraph((4, 3))
		P.forbid_subgraph(Flag("5:123124125345"))
		P.n = 7
		P.remove_types([3,4,5,6])
		C = BlowupConstruction(Flag("6:123234345451512136246356256146"))
		P.use_construction(C)
		P.calculate_product_densities()
		P.import_solution("../output/k4-f32")
		P.save("k4-f32")
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.save("k4-f32")
		P.make_exact(10000000,protect=[2])
		P.save("k4-f32")
		P.check_exact_bound()
		

	elif prob == "k4-cod-inexact":
	
		P = AxiomsProblem()
		C = None
		P.forbid_subgraph((4, 3))
		P.n = 6
		tg = Flag("2:")
		f1 = Flag("3:123", tg)
		f2 = Flag("2:", tg)
		P.add_axiom(tg, [f1, f2], [1, Rational("-264/1000")])
		
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()



	elif prob == "codr":
	
		P1 = AxiomsProblem()
		P1.forbid_subgraph((4, 3))
		P1.forbid_subgraph(Flag("6:612623634645651"))
		P1.n = 6
		P1.clear_axioms()
		P1.add_codegree_axiom(Rational("1/4"))
		P1._force_sharps = True
		C = RandomTournamentConstruction()
		P1.use_construction(C)
		P1.change_problem_bases()
		P1.calculate_product_densities()
		P1.write_sdp_input_file()
		P1.run_sdp_solver(True)
		P2 = P1.lose_small_densities(0.05)
		P2.write_sdp_input_file()
		P2.run_sdp_solver(True)
		P = P2.combine_densities(20)
		P.write_sdp_input_file()
		P.run_sdp_solver(True)
		P.check_floating_point_bound()
		P.make_exact(1024*1024)
		P.check_exact_bound()

	elif prob == "cod6":

		P = AxiomsProblem()
		P.forbid_subgraph((4, 3))
		P.forbid_subgraph(Flag("6:612623634645651"))
		P.n = 6
		P.clear_axioms()
		P.add_codegree_axiom(1/4)
		C = RandomTournamentConstruction()
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact(2^20)
		P.check_exact_bound()


	elif prob == "new-cod-1":
	
		P = AxiomsProblem()
		P.forbid_subgraph((4, 3))
		P.n = 7
		P.remove_types([3,4])
		P.clear_axioms()
		P.add_codegree_axiom(1/4)
		C = RandomTournamentConstruction()
		return P, C

		P.use_construction(C)
		P.calculate_product_densities()
		P.save("new-cod-1")
		P._force_sharps = True
		P.save("new-cod-1")
		P.write_sdp_input_file()
		P.run_sdp_solver(True)
		P.save("new-cod-1")


	elif prob == "k4-cod":
	
		P = AxiomsProblem()
		P.forbid_subgraph((4, 3))
		P.n = 7
		P.remove_types([3,4])
		P.save("k4-cod")
		#
		P.save("k4-cod")
		P.clear_axioms()
		P.add_codegree_axiom(Rational("1/4"))
		P._force_sharps = True
		C = RandomTournamentConstruction()
		P.use_construction(C)
		P.change_problem_bases()
		P.save("k4-cod")
		P.calculate_product_densities()
		P.save("k4-cod")
		P.write_sdp_input_file()
		P.save("k4-cod")
		P.run_sdp_solver(True)
		P.save("k4-cod")
		P.check_floating_point_bound()
		P.save("k4-cod")
		#P.make_exact(1024*1024)
		#P.check_exact_bound()


	elif prob == "marchant":
	
		P = AxiomsProblem()
		P.forbid_subgraph(Flag("5:123124125345"))
		P.n = 6
		P.remove_types([0,1,2,4])
		P.clear_axioms()
		P.add_codegree_axiom(Rational("1/3"))
		C = BlowupConstruction(Flag("3:112223331"))
		P.use_construction(C)
		P._sharp_graphs.extend([25, 27, 286, 289, 304, 389, 425])
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact(2**20)
		P.check_exact_bound()


	elif prob == "gammak4":
	
		P = AxiomsProblem()
		P.forbid_subgraph((4, 4))
		P.n = 6
		P.clear_axioms()
		P.add_codegree_axiom(Rational("1/2"))
		C = VariantRandomTournamentConstruction()
		P.use_construction(C)
		#._sharp_graphs.extend([25, 27, 286, 289, 304, 389, 425])
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver(True)
		#P.check_floating_point_bound()
		#P.change_solution_bases()
		#P.make_exact(2**20)
		#P.check_exact_bound()


	elif prob == "grzesik":

		P = Problem(2)
		P.forbid_subgraph((3, 3))
		P.n = 5
		P.set_density("5:1223344551")
		C = SymmetricBlowupConstruction(Flag("5:1223344551", 2))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "hirst":

		P = Problem(2)
		P.n = 7
		P.set_density("4:1223241314")
		C = BlowupConstruction(Flag("5:12131415232425343545", 2))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "paw":

		P = Problem(2)
		P.n = 5
		P.remove_types([4])
		
		P.set_density("4:12233114")
		C = BlowupConstruction(Flag("4:1223344111223344", 2))
		P.use_construction(C)
		P.change_problem_bases()
		
		P.add_zero_eigenvectors(0, matrix(QQ,[[1, 0, 0, 1/2, 31/70],
			[0, 1, 49/106, 7/108, -17/20]]))
		P.add_zero_eigenvectors(1, matrix(QQ,[[1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0]]))
		P.add_zero_eigenvectors(2, matrix(QQ,[[0,0,-5,0,8,0,0],[0,0,10,8,0,0,0]]))
		P.add_zero_eigenvectors(3, matrix(QQ,[[0,0,0,3,-1,0,0],[0,1,0,0,0,0,0]]))
		
		P._sharp_graphs.extend(eval("[0,4,11,18,19,24,27]"))
		P.change_problem_bases()
		
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver(False, sdpa="dd")
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "maxs3":

		P = Problem(2, oriented=True)
		P.n = 3
		P.set_density("3:1213")
		C = AdHocConstruction("maxs3")
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "maxs4":

		P = Problem(2, oriented=True)
		P.n = 4
		P.remove_types([0])
		P.set_density("4:121314")
		C = AdHocConstruction("maxs4")
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "34":

		P = Problem(2)
		P.forbid_induced_subgraph((4, 0))
		P.n = 5
		P.set_density((3, 3))
		C = BlowupConstruction(Flag("3:112233", 2))
		P.use_construction(C)
		P.change_problem_bases()
				
		P.add_zero_eigenvectors(2, matrix(QQ,[[0, 2, 1, 0, 0, 0, 0]]))
		P.add_zero_eigenvectors(3, matrix(QQ,[[1, 0, 1, 1, 0, 0, 0, 0]]))
		P.add_zero_eigenvectors(3, matrix(QQ,[[0, 0, 0, 0, 0, 0, 1, -1]]))
		P._sharp_graphs.extend(eval("[2, 3, 4, 16, 20]"))
		P.change_problem_bases()
		
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver(True)
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "43":


		P = Problem(2)
		P.forbid_induced_subgraph((3, 0))
		P.n = 6
		P.set_density((4, 6))
		C = BlowupConstruction(Flag("5:12233445511122334455", 2))
		P.use_construction(C)
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "44":

		P = Problem(2)
		P.forbid_subgraph((4, 6))
		P.n = 8
		P.remove_types([110, 132]) # these types have a full set of zero eigenvectors
		P.set_density_graph(Flag("4:",2))
		x = polygen(QQ)
		K = NumberField(x**3 - 2*x**2 + 2*x - Integer(2)/3, 'x', embedding=RDF(0.5))
		x = K.gen()
		#C = UnbalancedBlowupConstruction(Flag("8:122886633447755156781122334455667788",2),
		#	weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K)
		C = UnbalancedBlowupConstruction(Flag("8:131416171823242526273537384546485867",2),
			weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K)
		P.use_construction(C)
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver(True, "dd")
		return P, C

		#P.write_sdp_input_file()
		#P.run_sdp_solver()
		#P.check_floating_point_bound()
		#P.make_exact()
		#P.check_exact_bound()

		#P.remove_types([110,132])
		#P._sdp_Q_matrices = P._sdp_Q_matrices[0:110] + P._sdp_Q_matrices[111:132] + P._sdp_Q_matrices[133:]
		#P._solution_bases = P._solution_bases[0:110] + P._solution_bases[111:132] + P._solution_bases[133:]
		#P._inverse_solution_bases = P._inverse_solution_bases[0:110] + P._inverse_solution_bases[111:132] + P._inverse_solution_bases[133:]
		#P.make_exact(1024*1024)



	elif prob == "53":

		P = Problem(2)
		P.forbid_induced_subgraph((3, 0))
		P.n = 6
		P.set_density((5, binomial(5, 2)))
		C = BlowupConstruction(Flag("5:12233445511122334455", 2))
		P.use_construction(C)
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "63":

		P = Problem(2)
		P.forbid_induced_subgraph((3, 3))
		P.n = 7
		P.remove_types([1,2,4,7,9])
		P.set_density_graph(Flag("6:",2))
		C = SymmetricBlowupConstruction(ClebschGraph())
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver(True)
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "73":

		P = Problem(2)
		P.forbid_induced_subgraph((3, 3))
		P.n = 8
		P.remove_types([1,3,4,5,6,7,8,10])
		P.set_density_graph(Flag("7:",2))
		C = SymmetricBlowupConstruction(ClebschGraph())
 		P.use_construction(C)
 		P.change_problem_bases()
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver(True, True)
		P.check_floating_point_bound(tolerance=10e-11)
		P.make_exact(1024*1024,cholesky=range(40))
		#P.check_exact_bound()


	elif prob == "ch":
	
		P = AxiomsProblem(2, True)
		P.forbid_subgraph(Flag("3:122331",2,True))
		P.n = 6
		P.clear_axioms()
		P.add_out_degree_axiom(1/3)
		C = None
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver(True)

	return P,C
