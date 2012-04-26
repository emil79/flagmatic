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
		P.compute_products()
		P.write_sdp_input_file()
		P.run_sdp_solver()


def ClebschGraph():

	cleb = GraphFlag()
	cleb.n = 16
	edges = [(1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (2, 7), (2, 8), (2, 9), (2, 10),
		(3, 7), (3, 11), (3, 12), (3, 13), (4, 8), (4, 11), (4, 14), (4,15), (5, 9),
		(5, 12), (5, 14), (5, 16), (6, 10), (6, 13), (6, 15), (6, 16), (7, 14), (7, 15),
		(7, 16), (8, 12), (8, 13), (8, 16), (9, 11), (9, 13), (9, 15), (10, 11), (10, 12),
		(10, 14), (11, 16), (12, 15), (13, 14)]
	for e in edges:
		cleb.add_edge(e)
	return cleb



def example(prob):
	
	if prob == "ff83-cpb":
	
		P = ThreeGraphProblem()
		P.forbid_subgraph("5:123124345")
		P.n = 6
		C = BlowupConstruction(ThreeGraphFlag("3:123"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.change_problem_bases()
		P._approximate_field = RealField(113) # for testing - not needed.
		P.solve_sdp()
		P.make_exact()
		P.check_exact_bound()
	
	elif prob == "ff83":
	
		P = ThreeGraphProblem()
		P.forbid_subgraph("5:123124345")
		P.n = 6
		C = BlowupConstruction(ThreeGraphFlag("3:123"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
	
	elif prob == "k4-":

		P = ThreeGraphProblem()
		P.forbid_subgraph((4, 3))
		P.n = 6
		C = None
		P.compute_products()
		P.write_sdp_input_file()
		P.run_sdp_solver()

	elif prob == "c5axiom":

		P = ThreeGraphAxiomsProblem()
		P.forbid_subgraph("5:123234345451512")
		P.n = 6
		P.add_degree_axiom(2*3^(1/2)-3)
		C = None
		P.compute_products()
		P.write_sdp_input_file()
		P.run_sdp_solver(True)

	elif prob == "f32":

		P = ThreeGraphProblem()
		P.forbid_subgraph("5:123124125345")
		P.n = 6
		C = BlowupConstruction(ThreeGraphFlag("2:122"), weights=[1,2])
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "razb":
	
		P = ThreeGraphProblem()
		P.forbid_subgraph((4, 4))
		P.forbid_induced_subgraph((4, 1))
		P.n = 6
		C = BlowupConstruction(ThreeGraphFlag("3:112223331123"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()

	
	elif prob == "k5weak":
	
		P = ThreeGraphProblem()
		P.forbid_subgraph((5, 10))
		P.forbid_induced_subgraph((5, 8))
		P.n = 6
		C = BlowupConstruction(ThreeGraphFlag("2:112122"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.change_problem_bases()
		P.solve_sdp()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "turank4-":
	
		P = ThreeGraphProblem()
		P.forbid_subgraph((4, 4))
		P.n = 6
		P.set_density((4,3))
		C = BlowupConstruction(ThreeGraphFlag("3:112223331123"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "max59":
	
		P = ThreeGraphProblem()
		P.n = 6
		P.set_density((5, 9))
		C = BlowupConstruction(ThreeGraphFlag("2:112122"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "max56":
	
		P = ThreeGraphProblem()
		P.n = 6
		P.set_density_edge_number(5, 6)
		C = BlowupConstruction(ThreeGraphFlag("3:112223331122233311"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.change_problem_bases()
		P.solve_sdp()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "38":
	
		P = ThreeGraphProblem()
		P.forbid_subgraph("5:123124125345", "5:123124125134135145")
		P.n = 6
		C = BlowupConstruction(ThreeGraphFlag("4:123124134234"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "max42":

		P = ThreeGraphProblem()
		P.n = 5
		P.remove_types([0])
		P.set_density("4.2")
		C = AdHocConstruction("max42")
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
	
	elif prob == "k4-f32":

		P = ThreeGraphProblem()
		P.forbid_subgraph((4, 3))
		P.forbid_subgraph(ThreeGraphFlag("5:123124125345"))
		P.n = 7
		P.set_inactive_types(3, 4, 5, 6)
		C = BlowupConstruction(ThreeGraphFlag("6:123234345451512136246356256146"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.import_solution("../output/k4-f32")
		P.save("k4-f32")
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.save("k4-f32")
		P.make_exact(2^23, protect=[2])
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
		P.compute_products()
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
		P1.set_extremal_construction(C)
		P1.compute_products()
		P1.change_problem_bases()
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

		P = ThreeGraphAxiomsProblem()
		P.forbid_subgraph((4, 3), "6:612623634645651")
		P.n = 6
		P.clear_densities()
		P.add_codegree_axiom(1/4, False)
		C = RandomTournamentConstruction()
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
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

		P.set_extremal_construction(C)
		P.compute_products()
		P.save("new-cod-1")
		P._force_sharps = True
		P.save("new-cod-1")
		P.solve_sdp(True, sdpa="dd", tolerance=1e-10)
		P.save("new-cod-1")


	elif prob == "k4-cod":
	
		P = ThreeGraphAxiomsProblem()
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
		P.set_extremal_construction(C)
		P.compute_products()
		P.change_problem_bases()
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
	
		P = ThreeGraphAxiomsProblem()
		P.forbid_subgraph("5:123124125345")
		P.n = 6
		P.set_inactive_types(0, 1, 2, 4)
		P.clear_densities()
		P.add_codegree_axiom(1/3, False)
		C = BlowupConstruction(ThreeGraphFlag("3:112223331"))
		P.set_extremal_construction(C)
		P.add_sharp_graphs(25, 27, 286, 289, 304, 389, 425)
		P.compute_products()
		P.solve_sdp(sdpa="dd", tolerance=1e-10)
		P.change_solution_bases()
		P.make_exact(2^20)
		P.check_exact_bound()


	elif prob == "gammak4":
	
		P = ThreeGraphAxiomsProblem()
		P.forbid_subgraph((4, 4))
		P.n = 6
		P.add_codegree_axiom(1/2)
		C = RandomTournamentConstruction(True)
		P.set_extremal_construction(C)
		P.compute_products()
		P.write_sdp_input_file()
		P.run_sdp_solver(True)

		#P.check_floating_point_bound()
		#P.change_solution_bases()
		#P.make_exact(2^20)
		#P.check_exact_bound()


	elif prob == "grzesik":

		P = GraphProblem()
		P.forbid_subgraph((3, 3))
		P.n = 5
		P.set_density("5:1223344551")
		C = BlowupConstruction(GraphFlag("5:1223344551"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "hirst":

		P = GraphProblem()
		P.n = 7
		P.set_density("4:1223241314")
		C = BlowupConstruction(GraphFlag("5:12131415232425343545"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "paw":

		P = GraphProblem()
		P.n = 5
		P.set_inactive_types(4)
			
		P.set_density("4:12233114")
		C = BlowupConstruction(GraphFlag("4:1223344111223344"))
		P.set_extremal_construction(C)
		P.compute_products()
		
		P.compute_flag_bases()
		
		P.add_zero_eigenvectors(0, matrix(QQ,[[1, 0, 0, 1/2, 31/70],
			[0, 1, 49/106, 7/108, -17/20]]), use_bases=True)
		P.add_zero_eigenvectors(1, matrix(QQ,[[1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0]]), use_bases=True)
		P.add_zero_eigenvectors(2, matrix(QQ,[[0,0,-5,0,8,0,0],[0,0,10,8,0,0,0]]), use_bases=True)
		P.add_zero_eigenvectors(3, matrix(QQ,[[0,0,0,3,-1,0,0],[0,1,0,0,0,0,0]]), use_bases=True)
		
		P.change_problem_bases()
		
		P.add_sharp_graphs(0, 4, 11, 18, 19, 24, 27)
		
		P.solve_sdp()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "maxs3":

		P = OrientedGraphProblem()
		P.n = 3
		P.set_density("3:1213")
		C = AdHocConstruction("maxs3")
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "maxs4":

		P = OrientedGraphProblem()
		P.n = 4
		P.set_inactive_types(0)
		P.set_density("4:121314")
		C = AdHocConstruction("maxs4")
		P.set_extremal_construction(C)
		P.compute_products()
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "34old":

		P = GraphProblem()
		P.forbid_induced_subgraph((4, 0))
		P.n = 5
		P.set_density((3, 3))
		C = BlowupConstruction(GraphFlag("3:112233"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.compute_flag_bases()

		P.add_zero_eigenvectors(2, matrix(QQ,[[0, 2, 1, 0, 0, 0, 0]]), use_bases=True)
		P.add_zero_eigenvectors(3, matrix(QQ,[[1, 0, 1, 1, 0, 0, 0, 0]]), use_bases=True)
		P.add_zero_eigenvectors(3, matrix(QQ,[[0, 0, 0, 0, 0, 0, 1, -1]]), use_bases=True)
		P.add_sharp_graphs(2, 3, 4, 16, 20)
		P.change_problem_bases()
		
		P.minimize=True
		P.solve_sdp()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "34":

		P = GraphProblem()
		P.forbid_induced_subgraph((4, 0))
		P.n = 5
		P.set_density((3, 3))
		C = BlowupConstruction(GraphFlag("3:112233"))
		P.set_extremal_construction(C)
		P.compute_products()
		
		P.add_zero_eigenvectors(2, matrix(QQ, [[0, 1/3, 1/3, 1/3, 0, 0, 0, 0]]))
		P.add_zero_eigenvectors(3, matrix(QQ,
			[[1/3, 0, 0, 1/3, 1/3, 0, 0, 0],
			[1/3, 0, 1/3, 0, 0, 1/3, 0, 0]]))
		P.add_sharp_graphs(2, 3, 4, 16, 20)
				
		P.minimize = True
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "34pe":

		P = GraphProblem()
		P.forbid_induced_subgraph((4, 0))
		P.n = 5
		P.set_density((3, 3))
		C = BlowupConstruction(GraphFlag("3:112233"), phantom_edge=[1,2])
		P.set_extremal_construction(C)
		P.compute_products()
		
		P.add_sharp_graphs(2, 3, 4, 16, 20)
				
		P.minimize = True
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "43":


		P = GraphProblem()
		P.forbid_induced_subgraph((3, 0))
		P.n = 6
		P.set_density((4, 6))
		C = BlowupConstruction(GraphFlag("5:12233445511122334455"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.minimize=True
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "43i":


		P = GraphProblem()
		P.forbid_induced_subgraph((3, 0))
		P.n = 6
		P.set_density((4, 6))
		C = BlowupConstruction(GraphFlag("5:12233445511122334455"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.minimize=True
		P.import_solution("/Users/emil/Projects/flagmatic-1.5-mac64/tmp", complement=True)
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "44":

		P = GraphProblem()
		P.forbid_subgraph((4, 6))
		P.n = 8
		P.save("44")
		#P.set_inactive_types(110, 132) # these types have a full set of zero eigenvectors
		P.set_density((4, 0))
		x = polygen(QQ)
		K = NumberField(x**3 - 2*x**2 + 2*x - Integer(2)/3, 'x', embedding=RDF(0.5))
		x = K.gen()
		#C = BlowupConstruction(GraphFlag("8:122886633447755156781122334455667788"),
		#	weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K)
		C = BlowupConstruction(GraphFlag("8:131416171823242526273537384546485867"),
			weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K)
		P.set_extremal_construction(C)
		P.save("44")
		P.compute_products()
		P.minimize=True
		P.save("44")
		P._approximate_field = RealField(113)
		P.import_solution("/Users/emil/Projects/flagmatic/oleg/44")
		P.check_floating_point_bound(tolerance=1e-9)

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



	elif prob == "44i":

		P = GraphProblem()
		P.forbid_induced_subgraph((4, 0))
		P.n = 8
		P.set_inactive_types(110, 132) # these types have a full set of zero eigenvectors
		P.set_density((4, 6))
		P.minimize = True
		P.save("44i")
		
		x = polygen(QQ)
		K = NumberField(x^3 - 2*x^2 + 2*x - 2/3, 'x', embedding=RDF(0.5))
		x = K.gen()
		C = BlowupConstruction(GraphFlag("8:122886633447755156781122334455667788"),
			weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K)
		
		P.set_extremal_construction(C)
		P.save("44i")
		P.compute_products()

		P.save("44i")
		P._approximate_field = RealField(113)
		P.import_solution("/Users/emil/Projects/flagmatic/oleg/44", complement=True)
		P.save("44i")
		P.check_floating_point_bound(tolerance=1e-9)
		P.save("44i")
		

	elif prob == "53":

		P = GraphProblem()
		P.forbid_induced_subgraph((3, 0))
		P.n = 6
		P.set_density((5, binomial(5, 2)))
		C = BlowupConstruction(GraphFlag("5:12233445511122334455"))
		P.set_extremal_construction(C)
		P.compute_products()
		P.minimize = True
		P.solve_sdp()
		P.change_solution_bases()
		P.make_exact(2^20)
		P.check_exact_bound()


	elif prob == "63":

		P = GraphProblem()
		P.forbid_induced_subgraph((3, 3))
		P.n = 7
		P.set_inactive_types(1, 2, 4, 7, 9)
		P.set_density((6, 0))
		C = BlowupConstruction(ClebschGraph())
		P.set_extremal_construction(C)
		P.compute_products()
		P.minimize = True
		P.solve_sdp(show_output=True)
		P.change_solution_bases()
		P.make_exact(2^20)
		P.check_exact_bound()


	elif prob == "63-cpb":

		P = GraphProblem()
		P.forbid_induced_subgraph((3, 3))
		P.n = 7
		P.set_inactive_types(1, 2, 4, 7, 9)
		P.set_density((6, 0))
		C = SymmetricBlowupConstruction(ClebschGraph())
		P.set_extremal_construction(C)
		P.compute_products()
		P.change_problem_bases()
		P.minimize = True
		P.solve_sdp(show_output=True)
		P.make_exact()
		P.check_exact_bound()


	elif prob == "73":

		P = GraphProblem()
		P.forbid_induced_subgraph((3, 3))
		P.n = 8
		P.set_inactive_types(1, 3, 4, 5, 6, 7, 8, 10)
		P.set_density((7, 0))
		C = SymmetricBlowupConstruction(ClebschGraph())
 		P.set_extremal_construction(C)
		P.compute_products()
		P.minimize = True
		P.save("73")
		P.solve_sdp(show_output=True, sdpa="dd", tolerance=1e-12)
		P.save("73")
 		P.change_solution_bases()
		P.make_exact(2**30)
		P.check_exact_bound()


	elif prob == "ch":
	
		P = OrientedGraphAxiomsProblem()
		P.forbid_subgraph("3:122331")
		P.n = 6
		P.add_out_degree_axiom(1/3)
		C = None
		P.compute_products()
		P.write_sdp_input_file()
		P.run_sdp_solver(True)

	return P,C
