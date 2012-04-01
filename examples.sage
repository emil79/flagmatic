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


def example(prob):
	
	if prob == "ff83":
	
		P = Problem()
		P.forbid_subgraph(Flag("5:123124345"))
		P.n = 6
		C = BlowupConstruction(Flag("3:123"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P._approximate_field = RealField(113)
		P.write_sdp_input_file()
		P.run_sdp_solver(sdpa="qd")
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
	
	elif prob == "ff83n":
	
		P = Problem()
		P.forbid_subgraph(Flag("5:123124345"))
		P.n = 6
		C = BlowupConstruction(Flag("3:123"))
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact(1024)
		P.check_exact_bound()
	
	
	elif prob == "k4-":

		P = Problem()
		P.forbid_edge_number(4, 3)
		P.n = 7
		#
		C = None
		P.calculate_product_densities()
		#P.write_sdp_input_file()
		#P.run_sdp_solver()


	elif prob == "f32":

		P = Problem()
		P.forbid_subgraph(Flag("5:123124125345"))
		P.n = 6
		C = UnbalancedBlowupConstruction(Flag("2:122"), weights=[1,2])
		P.use_construction(C)
		P.change_problem_bases(use_blocks=True)
		P.calculate_product_densities()
		P._force_sharps = True
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()

	elif prob == "f32n":

		P = Problem()
		P.forbid_subgraph(Flag("5:123124125345"))
		P.n = 6
		C = UnbalancedBlowupConstruction(Flag("2:122"), weights=[1,2])
		P.use_construction(C)
		P.calculate_product_densities()
		P._force_sharps = True
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "razb":
	
		P = Problem()
		P.forbid_edge_number(4, 4)
		P.forbid_induced_edge_number(4, 1)
		P.n = 6
		
		C = BlowupConstruction(Flag("3:112223331123"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P._force_sharps = True
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()

	
	elif prob == "k5weak":
	
		P = Problem()
		P.forbid_edge_number(5, 10)
		P.forbid_induced_edge_number(5, 8)
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
		P.forbid_edge_number(4, 4)
		P.n = 6
		
		#P.set_density_graph(Flag("4:123124134"))
		P.set_density_edge_number(4, 3)
		C = BlowupConstruction(Flag("3:112223331123"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "max59":
	
		P = Problem()
		P.n = 6
		
		#P.set_density_graph(Flag("5:123124125134135145234235245"))
		P.set_density_edge_number(5, 9)
		C = BlowupConstruction(Flag("2:112122"))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
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
		P.forbid_edge_number(4, 3)
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
		P.forbid_edge_number(4, 3)
		P.n = 6
		tg = Flag("2:")
		f1 = Flag("3:123", tg)
		f2 = Flag("2:", tg)
		P.add_axiom(tg, [f1, f2], [1, Rational("-264/1000")])
		
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()

	elif prob == "codn":
	
		P = AxiomsProblem()
		P.forbid_edge_number(4, 3)
		P.forbid_subgraph(Flag("6:612623634645651"))
		P.n = 6
		P.clear_axioms()
		P.add_codegree_axiom(Rational("1/4"))
		C = RandomTournamentConstruction()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()


	elif prob == "codr":
	
		P1 = AxiomsProblem()
		P1.forbid_edge_number(4, 3)
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
		P.forbid_edge_number(4, 3)
		P.forbid_subgraph(Flag("6:612623634645651"))
		P.n = 6
		P.clear_axioms()
		P.add_codegree_axiom(Rational("1/4"))
		C = RandomTournamentConstruction()
		P.use_construction(C)
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.change_solution_bases()
		P.make_exact(1024*1024)
		P.check_exact_bound()



	elif prob == "k4-cod":
	
		P = AxiomsProblem()
		P.forbid_edge_number(4, 3)
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

#L=[0, 1, 2, 4, 5, 8, 9, 11, 16, 17, 18, 19, 25, 26, 28, 30, 33, 36, 37, 43, 44, 46, 49, 55, 56, 58, 61, 65, 66, 68, 73, 75, 79, 81, 83, 89, 90, 91, 94, 95, 99, 116, 118, 120, 127, 128, 131, 133, 136, 146, 168, 177, 183, 202, 209, 211, 214, 216, 218, 227, 229, 232, 245, 246, 247, 252, 260, 261, 262, 263, 268, 270, 277, 278, 285, 289, 301, 303, 314, 315, 322, 330, 331, 334, 337, 344, 346, 347, 462, 467, 469, 484, 488, 490, 491, 494, 503, 511, 520, 521, 529, 549, 565, 572, 575, 578, 579, 580, 582, 586, 591, 601, 607, 609, 615, 617, 618, 620, 621, 622, 624, 797, 815, 820, 821, 822, 825, 833, 834, 836, 838, 851, 852, 853, 872, 878, 879, 882, 883, 884, 887, 889, 975, 978, 981, 983, 1046, 1048, 1052, 1053, 1061, 1062, 1131, 1133, 1190, 1201, 1206, 1212, 1213, 1216, 1230, 1235, 1260, 1269, 1284, 1288, 1289, 1357, 1397, 1408, 1419, 1420, 1516, 1517, 1573, 1581, 1610, 1628, 1640]


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
		P.forbid_edge_number(4, 4)
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
		P.forbid_edge_number(3, 3)
		P.n = 5
		P.set_density_graph(Flag("5:1223344551", 2))
		C = SymmetricBlowupConstruction(Flag("5:1223344551", 2))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver(True)
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "hirst":

		P = Problem(2)
		P.n = 7
		P.set_density_graph(Flag("4:1223241314", 2))
		C = BlowupConstruction(Flag("5:12131415232425343545", 2))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "paw":

		P = Problem(2)
		P.n = 5
		P.remove_types([4])
		
		P.set_density_graph(Flag("4:12233114", 2))
		C = BlowupConstruction(Flag("4:1223344111223344", 2))
		P.use_construction(C)
		P.change_problem_bases()
		
		P.add_zero_eigenvectors(0, matrix(QQ,[[1, 0, 0, '1/2', '31/70'],
			[0, 1, '49/106', '7/108', '-17/20']]))
		
		#P.add_zero_eigenvectors(1, matrix(QQ,[[1,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0]]))
		#P.add_zero_eigenvectors(2, matrix(QQ,[[0,0,-5,0,8,0,0],[0,0,10,8,0,0,0]]))
		#P.add_zero_eigenvectors(3, matrix(QQ,[[0,0,0,5,-2,0,0],[0,1,0,0,0,0,0]]))
		P._sharp_graphs.extend([0,4,11,18,19,24, 27])
		P.change_problem_bases()
		
		P.calculate_product_densities()
		P._force_sharps = True
		P._approximate_field = RealField(113)
		P.write_sdp_input_file()
		P.run_sdp_solver(True, sdpa="qd")
		P.check_floating_point_bound()
		P.make_exact(2**30)
		P.check_exact_bound()

	elif prob == "paw2":
		
		P = Problem(2)
		P.n = 5
		P.remove_types([4])
		P.set_density_graph(Flag("4:12233114", 2))
		C = BlowupConstruction(Flag("4:1223344111223344", 2))
		P.use_construction(C)
		P.calculate_product_densities()
		P._approximate_field = RealField(113)
		P.write_sdp_input_file()
		P.run_sdp_solver(True, sdpa="qd")
		P._sharp_graphs.extend([24, 27])
		P.check_floating_point_bound()
		P.make_exact(2**30)
		

	elif prob == "maxs3":

		P = Problem(2, oriented=True)
		P.n = 3
		P.set_density_graph(Flag("3:1213", 2, oriented=True))
		C = AdHocConstruction("maxs3")
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver(True)
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()

	elif prob == "maxs4":

		P = Problem(2, oriented=True)
		P.n = 4
		P.remove_types([0])
		P.set_density_graph(Flag("4:121314", 2, oriented=True))
		C = AdHocConstruction("maxs4")
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "34":

		P = Problem(2)
		P.forbid_induced_edge_number(4, 0)
		P.n = 5
		
		P.set_density_graph(Flag("3:121323",2))
		C = BlowupConstruction(Flag("3:112233", 2))
		P.use_construction(C)
		P.change_problem_bases()
				
		P.add_zero_eigenvectors(2, matrix(QQ,[[0, 2, 1, 0, 0, 0, 0]]))
		P.add_zero_eigenvectors(3, matrix(QQ,[[1, 0, 1, 1, 0, 0, 0, 0]]))
		P.add_zero_eigenvectors(3, matrix(QQ,[[0, 0, 0, 0, 0, 0, 1, -1]]))
		P._sharp_graphs.extend([2, 3, 4, 16, 20])
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
		P.forbid_induced_edge_number(3, 0)
		P.n = 6
		
		P.set_density_graph(Flag("4:121314232434",2))
		C = BlowupConstruction(Flag("5:12233445511122334455", 2))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "44":

		P = Problem(2)
		P.forbid_edge_number(4, 6)
		P.n = 8
		#P.remove_types([110, 132]) # these types have a full set of zero eigenvectors
		
		P.set_density_graph(Flag("4:",2))
		x = polygen(QQ)
		K = NumberField(x**3 - 2*x**2 + 2*x - Integer(2)/3, 'x', embedding=RDF(0.5))
		x = K.gen()
		#C = UnbalancedBlowupConstruction(Flag("8:122886633447755156781122334455667788",2),
		#	weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K)
		C = UnbalancedBlowupConstruction(Flag("8:131416171823242526273537384546485867",2),
			weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K)
		P.use_construction(C)
		#P.calculate_product_densities()
		P._minimize=True
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
		P.forbid_induced_edge_number(3, 0)
		P.n = 6
		
		P.set_density_graph(Flag("5:12131415232425343545",2))
		C = BlowupConstruction(Flag("5:12233445511122334455", 2))
		P.use_construction(C)
		P.change_problem_bases()
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_sdp_solver()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "63":

		P = Problem(2)
		P.forbid_induced_edge_number(3, 3)
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
		P.forbid_induced_edge_number(3, 3)
		P.n = 8
		P.remove_types([1,3,4,5,6,7,8,10])
#		
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
		#P.forbid_subgraph(Flag("4:122334",2,True))
		P.n = 5
		P.clear_axioms()
		P.add_out_degree_axiom(Integer(34)/100)
		
		#C = BlowupConstruction(Flag("3:121323",2,True))
 		#P.use_construction(C)
 		P.change_problem_bases()
		P.calculate_product_densities()
		#P._force_sharps = True
		P.write_sdp_input_file()
		P.run_sdp_solver(True,sdpa="qd")
		C = None

	return P,C
