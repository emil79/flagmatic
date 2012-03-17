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
		P.forbidden_graphs=[Flag("5:123124345")]
		P.n = 6
		P.create_block_bases()
		C = BlowupConstruction(Flag("3:123"))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
	
	elif prob == "k4-":

		P = Problem()
		P.forbidden_edge_numbers={4:3}
		P.n = 7
		#P.create_block_bases()
		C = None
		P.calculate_product_densities()
		#P.write_sdp_input_file()
		#P.run_csdp()


	elif prob == "f32":

		P = Problem()
		P.forbidden_graphs=[Flag("5:123124125345")]
		P.n = 6
		P.create_block_bases()
		C = BlowupConstruction(Flag("3:122123133"))
		P.use_construction(C)
		P._force_sharps = True
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "razb":
	
		P = Problem()
		P.forbidden_edge_numbers={4:4}
		P.forbidden_induced_graphs=[Flag("4:123")]
		P.n = 6
		P.create_block_bases()
		C = BlowupConstruction(Flag("3:112223331123"))
		P.use_construction(C)
		P._force_sharps = True
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()

	
	elif prob == "k5weak":
	
		P = Problem()
		P.forbidden_edge_numbers={5:10}
		h1 = Flag("5:134234125135235145245345")
		h2 = Flag("5:124134234125135235145245")
		P.forbidden_induced_graphs=[h1,h2]
		P.n = 6
		P.create_block_bases()
		C = BlowupConstruction(Flag("2:112122"))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "turank4-":
	
		P = Problem()
		P.forbidden_graphs=[Flag("4:123124134234")]
		P.n = 6
		P.create_block_bases()
		P.density_graph =Flag("4:123124134")
		C = BlowupConstruction(Flag("3:112223331123"))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "max59":
	
		P = Problem()
		P.n = 6
		P.create_block_bases()
		P.density_graph = Flag("5:123124125134135145234235245")
		C = BlowupConstruction(Flag("2:112122"))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "38":

		P = Problem()
		P.forbidden_graphs=[Flag("5:123124125345"),Flag("5:123124125134135145")]
		P.n = 6
		P.create_block_bases()
		C = BlowupConstruction(Flag("4:123124134234"))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()
		

	elif prob == "max42":

		P = Problem()
		P.n = 5
		P.create_block_bases()
		P.density_graph = Flag("4:123124")
		#C = random_geometric_construction()
		C = None
		#P.use_construction(C)
		#P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
	
	
	elif prob == "k4-f32":

		P = Problem()
		P.forbidden_edge_numbers={4:3}
		P.forbidden_graphs=[Flag("5:123124125345")]
		P.n = 7
		P.create_block_bases()
		C = BlowupConstruction(Flag("6:123234345451512136246356256146"))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		

	elif prob == "cod6":
	
		P = AxiomsProblem()
		C = None
		P.forbidden_edge_numbers={4:3}
		P.n = 6
		tg = Flag("2:")
		f1 = Flag("3:123", tg)
		f2 = Flag("2:", tg)
		P.add_axiom(tg, [f1, f2], [1, Rational("-264/1000")])
		P.create_block_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()

	elif prob == "codn":
	
		P = AxiomsProblem()
		P.forbidden_edge_numbers={4:3}
		P.forbidden_graphs = [Flag("6:612623634645651")]
		P.n = 6
		tg = Flag("2:")
		f1 = Flag("3:123", tg)
		f2 = Flag("2:", tg)
		P.add_axiom(tg, [f1, f2], [1, Rational("-1/4")])
		C = RandomTournamentConstruction()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()


	elif prob == "codr":
	
		P = AxiomsProblem()
		P.forbidden_edge_numbers={4:3}
		P.forbidden_graphs = [Flag("6:612623634645651")]
		P.n = 6
		P.create_block_bases()
		P.clear_axioms()
		P.add_codegree_axiom(Rational("1/4"))
		P._force_sharps = True
		C = RandomTournamentConstruction()
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact(1024*1024)
		P.check_exact_bound()

	elif prob == "k4-cod":
	
		P = AxiomsProblem()
		P.forbidden_edge_numbers={4:3}
		P.n = 7
		P.save_json("k4-cod.js")
		P.create_block_bases()
		P.save_json("k4-cod.js")
		P.clear_axioms()
		P.add_codegree_axiom(Rational("1/4"))
		P._force_sharps = True
		C = RandomTournamentConstruction()
		P.use_construction(C)
		P.set_new_bases()
		P.save_json("k4-cod.js")
		P.calculate_product_densities()
		P.save_json("k4-cod.js")
		P.write_sdp_input_file()
		P.save_json("k4-cod.js")
		P.run_csdp(True)
		P.save_json("k4-cod.js")
		P.check_floating_point_bound()
		P.save_json("k4-cod.js")
		#P.make_exact(1024*1024)
		#P.check_exact_bound()


	elif prob == "marchant":
	
		P = AxiomsProblem()
		P.forbidden_graphs = [Flag("5:123124125345")]
		P.n = 6
		P.add_codegree_axiom(Rational("1/3"))
		P.create_block_bases()
		P._force_sharps = True
		C = BlowupConstruction(Flag("3:112223331"))
		P.use_construction(C)
		P._target_bound = 0
		P.set_new_bases()
		P.calculate_product_densities()
		#P.write_sdp_input_file()
		#P.run_csdp()
		#P.check_floating_point_bound()
		#P.make_exact(1024*1024)
		#P.check_exact_bound()


	elif prob == "grzesik":

		P = Problem(2)
		P.forbidden_edge_numbers = {3:3}
		P.n = 5
		P.create_block_bases()
		P.density_graph = Flag("5:1223344551", 2)
		C = SymmetricBlowupConstruction(Flag("5:1223344551", 2))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "hirst":

		P = Problem(2)
		P.n = 7
		P.create_block_bases()
		P.density_graph = Flag("4:1223241314", 2)
		C = BlowupConstruction(Flag("5:12131415232425343545", 2))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "paw":

		P = Problem(2)
		P.n = 5
		
		P._types = P._types[:4]
		P._flags = P._flags[:4]
		P._flag_bases = P._flag_bases[:4]
		
		P.create_block_bases()
		P.density_graph = Flag("4:12233114", 2)
		C = BlowupConstruction(Flag("4:1223344111223344", 2))
		P.use_construction(C)
		P.set_new_bases()
		
		P.add_zero_eigenvectors(0, 0, matrix(QQ,[[70,0,0,35,31]]))
		P.add_zero_eigenvectors(1, 0, matrix(QQ,[[1,0,0,0],[0,0,0,1]]))
		P.add_zero_eigenvectors(2, 0, matrix(QQ,[[0,0,-5,0,8],[0,0,10,8,0]]))
		P.add_zero_eigenvectors(3, 0, matrix(QQ,[[0,0,0,5,-2],[0,1,0,0,0]]))
		P._sharp_graphs.extend([0,4,11,18,19,24, 27])
		P.set_new_bases()
		
		P.calculate_product_densities()
		P._force_sharps = True
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "34":

		P = Problem(2)
		P.forbidden_induced_graphs = [Flag("4:",2)]
		P.n = 5
		P.create_block_bases()
		P.density_graph = Flag("3:121323",2)
		C = BlowupConstruction(Flag("3:112233", 2))
		P.use_construction(C)
		P.set_new_bases()

		P.add_zero_eigenvectors(2, 0, matrix(QQ,[[0, 2, 1, 0, 0]]))
		P.add_zero_eigenvectors(3, 0, matrix(QQ,[[1, 0, 1, 1, 0, 0]]))
		P.add_zero_eigenvectors(3, 1, matrix(QQ,[[1, -1]]))
		P._sharp_graphs.extend([2, 3, 4, 16, 20])
		P.set_new_bases()

		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "43":


		P = Problem(2)
		P.forbidden_induced_graphs = [Flag("3:",2)]
		P.n = 6
		P.create_block_bases()
		P.density_graph = Flag("4:121314232434",2)
		C = BlowupConstruction(Flag("5:12233445511122334455", 2))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "53":


		P = Problem(2)
		P.forbidden_induced_graphs = [Flag("3:",2)]
		P.n = 6
		P.create_block_bases()
		P.density_graph = Flag("5:12131415232425343545",2)
		C = BlowupConstruction(Flag("5:12233445511122334455", 2))
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_csdp()
		P.check_floating_point_bound()
		P.make_exact()
		P.check_exact_bound()


	elif prob == "63":


		P = Problem(2)
		P.forbidden_induced_graphs = [Flag("3:121323",2)]
		P.n = 7
		P.create_block_bases()
		P.density_graph = Flag("6:",2)
		C = SymmetricBlowupConstruction(ClebschGraph())
		P.use_construction(C)
		P.set_new_bases()
		P.calculate_product_densities()
		P._minimize=True
		P.write_sdp_input_file()
		P.run_csdp(True)
		P.check_floating_point_bound()
		P.make_exact(1024*1024)
		P.check_exact_bound()

	return P,C
