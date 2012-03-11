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
		P.set_inv_anti_inv_bases()
		C = BlowupConstruction(Flag("3:123"))
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
		P.find_sharps()
	
	elif prob == "k4-":

		P = Problem()
		P.forbidden_edge_numbers={4:3}
		P.n = 6
		P.set_inv_anti_inv_bases()
		C = None
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()

	elif prob == "f32":

		P = Problem()
		P.forbidden_graphs=[Flag("5:123124125345")]
		P.n = 6
		P.set_inv_anti_inv_bases()
		C = BlowupConstruction(Flag("3:122123133"))
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()

	elif prob == "razb":
	
		P = Problem()
		P.forbidden_edge_numbers={4:4}
		P.forbidden_induced_graphs=[Flag("4:123")]
		P.n = 6
		P.set_inv_anti_inv_bases()
		C = BlowupConstruction(Flag("3:112223331123"))
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
	
	elif prob == "k5weak":
	
		P = Problem()
		P.forbidden_edge_numbers={5:10}
		h1 = Flag("5:134234125135235145245345")
		h2 = Flag("5:124134234125135235145245")
		P.forbidden_induced_graphs=[h1,h2]
		P.n = 6
		P.set_inv_anti_inv_bases()
		C = BlowupConstruction(Flag("2:112122"))
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
	
	elif prob == "38":

		P = Problem()
		P.forbidden_graphs=[Flag("5:123124125345"),Flag("5:123124125134135145")]
		P.n = 6
		P.set_inv_anti_inv_bases()
		C = BlowupConstruction(Flag("4:123124134234"))
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()

	elif prob == "max42":

		P = Problem()
		P.n = 5
		P.set_inv_anti_inv_bases()
		P.density_graph = Flag("4:123124")
		#C = random_geometric_construction()
		C = None
		#P.construction = C
		#P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()
	
	
	elif prob == "k4-f32":

		P = Problem()
		P.forbidden_edge_numbers={4:3}
		P.forbidden_graphs=[Flag("5:123124125345")]
		P.n = 7
		P.set_inv_anti_inv_bases()
		C = BlowupConstruction(Flag("6:123234345451512136246356256146"))
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		

	elif prob == "cod":
	
		P = AxiomsProblem()
		C = None
		P.forbidden_edge_numbers={4:3}
		P.n = 6
		#P.set_inv_anti_inv_bases()
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
		P.set_inv_anti_inv_bases()
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
		tg = Flag("2:")
		f1 = Flag("3:123", tg)
		f2 = Flag("2:", tg)
		P.add_axiom(tg, [f1, f2], [1, Rational("-1/4")])
		P.set_inv_anti_inv_bases()
		C = RandomTournamentConstruction()
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()


	elif prob == "marchant":
	
		P = axioms_problem()
		P.forbidden_graphs = [Flag("5:123124125345")]
		P.n = 6
		tg = Flag("2:")
		f1 = Flag("3:123", tg)
		f2 = Flag("2:", tg)
		P.add_axiom(tg, [f1, f2], [1, Rational("-1/3")])
		P.set_inv_anti_inv_bases()
		C = BlowupConstruction(Flag("3:112223331"))
		P.construction = C
		P.set_new_bases()
		P.calculate_product_densities()
		P.write_sdp_input_file()
		P.run_csdp()


	return P,C
