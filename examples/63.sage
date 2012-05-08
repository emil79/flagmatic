P = GraphProblem()
P.forbid_induced_subgraph((3, 3))
P.n = 7
P.set_inactive_types(1, 2, 4, 7, 9)
P.set_density((6, 0))
C = BlowupConstruction(GraphFlag("g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde"))
P.set_extremal_construction(C)
P.compute_products()
P.minimize = True
P.solve_sdp(show_output=True, sdpa="dd")
P.change_solution_bases()
P.make_exact(2^30)
P.check_exact_bound()
