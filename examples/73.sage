P = GraphProblem()
P.forbid_induced_subgraph((3, 3))
P.n = 8
P.set_inactive_types(1, 3, 4, 5, 6, 7, 8, 10)
P.set_density((7, 0))
C = BlowupConstruction(GraphFlag("g:12131415162728292a373b3c3d484b4e4f595c5e5g6a6d6f6g7e7f7g8c8d8g9b9d9fabacaebgcfde"))
P.set_extremal_construction(C)
#P.add_zero_eigenvectors(0, matrix(QQ,[[611/4096, 225/1024, 135/512, 125/1024, 45/2048, 45/256, 195/4096]]), use_bases=False)
P.compute_products()
P.minimize = True
P.solve_sdp(show_output=True, sdpa="dd")
P.change_solution_bases()
P.make_exact(2^40)
P.check_exact_bound()