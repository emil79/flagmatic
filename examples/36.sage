P = GraphProblem()
P.forbid_induced_subgraph((6, 0))
P.n = 7
P.set_density((3, 3))
C = BlowupConstruction(GraphFlag("5:1122334455"), phantom_edge=[1,2])
P.set_extremal_construction(C)
P.set_inactive_types(0, 2)
P.compute_products()
P.minimize = True
P.solve_sdp(tolerance=3e-5)
P.change_solution_bases()
P.make_exact(2^40)
P.check_exact_bound()