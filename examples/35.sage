P = GraphProblem()
P.forbid_induced_subgraph((5, 0))
P.n = 6
P.set_density((3, 3))
C = BlowupConstruction(GraphFlag("4:11223344"), phantom_edge=[1,2])
P.set_extremal_construction(C)
P.set_inactive_types(0, 1, 2)
P.compute_products()
P.minimize = True
P.solve_sdp()
P.change_solution_bases()
P.make_exact(2^20)
P.check_exact_bound()