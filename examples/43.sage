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