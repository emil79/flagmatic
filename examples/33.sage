P = GraphProblem()
P.forbid_induced_subgraph((3, 0))
P.n = 5
P.set_density((3, 3))
C = BlowupConstruction(GraphFlag("2:1122"), phantom_edge=[1,2])
P.set_extremal_construction(C)
P.minimize = True
P.solve_sdp()
P.change_solution_bases()
P.make_exact()
P.check_exact_bound()