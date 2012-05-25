P = GraphProblem()
P.forbid_induced_subgraph((7, 0))
P.generate_flags(8, type_orders=[6])
P.set_density((3, 3))
C = BlowupConstruction(GraphFlag("6:112233445566"), phantom_edge=[1,2])
P.set_extremal_construction(C)
P.minimize = True
P.solve_sdp(True)
