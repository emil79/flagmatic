P = GraphProblem(6, forbid_induced=(5,0), density=(3,3), minimize=True, type_orders=[3,4])
C = BlowupConstruction(GraphFlag("4:11223344"), phantom_edge=[1,2])
P.set_extremal_construction(C)
P.solve_sdp(solver="sdpa_dd")
P.make_exact(2^20)
graphs=[x[0] for x in C.subgraph_densities(6)]
[g.has_forbidden_graphs([GraphFlag("4:1223344113")],induced=True) for g in graphs]