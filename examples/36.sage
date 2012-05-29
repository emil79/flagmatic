problem = GraphProblem(7, forbid_induced=(6,0), density=(3,3), minimize=True, type_orders=[3,5])
construction = GraphBlowupConstruction("5:1122334455", phantom_edge=(1,2))
problem.set_extremal_construction(construction)
problem.solve_sdp(True, solver="sdpa_dd")
problem.make_exact(2^20)
graphs=[x[0] for x in construction.subgraph_densities(7)]
[g.has_forbidden_graphs([GraphFlag("4:1223344113")],induced=True) for g in graphs]