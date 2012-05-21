P = GraphProblem()
P.forbid_induced_subgraph((6, 0))
P.generate_flags(7, type_orders=[3,5])
P.set_density((3, 3))
C = BlowupConstruction(GraphFlag("5:1122334455"), phantom_edge=[1,2])
P.set_extremal_construction(C)
P.minimize = True
P.solve_sdp(True, solver="sdpa_dd")
P.change_solution_bases()
P.make_exact(2^20)
P.check_exact_bound()
graphs=[x[0] for x in C.subgraph_densities(7)]
[g.has_forbidden_graphs([GraphFlag("4:1223344113")],induced=True) for g in graphs]