problem = GraphProblem(7, forbid=(3,3), type_orders=[3,5], density=["7:12132435465767", "7:1213142526354767", "7:121314252635374667", (2, "7:121314252635364757"), (4, "7:12131425263536454767")])
construction = GraphBlowupConstruction("5:1223344551")
problem.set_extremal_construction(construction)
problem.solve_sdp()
problem.make_exact(2^30)