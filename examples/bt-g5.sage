P = ThreeGraphProblem()
P.forbid_subgraph("6:123124135145346256", "6:123124134125345136246", "6:123124134125135126136456")
P.forbid_homomorphic_images()
P.n = 6
C = BlowupConstruction(ThreeGraphFlag("6:123124125126134135146235246256345346356456"))
P.set_extremal_construction(C)
P.compute_products()
P.solve_sdp(True)
P.change_solution_bases()
P.make_exact(2^30)
P.check_exact_bound()


