problem = ThreeGraphProblem(6, forbid=["6:123124135145346256", "6:123124134125345136246", "6:123124134125135126136456"], forbid_homomorphic_images=True)
construction = ThreeGraphBlowupConstruction("6:123124125126134135146235246256345346356456")
problem.set_extremal_construction(construction)
problem.solve_sdp(True)
problem.make_exact(2^30)

