P = ThreeGraphProblem()
P.forbid("6:123124135145346256", "6:123124134125345136246", "6:123124134125135126136456")
P.forbid_homomorphic_images()
P.generate_flags(6)
C = BlowupConstruction(ThreeGraphFlag("6:123124125126134135146235246256345346356456"))
P.set_extremal_construction(C)
P.solve_sdp(True)
P.make_exact(2^30)

