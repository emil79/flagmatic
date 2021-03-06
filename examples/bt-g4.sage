problem = ThreeGraphProblem(7, forbid=["6:123124134156256", "7:123124134125126357367457467567", "7:123124345156257"], forbid_homomorphic_images=True, max_flags=190)
construction = ThreeGraphBlowupConstruction("5:123124125134135145", weights=[1/3,1/6,1/6,1/6,1/6])
problem.set_extremal_construction(construction)
problem.solve_sdp(True)
problem.make_exact(2^30)