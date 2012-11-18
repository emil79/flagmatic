problem = ThreeGraphProblem(7, forbid=["5:123124125345", "4.3"], max_flags=150)
construction = ThreeGraphBlowupConstruction("6:123234345451512136246356256146")
problem.set_extremal_construction(construction)
problem.solve_sdp(True)
problem.make_exact(2^20)