problem = ThreeGraphProblem(6, forbid=["6:123124135345146256346", "6:123124134125135126136456", "6:123124134125136256356456", "6:123124134125135145126136146156", "6:123124134234125135245236146346"], forbid_homomorphic_images=True)
x = polygen(QQ)
K.<x> = NumberField(x^2 - 13, embedding=3.6)
w1 = (5 - x)/6
w2 = (x - 2)/9
construction = ThreeGraphBlowupConstruction("5:123124134234125135235145245", weights=[w1, w1, w2, w2, w2], field=K)
problem.set_extremal_construction(construction)
problem.solve_sdp(True)
problem.make_exact(2^30)