problem = ThreeGraphProblem(7, forbid=["6:123124135146156", "7:123124156346257", "7:123124156347567"], forbid_homomorphic_images=True, type_orders=[1,2,3])
x = polygen(QQ)
K = NumberField(x^2 - 5, 'x', embedding=2.2)
x = K.gen()
w1 = (13 + 3*x)/62
w2 = (6 - x)/31
construction = BlowupConstruction(ThreeGraphFlag("5:123124125345"), weights=[w1,w1,w2,w2,w2], field=K)
problem.set_extremal_construction(construction)
problem.solve_sdp(True, tolerance=1e-8, solver="sdpa_dd")
problem.make_exact(2^30)