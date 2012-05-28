P = ThreeGraphProblem()
P.forbid("6:123124135146156", "7:123124156346257", "7:123124156347567")
P.forbid_homomorphic_images()
P.generate_flags(7, type_orders=[1,2,3])
x = polygen(QQ)
K = NumberField(x^2 - 5, 'x', embedding=2.2)
x = K.gen()
w1 = (13 + 3*x)/62
w2 = (6 - x)/31
C = BlowupConstruction(ThreeGraphFlag("5:123124125345"), weights=[w1,w1,w2,w2,w2], field=K)
P.set_extremal_construction(C)
P.solve_sdp(True, tolerance=1e-8, solver="sdpa_dd")
P.make_exact(2^30)