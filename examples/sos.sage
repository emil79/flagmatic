problem = ThreeGraphAxiomsProblem(7, forbid=(4, 3))
problem.clear_densities()
tg = ThreeGraphFlag("2:")
ff = ThreeGraphFlag.generate_flags(5, tg, forbidden_edge_numbers=[(4,3), (4,4)])
ff1 = [g for g in ff if not ((1,2,3) in g.edges or (1,2,4) in g.edges or (1,2,5) in g.edges) and (3,4,5) in g.edges]
ff2 = [g for g in ff if not ((1,2,3) in g.edges or (1,2,4) in g.edges or (1,2,5) in g.edges or (3,4,5) in g.edges)]
axiom = [(g, 3) for g in ff1] + [(g, -1) for g in ff2]
problem.add_axiom(tg, axiom, False)
construction = RandomTournamentConstruction()
problem.set_extremal_construction(construction)
problem.solve_sdp(True, use_initial_point=True, force_sharp_graphs=True)
#problem.make_exact(2^20)
