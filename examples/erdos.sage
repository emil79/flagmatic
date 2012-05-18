P = GraphAxiomsProblem()
n = 6
target = Integer(1)/14       # can do 1/14 with 6 and 1/15 with 7

P.forbid_subgraph((3, 3))
P.generate_flags(n)
P.clear_densities()

tgs = []
for s in range(0, n-1): # n-2
	tgs.extend(GraphFlag.generate_graphs(s, forbidden_edge_numbers=[(3, 3)]))

for tg in tgs:

	s = tg.n
	ctg = copy(tg)
	ctg.t = s
			
	vertices = range(1, s + 1)
	fgs = GraphFlag.generate_flags(s + 2, tg, forbidden_edge_numbers=[(3, 3)])
	fgs = [fg for fg in fgs if (s + 1, s + 2) in fg.edges]

	mfgs = GraphFlag.generate_flags(s + 3, tg, forbidden_edge_numbers=[(3, 3)])
	mfgs = [fg for fg in fgs if (s + 1, s + 2) in fg.edges and (s + 1, s + 3) in fg.edges]

	partitions = [[vertices, []]]
	if len(vertices) > 1:
		partitions.extend(list(SetPartitions(vertices, 2)))

	#print partitions

	for partition in partitions:
		red_vertices = list(partition[0])
		blue_vertices = list(partition[1])

		if all((e[0] in red_vertices and e[1] in blue_vertices) or
			(e[0] in blue_vertices and e[1] in red_vertices) for e in tg.edges) and all(
			(x, y) in tg.edges or (y, x) in tg.edges for x in red_vertices for y in blue_vertices):
		
			use_flags = []
			for fg in fgs:
				rn1 = [v for v in red_vertices if (v, s + 1) in fg.edges]
				rn2 = [v for v in red_vertices if (v, s + 2) in fg.edges]
				bn1 = [v for v in blue_vertices if (v, s + 1) in fg.edges]
				bn2 = [v for v in blue_vertices if (v, s + 2) in fg.edges]
				if ((len(rn1) > 0 and len(rn2) > 0) or (len(bn1) > 0 and len(bn2) > 0) or
					(len(rn1) == len(bn1) == len(rn2) == len(bn2) == 0)):
					use_flags.append(fg)

			#print tg, list(red_vertices), blue_vertices, use_flags
		
			sys.stdout.write("Axiom: ")
			for fg in use_flags:
				sys.stdout.write("%s + " % fg)
			sys.stdout.write(" >= %s * %s\n" % (target, ctg))
			
			axiom = [(ctg, -target)]
			for fg in use_flags:
				axiom.append((fg, 1))
			P.add_axiom(tg, axiom, False)
			
			use_flags = []
			for fg in mfgs:
				if (any((v, s + 1) in fg.edges for v in red_vertices) and any(
					(v, s + 2) in fg.edges for v in red_vertices)) or (
					any((v, s + 1) in fg.edges for v in blue_vertices) and any(
					(v, s + 2) in fg.edges for v in blue_vertices)) or (
					not any((v, s + 1) in fg.edges or (v, s + 2) in fg.edges for v in vertices)):
					use_flags.append(fg)

			#print tg, list(red_vertices), blue_vertices, use_flags
		
			sys.stdout.write("Axiom: ")
			for fg in use_flags:
				sys.stdout.write("%s + " % fg)
			sys.stdout.write(" >= %s * %s\n" % (target, ctg))
			
			axiom = [(ctg, -target)]
			for fg in use_flags:
				axiom.append((fg, 1))
			P.add_axiom(tg, axiom, False)

#P.add_axiom(GraphFlag("1:"), [(GraphFlag("2:12(1)"), 1), (GraphFlag("1:(1)"), -target)], False)

#P.add_axiom(GraphFlag("2:12"), [(GraphFlag("4:1234(2)"), 1), (GraphFlag("2:12(2)"), -target)], False)
#P.add_axiom(GraphFlag("1:"), [(GraphFlag("3:23(1)"), 1), (GraphFlag("1:(1)"), -target)], False)
#P.add_axiom(GraphFlag("1:"), [(GraphFlag("2:12(1)"), 1), (GraphFlag("1:(1)"), -2*target)], False)

P.solve_sdp(True)
