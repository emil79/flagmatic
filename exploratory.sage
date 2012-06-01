


def write_certificate(problem, filename):

	problem.state("check_exact", "ensure_yes")

	description = E.problem._flag_cls.description() + "; "
	description += "minimize " if problem._minimize else "maximize "
	description += ", ".join(str(g) for g in problem._density_graphs) + " density"
	forbidden = []
	for g in problem._forbidden_graphs:
		forbidden.append(str(g))
	for g in problem._forbidden_induced_graphs:
		forbidden.append("induced %s" % g)
	for pair in problem._forbidden_edge_numbers:
		forbidden.append("induced %d.%d" % pair)
	if len(forbidden) > 0:
		description += "; forbid " + ", ".join(forbidden)

	def upper_triangular_matrix_to_list (M):
		return [list(M.row(i))[i:] for i in range(M.nrows())]

	def matrix_to_list (M):
		return [list(M.row(i)) for i in range(M.nrows())]

	if problem.state("diagonalize") == "yes":
		qdash_matrices = problem._exact_diagonal_matrices
		r_matrices = problem._exact_r_matrices
	else:
		qdash_matrices = problem._exact_Qdash_matrices
		r_matrices = problem._inverse_flag_bases

	data = {
		"description" : description,
		"bound" : problem._bound,
		"order_of_admissible_graphs" : problem._n,
		"number_of_admissible_graphs" : len(problem._graphs),
		"admissible_graphs" : problem._graphs,
		"admissible_graph_densities" : problem._densities[0], # TODO: multiple densities
		"number_of_types" : len(problem._types),
		"types" : problem._types,
		"numbers_of_flags" : [len(L) for L in problem._flags],
		"flags" : problem._flags,
		"qdash_matrices" : [upper_triangular_matrix_to_list(M) for M in qdash_matrices],
		"r_matrices" : [matrix_to_list(M) for M in r_matrices]
	}

	def default_handler (O):
		if O in ZZ:
			return int(Integer(O))
		return repr(O)

	try:
		with open(filename, "w") as f:
			json.dump(data, f, indent=4, default=default_handler)
		sys.stdout.write("Written certificate.\n")
	
	except IOError:
		sys.stdout.write("Cannot open file for writing.\n")


def show_eigenvalues(problem, ti):

	eigvals = sorted(numpy.linalg.eigvalsh(problem._sdp_Q_matrices[ti]))
	for i in range(len(eigvals)):
		sys.stdout.write("%d. %g\n" % (i + 1, eigvals[i]))


def show_zero_eigenvalues(problem, tolerance=1e-5, types=None):

	if types is None:
		types = problem._active_types

	for ti in types:
		eigvals = numpy.linalg.eigvalsh(problem._sdp_Q_matrices[ti])
		if tolerance is None:
			zero_eigvals = sorted(eigvals)
		else:
			zero_eigvals = sorted([e for e in eigvals if e < tolerance])
		if len(zero_eigvals) == 0:
			sys.stdout.write("Type %d. None.\n" % ti)
		else:
			sys.stdout.write("Type %d. %d possible: %s.\n" % (ti, len(zero_eigvals),
				" ".join("%s" % e for e in zero_eigvals)))


# TODO: transform with _flag_bases if present.

def check_construction(problem, C, types=None):

	if types is None:
		types = problem._active_types
	
	for ti in types:
		M = C.zero_eigenvectors(problem._types[ti], problem._flags[ti]) 
		if M.nrows() == 0:
			sys.stdout.write("Type %d. None.\n" % ti)
		else:
			R = M * problem._sdp_Q_matrices[ti]
			norms = [R[i,:].norm() for i in range(R.nrows())]
			sys.stdout.write("Type %d. %d possible: %s.\n" % (ti, len(norms),
				" ".join("%s" % e for e in norms)))


def find_extra_zero_eigenvectors(problem, ti, tolerance=1e-5, threshold=1e-10, echelonize=True, denominator=None):

	if not ti in problem._active_types:
		raise ValueError("Type is not active.")

	M = problem._zero_eigenvectors[ti].echelon_form()
	E = copy(get_zero_eigenvectors(problem, ti, tolerance=tolerance, use_bases=False))
	
	pivots = M.pivots()
	for i in range(len(pivots)):
		c = pivots[i]
		for r in range(E.nrows()):
			E[r, :] -= E[r, c] * M[i, :]

	if not echelonize:
		return E

	for r in range(E.nrows()):
		for c in range(E.ncols()):
			if abs(E[r, c]) > threshold:
				E[r, :] /= E[r, c]
				for s in range(E.nrows()):
					if s != r:
						E[s, :] -= E[s, c] * E[r, :]
				break

	for r in range(E.nrows()):
		for c in range(E.ncols()):
			if abs(E[r, c]) < threshold:
				 E[r, c] = 0

	r = E.nrows() - 1
	while r >= 0 and E[r, :].is_zero():
		E = E[:r, :]
		r -= 1

	if denominator is None:
		return E
	
	ER = matrix(QQ, E.nrows(), E.ncols())

	def rationalize (f):
		return Integer(round(f * denominator)) / denominator

	for r in range(E.nrows()):
		for c in range(E.ncols()):
			ER[r, c] = rationalize(E[r, c])

	return ER


def get_zero_eigenvectors(problem, ti, tolerance=1e-5, use_bases=True):

	if not ti in problem._active_types:
		raise ValueError("Type is not active.")

	if use_bases and problem.state("transform_solution") == "yes":
		QM = problem._sdp_Qdash_matrices[ti]
	else:
		QM = problem._sdp_Q_matrices[ti]

	ns = len(QM.subdivisions()[0]) + 1

	B = []
	for i in range(ns):
		QB = QM.subdivision(i, i)
		eigvals, T = numpy.linalg.eigh(QB)
		M = matrix(problem._approximate_field, 0, QB.ncols())
		for ei in range(len(eigvals)):
			if eigvals[ei] < tolerance:
				M = M.stack(matrix(numpy.matrix(T[:, ei])))
		B.append(M)
	return block_diagonal_matrix(B)


def easy_guess(problem, ti, guesses, target=None, tolerance=1e-8, column_threshold=1e-6):

	K = problem._field
	Q = matrix(RDF, problem._sdp_Q_matrices[ti])
	NZ = problem.get_zero_eigenvectors(ti, tolerance, use_bases=False)
	
	good_cols = [i for i in range(NZ.ncols()) if NZ[:,i].norm() > column_threshold]
	ngc = len(good_cols)
	sys.stdout.write("Using columns %s.\n" % good_cols)
	
	FZ = matrix(K, 0, Q.nrows())
	Z = matrix(K, 1, Q.nrows())
	
	for guess in guesses:
		l = len(guess)
		for tup in Tuples(good_cols, l):
			for i in range(l):
				Z[0, tup[i]] = guess[i]
			R = Z * Q
			norm = R.norm()
			if norm < tolerance:
				FZ = FZ.stack(Z)
				FZ.echelonize()
				if FZ[-1,:].is_zero():
					FZ = FZ[:-1,:]
				else:
					sys.stdout.write("%d. %.3g %s\n" % (FZ.nrows(), norm, Z.list()))
					if FZ.nrows() == target:
						return FZ
			for i in range(l):
				Z[0, tup[i]] = 0

	return FZ
	


def guess_zero_eigenvectors(P, ti, target=None):

	#x = polygen(QQ)
	#K = NumberField(x^3 - 2*x^2 + 2*x - 2/3, 'x', embedding=RDF(0.5))
	x = P._field.gen()
	a = x/4
	b = (1-x)/4
	af = RDF(a)
	bf = RDF(b)
	
	nf = len(P._flags[ti])
	Z = matrix(RDF, [[0.0 for i in range(nf)]]).T
	M = matrix(RDF, P._sdp_Q_matrices[ti])
	found_zero_eigenvectors = matrix(P._field, 0, nf)
	echelon_found_zero_eigenvectors = matrix(P._field, 0, nf)

	MZ = P.get_zero_eigenvectors(ti, tolerance=1e-12, use_bases=False)
	good_cols = [i for i in range(MZ.ncols()) if MZ[:,i].norm() > 1e-7]

	KZ = P._zero_eigenvectors[ti]
	known_cols = [i for i in range(KZ.ncols()) if not KZ[:,i].is_zero()]
	good_cols = known_cols

	ngc = len(good_cols)
	print nf, ngc

	
	total = IntegerVectors(4, max_length=ngc, min_length=ngc).cardinality() ** 2
	sys.stdout.write("%d possibilities to check.\n" % total)
	count = 0
	start_time = time.time()
	
	for p1 in IntegerVectors(4, max_length=ngc, min_length=ngc):
		for p2 in IntegerVectors(4, max_length=ngc, min_length=ngc):
			for i in range(ngc):
				Z[good_cols[i], 0] = (p1[i] * af) + (p2[i] * bf)
			R = M * Z
			norm = R.norm()
			if norm < 1e-10:
				NZ = matrix(P._field, 1, nf)
				for i in range(ngc):
					NZ[0, good_cols[i]] = (p1[i] * a) + (p2[i] * b)
				echelon_found_zero_eigenvectors = echelon_found_zero_eigenvectors.stack(NZ)
				echelon_found_zero_eigenvectors.echelonize()
				if echelon_found_zero_eigenvectors[-1, :].is_zero():
					echelon_found_zero_eigenvectors = echelon_found_zero_eigenvectors[:-1, :]
				else:
					found_zero_eigenvectors = found_zero_eigenvectors.stack(NZ)
					if not target is None and found_zero_eigenvectors.nrows() == target:
						return found_zero_eigenvectors
					print norm
				
				sys.stdout.write("+")
				sys.stdout.flush()
			count += 1
			#if count % 1000 == 0:
			#	sys.stdout.write(".")
			#	sys.stdout.flush()
			if count % 10000 == 0:
				time_elapsed = time.time() - start_time
				proportion_done = count * 1.0 / total
				sys.stdout.write("%3f%% done. (ETA: %.0f seconds)\n" % (proportion_done * 100, time_elapsed * (1 / proportion_done - 1)))
				
			
	return found_zero_eigenvectors


def simple_guess_zero_eigenvectors(P, ti, target=None):

	x = P._field.gen()
	a = x/4
	b = (1-x)/4
	af = RDF(a)
	bf = RDF(b)
	
	nf = len(P._flags[ti])
	Z = matrix(RDF, [[0.0 for i in range(nf)]]).T
	M = matrix(RDF, P._sdp_Q_matrices[ti])
	found_zero_eigenvectors = matrix(P._field, 0, nf)
	echelon_found_zero_eigenvectors = matrix(P._field, 0, nf)

	MZ = P.get_zero_eigenvectors(ti, tolerance=1e-12, use_bases=False)
	good_cols = [i for i in range(MZ.ncols()) if MZ[:,i].norm() > 1e-7]
	ngc = len(good_cols)
	print nf, ngc
	
	total = binomial(ngc, 8) * binomial(8, 4)
	sys.stdout.write("%d possibilities to check.\n" % total)
	count = 0
	start_time = time.time()
	
	for p1 in Combinations(range(ngc), 8):
		for p2 in Combinations(p1, 4):
			for i in range(ngc):
				if i in p2:
					Z[good_cols[i], 0] = af
				elif i in p1:
					Z[good_cols[i], 0] = bf
				else:
					Z[good_cols[i], 0] = 0
			R = M * Z
			norm = R.norm()
			if norm < 1e-10:
				NZ = matrix(P._field, 1, nf)
				for i in range(ngc):
					if i in p2:
						NZ[0, good_cols[i]] = a
					elif i in p1:
						NZ[0, good_cols[i]] = b
				echelon_found_zero_eigenvectors = echelon_found_zero_eigenvectors.stack(NZ)
				echelon_found_zero_eigenvectors.echelonize()
				if echelon_found_zero_eigenvectors[-1, :].is_zero():
					echelon_found_zero_eigenvectors = echelon_found_zero_eigenvectors[:-1, :]
				else:
					found_zero_eigenvectors = found_zero_eigenvectors.stack(NZ)
					if not target is None and found_zero_eigenvectors.nrows() == target:
						return found_zero_eigenvectors
					print norm
				
				sys.stdout.write("+")
				sys.stdout.flush()
			count += 1
			#if count % 1000 == 0:
			#	sys.stdout.write(".")
			#	sys.stdout.flush()
			if count % 100000 == 0:
				time_elapsed = time.time() - start_time
				proportion_done = count * 1.0 / total
				sys.stdout.write("%3f%% done. (ETA: %.0f seconds)\n" % (proportion_done * 100, time_elapsed * (1 / proportion_done - 1)))
				
			
	return found_zero_eigenvectors

