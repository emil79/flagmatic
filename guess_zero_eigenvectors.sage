
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

