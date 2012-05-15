"""

flagmatic 2

Copyright (c) 2012, E. R. Vaughan. All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1) Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

from flagmatic.all import *
import time

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

import cStringIO, datetime, os, sys, time, traceback

class Tee:

	def __init__(self, *args):
		self.files = list(args)
		self.start_time = time.time()
		self.needs_stamp = True

	def write(self, text):
		if self.needs_stamp:
			ts = "[%s] " % datetime.timedelta(seconds=int(time.time() - self.start_time))
			for f in self.files:
				f.write(ts)
		for f in self.files:
			f.write(text)
		self.needs_stamp = len(text) == 0 or text[-1] == "\n"

	def flush(self):
		for f in self.files:
			f.flush()


class Example:

	def __init__(self, name):

		P = None
		C = None
		base_filename = os.path.join("examples", name)
		
		buffer = cStringIO.StringIO()
		saved_stdout = sys.stdout
		sys.stdout = Tee(sys.stdout, buffer)
	
		with open(base_filename + ".sage", "r") as f:
			lines = f.readlines()
		
		# last line of file might not end in a newline
		if lines[-1][-1] != "\n":
			lines[-1] += "\n"
		
		lines.append('P.save("' + base_filename + '")\n')
		
		for line in lines:
			sys.stdout.write("sage: %s" % line)
			try:
				exec(preparse(line))
			except:
				traceback.print_exc(None, sys.stdout)
				break
		
		sys.stdout = saved_stdout
	
		self.problem = P
		self.construction = C
		self.output = buffer.getvalue()
	
		with open(base_filename + ".txt", "w") as f:
			f.write(self.output)
