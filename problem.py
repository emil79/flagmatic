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

import sys

class flagmatic_problem(object):

	forbidden_edge_numbers = {}
	forbidden_graphs = []
	forbidden_induced_graphs = []

	_flag_products = {}

	def __init__(self):
		self._n = 0
		self._graphs = []
		self._types = []
		self._flags = []
		
	@property
	def n(self):
		return self._n

	@n.setter
	def n(self, n):
		self._n = n

		sys.stdout.write("Generating graphs...\n")
		self._graphs = generate_graphs(n, forbidden_edge_numbers=self.forbidden_edge_numbers)
		sys.stdout.write("Generated %d graphs.\n" % len(self._graphs))
	
		sys.stdout.write("Generating types and flags...\n")
		num_types = 0
		self._flag_products = {}
	
		for s in range(n % 2, n - 1, 2):
			
			these_types = generate_graphs(s, forbidden_edge_numbers=self.forbidden_edge_numbers)
			sys.stdout.write("Generated %d types of order %d, " % (
				len(these_types), s))

			m = (n + s) / 2
			these_flags = []
			for tg in these_types:
				these_flags.append(generate_flags(m, tg, forbidden_edge_numbers=self.forbidden_edge_numbers))
			sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))

			for gi in range(len(self._graphs)):
				these_flags_products = flag_products(self._graphs[gi], s, m, these_types, these_flags)
				for key in these_flags_products.iterkeys():
					ti, fai, fbi = key
					self._flag_products[(gi, ti + num_types, fai, fbi)] = these_flags_products[key]
			
			self._types.extend(these_types)
			self._flags.extend(these_flags)
			num_types += len(these_types)

	@property
	def graphs(self):
		return self._graphs

	@property
	def types(self):
		return self._types

	@property
	def flags(self):
		return self._flags
	
	def write_dats(self):
	
		for key in sorted(self._flag_products.keys()):
			gi, ti, fai, fbi = key
			value = self._flag_products[key]
			sys.stdout.write("%d %d %d %d %s\n" % (gi + 1, ti + 2, fai + 1, fbi + 1,
				value.n(digits=64)))
	
