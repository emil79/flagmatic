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
		self._graph_block = make_graph_block(self._graphs, n)
		sys.stdout.write("Generated %d graphs.\n" % len(self._graphs))
	
		sys.stdout.write("Generating types and flags...\n")
		self._flag_products = {}
		self._types = []
		self._flags = []
		num_types = 0
	
		for s in range(n % 2, n - 1, 2):
			
			these_types = generate_graphs(s, forbidden_edge_numbers=self.forbidden_edge_numbers)
			sys.stdout.write("Generated %d types of order %d, " % (
				len(these_types), s))

			m = (n + s) / 2
			these_flags = []
			for tg in these_types:
				these_flags.append(generate_flags(m, tg, forbidden_edge_numbers=self.forbidden_edge_numbers))
			sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))
			
			#for gi in range(len(self._graphs)):
			#	flag_products = multiple_equal_flag_products(self._graphs[gi],
			#			these_types, s, these_flags, m)
			#	for mi in range(len(these_types)):
			#		self._flag_products[(gi, mi + num_types)] = flag_products[mi]
			
			self._types.extend(these_types)
			self._flags.extend(these_flags)
			num_types += len(these_types)

 	def calculate_flag_products(self):
 	
 		sys.stdout.write("Averaging flag products...\n")
 		for ti in range(len(self._types)):
 			tg = self._types[ti]
 			s = tg[0]
 			m = (self._n + s) / 2
 			print ti, m, s
			self._flag_products[ti] = gb_flag_products(self._graph_block, tg, self._flags[ti], m, self._flags[ti], m)
			
# 			for gi in range(len(self._graphs)):
# 				self._flag_products[(gi, ti)] = multiple_equal_flag_products(self._graphs[gi],
# 						[tg], s, [self._flags[ti]], m)


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
	
