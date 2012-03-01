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

class flagmatic_problem(object):

	forbidden_edge_numbers = {}
	forbidden_graphs = []
	forbidden_induced_graphs = []

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

		print "Generating graphs..."
		self._graphs = generate_graphs(n, forbidden_edge_numbers=self.forbidden_edge_numbers)
		print "Generated %d graphs." % len(self._graphs)
	
		print "Generating types and flags..."
		for s in range(n % 2, n - 1, 2):
			m = (n + s) / 2
			these_types = generate_graphs(s, forbidden_edge_numbers=self.forbidden_edge_numbers)
			these_flags = []
			for tg in these_types:
				these_flags.append(generate_flags(m, tg, forbidden_edge_numbers=self.forbidden_edge_numbers))
			for g in self._graphs:
				these_flags_products = flag_products(g, s, m, these_types, these_flags)
			self._types.extend(these_types)
			self._flags.extend(these_flags)

		print "Generating flags..."
		for tg in self._types:
			self._flags.append(generate_flags((n + tg[0]) / 2, tg, forbidden_edge_numbers=self.forbidden_edge_numbers))
	
	@property
	def graphs(self):
		return self._graphs

	@property
	def types(self):
		return self._types

	@property
	def flags(self):
		return self._flags
	

