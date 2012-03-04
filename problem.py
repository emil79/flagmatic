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

import os
import pexpect
import sys

cdsp_binary = "/Users/emil/Projects/flagmatic/csdp"

class flagmatic_problem(object):

	forbidden_edge_numbers = {}
	forbidden_graphs = []
	forbidden_induced_graphs = []

	_flag_products = []

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
	
		self._graph_densities = []
		for g in self._graphs:
			self._graph_densities.append(len(g[1]) / binomial(n, 3))
	
		sys.stdout.write("Generating types and flags...\n")
		self._flag_products = []
		self._types = []
		self._flags = []
	
		for s in range(n % 2, n - 1, 2):
			
			these_types = generate_graphs(s, forbidden_edge_numbers=self.forbidden_edge_numbers)
			sys.stdout.write("Generated %d types of order %d, " % (
				len(these_types), s))

			m = (n + s) / 2
			these_flags = []
			for tg in these_types:
				these_flags.append(generate_flags(m, tg, forbidden_edge_numbers=self.forbidden_edge_numbers))
			sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))
						
			self._types.extend(these_types)
			self._flags.extend(these_flags)


 	def calculate_flag_products(self):
 	
		self._flag_products = []
 		sys.stdout.write("Averaging flag products...\n")
 		for ti in range(len(self._types)):
 			tg = self._types[ti]
 			s = tg[0]
 			m = (self._n + s) / 2
 			sys.stdout.write("Doing type %d (order %d; flags %d).\n" % (ti + 1, s, m))
 			flags_block = make_graph_block(self._flags[ti], m)
			self._flag_products.append(flag_products(self._graph_block,
				make_graph_block([tg], s), flags_block, None))

	@property
	def graphs(self):
		return self._graphs

	@property
	def types(self):
		return self._types

	@property
	def flags(self):
		return self._flags
	
	def write_sdp_input_file(self):
	
		num_graphs = len(self._graphs)
		num_types = len(self._types)
	
		self._sdp_input_filename = os.path.join(SAGE_TMP, "sdp.dat-s")
	
		with open(self._sdp_input_filename, "w") as f:
	
			f.write("%d\n" % (num_graphs + 1,))
			f.write("%d\n" % (num_types + 3,))
			f.write("1 %s -%d -1\n" % (" ".join([str(len(P._flags[i])) for i in range(num_types)]),
				num_graphs))
			f.write("%s1.0\n" % ("0.0 " * num_graphs,))
			f.write("0 1 1 1 -1.0\n")
	
			for i in range(num_graphs):
				f.write("%d 1 1 1 -1.0\n" % (i + 1,))
				f.write("%d %d %d %d 1.0\n" % (i + 1, num_types + 2, i + 1, i + 1))
	
			for i in range(num_graphs):
				d = self._graph_densities[i]
				if d != 0:
					f.write("%d %d 1 1 %s\n" % (i + 1, num_types + 3, d.n(digits=64)))
			
			f.write("%d %d 1 1 1.0\n" % (num_graphs + 1, num_types + 3))
		
			for i in range(num_graphs):
				for j in range(num_types):
					for key, value in self._flag_products[j][i].dict().iteritems():
						f.write("%d %d %d %d %s\n" % (i + 1, j + 2, key[0] + 1, key[1] + 1,
							value.n(digits=64)))

	def run_csdp(self):
	
		self._sdp_output_filename = os.path.join(SAGE_TMP, "sdp.out")
	
		cmd = "%s %s %s" % (cdsp_binary, self._sdp_input_filename, self._sdp_output_filename)
		print cmd
		p = pexpect.spawn(cmd, timeout=60*60*24*7)
		
		while True:
			if p.eof():
				break
			try:
				p.expect("\r\n")
				line = p.before.strip() + "\n"
				sys.stdout.write(line)
			except pexpect.EOF:
				break
		p.close()
		returncode = p.exitstatus
		