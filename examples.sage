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

		problem = None
		construction = None

		base_filename = os.path.join("examples", name)	
		with open(base_filename + ".sage", "r") as f:
			lines = f.readlines()
		
		# last line of file might not end in a newline
		if lines[-1][-1] != "\n":
			lines[-1] += "\n"
		lines.append('problem.write_certificate("' + base_filename + '.js")\n')
		lines.append('problem.save("' + base_filename + '")\n')

		buffer = cStringIO.StringIO()
		saved_stdout = sys.stdout
		sys.stdout = Tee(sys.stdout, buffer)
		
		for line in lines:
			sys.stdout.write("sage: %s" % line)
			try:
				exec(preparse(line))
			except:
				traceback.print_exc(None, sys.stdout)
				break
		
		sys.stdout = saved_stdout
	
		self.problem = problem
		self.construction = construction
		self.output = buffer.getvalue()
	
		with open(base_filename + ".txt", "w") as f:
			f.write(self.output)
