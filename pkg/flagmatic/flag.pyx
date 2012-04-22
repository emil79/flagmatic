from sage.structure.sage_object import SageObject       
from sage.structure.sage_object cimport SageObject

cdef class Flag (SageObject):

	def hello(self):
	
		print "Hello"