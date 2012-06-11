from sage.structure.sage_object cimport SageObject

cdef class Flag (SageObject):

	cdef public bint _is_immutable
