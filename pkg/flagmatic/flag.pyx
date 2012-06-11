from sage.structure.sage_object cimport SageObject

cdef class Flag (SageObject):

	def _require_mutable(self):
		"""
		EXAMPLES::

			sage: G = Flag()
			sage: G._require_mutable()
			sage: G.set_immutable()
			sage: G._require_mutable()
			Traceback (most recent call last):
			...
			ValueError: object is immutable; please change a copy instead.
		"""
		if self._is_immutable:
			raise ValueError("object is immutable; please change a copy instead.")

	def set_immutable(self):
		"""
		Make this object immutable, so it can never again be changed.

		EXAMPLES::
		
			sage: G = GraphFlag(2)
			sage: G.add_edge((1, 2))
			sage: G
			2:12
			sage: G.set_immutable()
			sage: G.delete_edge((1, 2))
			Traceback (most recent call last):
			...
			ValueError: object is immutable; please change a copy instead.
		"""
		self._is_immutable = True

	def is_immutable(self):
		"""
		Return True if this object is immutable (can not be changed)
		and False if it is not.

		To make this object immutable use :meth:`set_immutable`.

		EXAMPLE::
		
			sage: G=Flag()
			sage: G.is_immutable()
			False
			sage: G.set_immutable()
			sage: G.is_immutable()
			True
		"""
		try:
			return self._is_immutable
		except AttributeError:
			return False

	def is_mutable(self):
		"""
		EXAMPLES::
		
			sage: G=Flag()
			sage: G.is_mutable()
			True
			sage: G.set_immutable()
			sage: G.is_mutable()
			False
		"""
		try:
			return not self._is_immutable
		except AttributeError:
			return True
