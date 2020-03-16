# -*- coding: utf-8 -*-
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

import gzip, json, os, sys
import numpy
import pexpect

from sage.structure.sage_object import SageObject
from sage.rings.all import Integer, QQ, ZZ, RDF
from sage.functions.other import floor
from sage.matrix.all import matrix, identity_matrix, block_matrix, block_diagonal_matrix
from sage.modules.misc import gram_schmidt
from sage.misc.misc import SAGE_TMP 
from copy import copy

from hypergraph_flag import make_graph_block
from flag import *
from three_graph_flag import *
from graph_flag import *
from oriented_graph_flag import *
from multigraph_flag import *
from construction import *

# pexpect in Sage 4.8 has a bug, which prevents it using commands with full paths.
# So for now, CSDP has to be in a directory in $PATH.

cdsp_cmd = "csdp"
sdpa_cmd = "sdpa"
sdpa_dd_cmd = "sdpa_dd"
sdpa_qd_cmd = "sdpa_qd"
dsdp_cmd = "dsdp"


def block_structure(M):
    r"""
    Given a matrix, this function returns a tuple. The first entry is the number of
    row subdivisions that the matrix has. The second entry is a list of the sizes of the
    row subdivisions, and the third entry is a list containing the rows at which each
    subdivision begins. (Note that column subdivisions are ignored.)
    """
    row_div = M.subdivisions()[0]
    div_offsets = [0] + row_div
    div_sizes = row_div + [M.nrows()]
    num_blocks = len(div_sizes)

    for i in range(1, num_blocks):
        div_sizes[i] -= div_sizes[i - 1]

    return num_blocks, div_sizes, div_offsets


def safe_gram_schmidt(M):
    r"""
    Performs Gram Schmidt orthogonalization using Sage functions. The returned matrix
    will have the same subdivisions are the input matrix, and it will be sparse if and
    only if the input matrix is sparse.

    In addition, the new matrix.gram_schmidt method appears to have certain issues, at
    least under Sage 4.8, when dealing with number fields (in particular the ones used
    by the 'maxs3' and 'maxs4' problems.

    So, this function uses the old gram_schmidt function, whenever the base ring of the
    matrix is not the rational field.
    """

    subdivision = M.subdivisions()
    BR = M.base_ring()
    sparse = M.is_sparse()

    if BR == QQ:
        M, mu = M.gram_schmidt()
        # .gram_schmidt doesn't appear to preserve sparsity, so recreate matrix from rows.
        M = matrix(QQ, M.rows(), sparse=sparse)

    else:  # .gram_schmidt is broken for some number fields in 4.8.
        rows, mu = gram_schmidt(M.rows())
        M = matrix(BR, rows, sparse=sparse)

    M.subdivide(subdivision)
    return M


def LDLdecomposition(M):  # TODO: does this handle matrices with zero eigenvalues?
    MS = M.parent()
    D = MS.matrix()
    if M.is_zero():
        D.set_immutable()
        return D, D
    L = copy(MS.identity_matrix())
    for i in xrange(M.nrows()):
        for j in xrange(i):
            L[i, j] = (Integer(1) / D[j, j]) * (M[i, j] - sum(L[i, k] * L[j, k] * D[k, k] for k in xrange(j)))
        D[i, i] = M[i, i] - sum(L[i, k] ** 2 * D[k, k] for k in xrange(i))
    L.set_immutable()
    D.set_immutable()
    return L, D


class Problem(SageObject):
    r"""
    This is the principal class of flagmatic. Objects of this class represent Turán-type
    problems.
    """

    def __init__(self, flag_cls, order=None, forbid_induced=None, forbid=None,
                 forbid_homomorphic_images=False, density=None, minimize=False,
                 type_orders=None, types=None, max_flags=None, compute_products=True):
        r"""
        Creates a new Problem object. Generally it is not necessary to call this method
        directly, as Problem objects are more easily created using the helper functions:

        sage: problem = GraphProblem()
        sage: problem = ThreeGraphProblem()
        sage: problem = OrientedGraphProblem()
        sage: problem = TwoMultigraphProblem()
        sage: problem = ThreeMultigraphProblem()

        If this method is called directory, then a class that inherits from Flag should be
        provided, for example:

        sage: problem = Problem(GraphFlag)

        """

        self._flagmatic_version = "2.0"

        if issubclass(flag_cls, Flag):
            self._flag_cls = flag_cls
        else:
            raise ValueError

        self._n = 0

        self._field = QQ
        self._approximate_field = RDF

        self._forbidden_edge_numbers = []
        self._forbidden_graphs = []
        self._forbidden_induced_graphs = []

        self.state("specify", "yes")
        self.set_objective(minimize=minimize)

        if density is None:
            self.set_density(flag_cls.default_density_graph())
        else:
            self.set_density(density)

        if not forbid_induced is None:
            self.forbid_induced(forbid_induced)

        if not forbid is None:
            self.forbid(forbid)

        if forbid_homomorphic_images:
            self.forbid_homomorphic_images()

        if not order is None:
            self.generate_flags(order, type_orders=type_orders, types=types, max_flags=max_flags,
                                compute_products=compute_products)

    def state(self, state_name=None, action=None):
        r"""
        Keeps track of which things have been done. To get a list of all the states, enter

        sage: problem.state()

        If the argument state_name is supplied, then the value of the state called
        state_name is returned. This will be either "yes", "no" or "stale". "yes" means
        that the thing has been done, "no" means that it has not. "stale" means that it
        has been done, but subsequent actions mean that it needs to be re-done.

        If an action is supplied, it will change the value of the state state_name. The
        action must be one of "yes", "no" or "stale". It is not recommended that the
        value of any states be changed by the user.

        Setting ``force`` to True allows states to be set irrespective of other
        states.

        """

        # We use tuples here because there is an order to the states.
        state_list = [
            ("specify", {
                "requires": [],
                "depends": []
            }),
            ("set_objective", {
                "requires": [],
                "depends": []
            }),
            ("compute_flags", {
                "requires": [],
                "depends": ["specify"]
            }),
            ("set_construction", {
                "requires": ["compute_flags"],
                "depends": ["set_objective"]
            }),
            ("add_zero_eigenvectors", {
                "requires": ["set_construction"],
                "depends": []
            }),
            ("compute_block_bases", {
                "requires": ["compute_flags"],
                "depends": []
            }),
            ("compute_flag_bases", {
                "requires": ["set_construction"], # makes block bases if required
                "depends": ["add_zero_eigenvectors"]
            }),
            ("compute_products", {
                "requires": ["compute_flags"],
                "depends": []
            }),
            ("set_active_types", {
                "requires": ["compute_flags"],
                "depends": []
            }),
            ("set_block_matrix_structure", {
                "requires": ["compute_flags"],
                "depends": ["set_active_types"]
            }),
            ("write_sdp_input_file", {
                "requires": ["set_block_matrix_structure"],
                "depends": ["set_objective", "set_active_types"]
            }),
            ("write_sdp_initial_point_file", {
                "requires": ["set_block_matrix_structure"],
                "depends": ["set_objective", "set_active_types"]
            }),
            ("run_sdp_solver", {
                "requires": ["write_sdp_input_file"],
                "depends": ["write_sdp_input_file", "write_sdp_initial_point_file"]
            }),
            ("read_solution", {
                "requires": ["run_sdp_solver"],
                "depends": []
            }),
            ("check_solution", {
                "requires": ["read_solution"],
                "depends": []
            }),
            ("transform_solution", {
                "requires": ["compute_flag_bases", "check_solution"],
                "depends": []
            }),
            ("add_sharp_graphs", {
                "requires": ["set_construction"],
                "depends": []
            }),
            ("make_exact", {
                "requires": ["check_solution"],
                "depends": ["transform_solution", "add_sharp_graphs"]
            }),
            ("meet_target_bound", {
                "requires": ["transform_solution", "make_exact"],
                "depends": []
            }),
            ("check_exact", {
                "requires": ["make_exact"],
                "depends": ["meet_target_bound"]
            }),
            ("diagonalize", {
                "requires": ["meet_target_bound", "check_exact"],
                "depends": []
            }),
            ("write_certificate", {
                "requires": ["check_exact"],
                "depends": ["diagonalize"]
            })
        ]

        state_names = [s[0] for s in state_list]
        states = dict(state_list)

        if not hasattr(self, "_states"):
            self._states = dict((sn, "no") for sn in state_names)

        if state_name is None:
            width = max(len(sn) for sn in state_names)
            for sn in state_names:
                sys.stdout.write(("%%-%ds : %%s\n" % width) % (sn, self._states[sn]))
            return

        if not state_name in state_names:
            raise ValueError("unknown state.")

        if action == "yes" or action == "force_yes":
            if action == "yes" and not all(self._states[sn] == "yes" for sn in states[state_name]["requires"]):
                raise NotImplementedError("not ready for this yet!")
            self._states[state_name] = "yes"
            for sn in states:
                if state_name in states[sn]["depends"] and self._states[sn] == "yes":
                    self.state(sn, "stale")

        elif action == "stale":
            self._states[state_name] = "stale"
            for sn in states:
                if state_name in states[sn]["depends"]:
                    self.state(sn, "stale")

        elif action == "ensure_yes":
            if self._states[state_name] != "yes":
                raise NotImplementedError("not ready for this yet!")

        elif action == "ensure_yes_or_stale":
            if self._states[state_name] not in ["yes", "stale"]:
                raise NotImplementedError("not ready for this yet!")

        elif action is None:
            pass

        else:
            raise ValueError("unknown action.")

        return self._states[state_name]


    @property
    def flag_cls(self):
        return self._flag_cls


    @property
    def order(self):
        return self._n

    # Deprecated - use order instead.
    @property
    def n(self):
        return self._n

    # TODO: sanity checking of type orders

    def generate_flags(self, order, type_orders=None, types=None, max_flags=None, compute_products=True):
        r"""
        Generates the types and flags that will be used in the problem.

        INPUT:

         - ``order`` -- an integer. The order for the unlabelled flags (the admissible
           graphs).

         - ``type_orders`` -- (default: None) None, or a list of integers. If a list of
           integers is given, then types of these orders will be generated. If None is
           given then types of all orders less than ``order``, and congruent modulo 2 to
           ``order`` will be generated.

         - ``types`` - (default: None) A list of flags to use as types.

         - ``max_flags`` -- (default: None) None, or an integer. If an integer is given,
           then it provides an upper bound on the number of flags each type can have. If a
           type has more than this many flags, it will be removed.

         - ``compute_products`` -- (default: True) Boolean. If True then the flag products
           will be computed. For some large problems this may take a long time. If False,
           then the flag products must be computed later using the ``compute_products``
           method.
        """

        n = order
        if type_orders is None:
            orders = [(s, floor((n + s) / 2)) for s in range(n % 2, n - 1, 2)]
        else:
            orders = []
            for to in type_orders:
                if type(to) is tuple:
                    orders.append(to)
                else:
                    orders.append((to, floor((n + to) / 2)))

        self.state("compute_flags", "yes")
        self._n = n

        sys.stdout.write("Generating graphs...\n")

        import cProfile, pstats, StringIO
        pr = cProfile.Profile()
        pr.enable()

        self._graphs = self._flag_cls.generate_graphs(n, forbidden_edge_numbers=self._forbidden_edge_numbers,
                                                      forbidden_graphs=self._forbidden_graphs, forbidden_induced_graphs=self._forbidden_induced_graphs)

        pr.disable()
        ps = pstats.Stats(pr).sort_stats('cumulative')
        ps.print_stats()

        sys.stdout.write("x Generated %d graphs.\n" % len(self._graphs))

        for g in self._graphs:    # Make all the graphs immutable
            g.set_immutable()

        self._compute_densities()

        sys.stdout.write("Generating types and flags...\n")
        self._types = []
        self._flags = []

        allowed_types = []
        if types:
            for h in types:
                if isinstance(h, basestring):
                    h = self._flag_cls(h)
                if not isinstance(h, self._flag_cls):
                    raise ValueError
                allowed_types.append(h)

        for s, m in orders:

            these_types = self._flag_cls.generate_graphs(s, forbidden_edge_numbers=self._forbidden_edge_numbers,
                                                         forbidden_graphs=self._forbidden_graphs,
                                                         forbidden_induced_graphs=self._forbidden_induced_graphs)

            if types:
                these_types = [h for h in these_types if h in allowed_types]

            sys.stdout.write("Generated %d types of order %d, " % (len(these_types), s))

            these_flags = []
            for tg in these_types:
                these_flags.append(self._flag_cls.generate_flags(m, tg, forbidden_edge_numbers=self._forbidden_edge_numbers,
                                                                 forbidden_graphs=self._forbidden_graphs,
                                                                 forbidden_induced_graphs=self._forbidden_induced_graphs))
            sys.stdout.write("with %s flags of order %d.\n" % ([len(L) for L in these_flags], m))

            self._types.extend(these_types)
            self._flags.extend(these_flags)

        num_types = len(self._types)

        if not max_flags is None:
            bad_indices = [i for i in range(num_types) if len(self._flags[i]) > max_flags]
            if len(bad_indices) > 0:
                good_indices = [i for i in range(num_types) if not i in bad_indices]
                self._types = [self._types[i] for i in good_indices]
                self._flags = [self._flags[i] for i in good_indices]
                sys.stdout.write("Removed types %s as they have too many flags.\n" % bad_indices)

        num_types = len(self._types)  # may have changed!

        self._active_types = range(num_types)

        for ti in range(num_types):			  # Make everything immutable!
            self._types[ti].set_immutable()
            for g in self._flags[ti]:
                g.set_immutable()

        if compute_products:
            self.compute_products()


    @property
    def graphs(self):
        r"""
        Read-only. A (copy) of the list of admissible graphs. Modifying the list will have
        no effect on the Problem.

        """
        return copy(self._graphs)

    @property
    def types(self):
        """
        Read-only. A (copy) of the list of types. Modifying the list will have
        no effect on the Problem.

        """
        return copy(self._types)

    @property
    def flags(self):
        """
        Read-only. A (copy) of the flags, as a list of lists. The flags are listed by
        type, in the same order as the types appear in the types list. Modifying the lists
        will have no effect on the Problem.

        """
        return copy(self._flags)

    @property
    def density_graphs(self):
        """
        Read-only. A (copy) of the list of density graphs. Modifying the lists
        will have no effect on the Problem.

        """
        return copy(self._density_graphs)

    def set_objective(self, minimize=False):
        r"""
        Sets the Problem to be a "maximization" or a "minimization" problem.

        INPUT:

         - ``minimize`` - boolean (default: False) sets whether the problem is a
           minimization problem (or a maximization problem).

        In a maximization problem, the objective is to find a good (i.e. low) upper bound
        on a density. In a minimization problem the objective is to find a good (i.e.
        high) lower bound on a density.
        """
        if not type(minimize) is bool:
            raise ValueError

        self.state("set_objective", "yes")

        self._minimize = minimize

    def _compute_densities(self):

        self._densities = []
        for dg in self._density_graphs:
            density_values = []
            for g in self._graphs:
                dv = 0
                for h, coeff in dg:
                    if h.n == g.n:
                        # comparison will be fast, as both g and h should have
                        # _certified_minimal_isomorph set to True
                        if g == h:
                            dv += coeff
                    else:
                        dv += coeff * g.subgraph_density(h)
                density_values.append(dv)
            self._densities.append(density_values)

    def set_density(self, *args):

        self.state("set_objective", "yes")

        flattened_args = []
        for arg in args:
            if isinstance(arg, list):
                flattened_args.extend(arg)
            else:
                flattened_args.append(arg)

        density_graphs = []

        for arg in flattened_args:

            if isinstance(arg, basestring) and "." in arg:
                arg = tuple(map(int, arg.split(".")))

            if isinstance(arg, tuple):

                if len(arg) != 2:
                    raise ValueError

                if arg[0] in ZZ:

                    k, ne = arg
                    if k < self._flag_cls().r:
                        raise ValueError
                    max_e = self._flag_cls.max_number_edges(k)
                    if not ne in range(max_e + 1):
                        raise ValueError

                    # Don't forbid anything - if we do, we'll have to keep list
                    # updated whenever forbidden things change. Also it seems the most
                    # appropriate thing to do.
                    graphs = self._flag_cls.generate_graphs(k)
                    for g in graphs:
                        if g.ne == ne:
                            density_graphs.append((g, Integer(1)))
                    continue

                else:
                    h, coeff = arg

            else:
                h, coeff = arg, Integer(1)

            if isinstance(h, basestring):
                h = self._flag_cls(h)

            if not isinstance(h, self._flag_cls):
                raise ValueError

            h = copy(h)
            h.make_minimal_isomorph()
            density_graphs.append((h, coeff))

        if len(density_graphs) == 0:
            raise ValueError

        # Note that this function only sets one of the densities.
        self._density_graphs = [density_graphs]
        self._active_densities = [0]
        self._density_coeff_blocks = [[0]]

        if self.state("compute_flags") == "yes":
            self._compute_densities()

    def _forbid(self, h, induced):

        self.state("specify", "yes")

        if isinstance(h, basestring) and "." in h:
            h = tuple(map(int, h.split(".")))

        if isinstance(h, tuple):
            k, ne = h
            if k < self._flag_cls().r:
                raise ValueError
            max_e = self._flag_cls.max_number_edges(k)
            if not ne in range(max_e + 1):
                raise ValueError
            if induced:
                self._forbidden_edge_numbers.append((k, ne))
                sys.stdout.write("Forbidding %d-sets from spanning exactly %d edges.\n" % (k, ne))
            else:
                for i in range(ne, max_e + 1):
                    self._forbidden_edge_numbers.append((k, i))
                sys.stdout.write("Forbidding %d-sets from spanning at least %d edges.\n" % (k, ne))
            return

        if isinstance(h, basestring):
            h = self._flag_cls(h)

        if not isinstance(h, self._flag_cls):
            raise ValueError

        if induced:
            self._forbidden_induced_graphs.append(copy(h))
            self._forbidden_induced_graphs.sort(key=lambda g: (g.n, g.ne))
            sys.stdout.write("Forbidding %s as an induced subgraph.\n" % h)
        else:
            self._forbidden_graphs.append(copy(h))
            self._forbidden_graphs.sort(key=lambda g: (g.n, g.ne))
            sys.stdout.write("Forbidding %s as a subgraph.\n" % h)

    def forbid(self, *args):
        r"""
        Sets the problem to be in the theory of graphs that do not have contain the given subgraphs.

        INPUT:

         - arguments must be Flags, strings, or lists of these things. Flags must be of the
           appropriate class for the problem. Strings will be interpreted as representing flags
           of the appropriate class.

        """
        for h in args:
            if isinstance(h, list):
                for x in h:
                    self._forbid(x, False)
            else:
                self._forbid(h, False)


    def forbid_induced(self, *args):
        r"""
        Sets the problem to be in the theory of graphs that do not have contain the given
        graphs as induced subgraphs.

        INPUT:

         - arguments must be Flags, strings, or lists of these things. Flags must be of the
           appropriate class for the problem. Strings will be interpreted as representing flags
           of the appropriate class.

        """
        for h in args:
            if isinstance(h, list):
                for x in h:
                    self._forbid(x, True)
            else:
                self._forbid(h, True)


    def forbid_homomorphic_images(self):
        r"""
        Restricts the problem to be in the theory of graphs that do not contain homomorphic images
        of graphs already specified using ``forbid``. For certain problems this will make the
        computation simpler, without affecting the result.
        """
        L = sum([g.homomorphic_images() for g in self._forbidden_graphs], [])
        LM = self._flag_cls.minimal_by_inclusion(L)
        if len(LM) == 0:
            return
        #sys.stdout.write("Forbidding")
        for g in LM:
            #sys.stdout.write(" %s" % repr(g))
            self._forbid(g, False)
        #sys.stdout.write("\n")


    # TODO: warn if already solved

    def set_inactive_types(self, *args):
        r"""
        Specifies that the Q matrices for certain types should be zero matrices.

        INPUT:

        - arguments should be integers, specifying the indices of types in ``types``
          that should be marked as being "inactive".
        """
        self.state("set_active_types", "yes")

        num_types = len(self._types)

        for arg in args:
            ti = int(arg)
            if not ti in range(num_types):
                raise ValueError
            if ti in self._active_types:
                self._active_types.remove(ti)
            else:
                sys.stdout.write("Warning: type %d is already inactive.\n" % ti)


    def set_approximate_field(self, field):
        r"""
        Specifies the approximate field used when reading in and transforming the solution
        matrices.

        INPUT:

        - ``field`` - a field that uses finite precision, or in other words, a field that uses
          floating point arithmetic.

        EXAMPLES:

        sage: problem.set_approximate_field(RealField(113))

        This specifies that quadruple precision should be used. The default approximate field is
        the real dense field, RDF(), which uses 53 bits of precision.
        """
        if not field.is_field():
            raise ValueError("not a field.")

        if field.is_exact():
            raise ValueError("field must be floating point.")

        self._approximate_field = field


    # TODO: sanity check target_bound ?

    def set_extremal_construction(self, construction=None, field=None, target_bound=None):
        r"""
        Sets the extremal construction. This will be used to determine sharp graphs and forced
        zero eigenvectors.

        INPUT:

         - ``construction`` - a Construction object (default: None). None can be specified, in
           which case sharp graphs and forced zero eigenvectors can be added manually by using
           ``add_sharp_graphs`` and ``add_zero_eigenvectors``.

         - ``field`` - a field object (default: None). If ``construction`` is None, then this
           argument can be used to specify the exact field to use for computing the bound. This
           argument must be None if ``construction`` is not None, as the field will be taken from
           the Construction object.

         - ``target_bound`` - a number (default: None). If ``construction`` is None, then this
           argument can be used to specify the bound that should be aimed for. This argument must
           be None if ``construction`` is not None, as the target bound will be taken from
           the Construction object.
        """
        num_types = len(self._types)

        if construction is None:

            if not field.is_field():
                raise ValueError("not a valid field.")

            if not field.is_exact():
                raise ValueError("field must be an exact field (not floating point).")

            self.state("set_construction", "yes")

            self._construction = None
            self._field = field
            sys.stdout.write("Field is \"%s\" with embedding x=%s.\n" %
                (str(self._field), self._field.gen_embedding().n(digits=10)))
            self._target_bound = target_bound
            sys.stdout.write("Set target bound to be %s (%s).\n" %
                (self._target_bound, self._target_bound.n(digits=10)))

            self._zero_eigenvectors = []
            for ti in range(num_types):
                M = matrix(self._field, 0, len(self._flags[ti]), sparse=True)
                M.set_immutable()
                self._zero_eigenvectors.append(M)

            self._sharp_graphs = []
            self._sharp_graph_densities = []

            return

        if not isinstance(construction, Construction):
            raise ValueError("not a valid construction.")

        if not field is None:
            raise ValueError("field should be None if construction is given.")

        if not target_bound is None:
            raise ValueError("target_bound should be None if construction is given.")

        self.state("set_construction", "yes")

        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        self._construction = construction
        self._field = construction.field

        sys.stdout.write("Determining which graphs appear in construction...\n")

        sharp_graphs = construction.subgraph_densities(self._n)
        target_densities = [0 for j in range(num_densities)]
        self._sharp_graphs = []
        self._sharp_graph_densities = []

        for pair in sharp_graphs:
            g, den = pair
            for gi in range(num_graphs):
                if g.is_labelled_isomorphic(self._graphs[gi]):
                    self._sharp_graphs.append(gi)
                    self._sharp_graph_densities.append(den)
                    for j in range(num_densities):
                        target_densities[j] += self._densities[j][gi] * den
                    break
            else:
                sys.stdout.write("Warning: non-admissible graph %s appears in construction!\n" % g)

        # set target_bound to equal the maximum - probably this will always be what is wanted...
        self._target_bound = max(target_densities)

        sys.stdout.write("Density of construction is %s.\n" % self._target_bound)

        self._zero_eigenvectors = []

        for ti in range(len(self._types)):

            self._zero_eigenvectors.append(construction.zero_eigenvectors(self._types[ti], self._flags[ti]))

            sys.stdout.write("Found %d zero eigenvectors for type %d.\n" % (
                self._zero_eigenvectors[ti].nrows(), ti))

        for ti in range(len(self._types)):
            self._zero_eigenvectors[ti].set_immutable()


    # TODO: reinstate the per-block option?

    def add_zero_eigenvectors(self, ti, eigenvectors, use_bases=False):
        r"""
        Adds a zero eigenvector. This method is necessary when not all the zero eigenvectors
        can be determined from the construction.

        INPUT:

         - ``ti`` - integer specifying which type the eigenvector is for.

         - ``eigenvectors`` - a vector, or matrix, of zero eigenvector(s) to add for type
           ``ti``. If adding more than one vector, the vectors can be given as the rows of a
           matrix. (Alternatively, they can be added one at a time with multiple calls to
           this method.)

         - ``use_bases`` - specifies that the vector is given in the basis of the Q' matrix, as
           opposed to the standard basis. The vector will be tranformed before being added to the
           zero eigenvectors.
        """
        self.state("add_zero_eigenvectors", "yes")

        if use_bases:
            self.state("compute_flag_bases", "ensure_yes_or_stale")
            NZ = (self._flag_bases[ti].T).solve_left(eigenvectors)
        else:
            NZ = eigenvectors

        self._zero_eigenvectors[ti] = self._zero_eigenvectors[ti].stack(NZ)
        self._zero_eigenvectors[ti].set_immutable()

    # TODO: is good idea to assume zero densities?

    def add_sharp_graphs(self, *args):
        r"""
        Adds a sharp graph. This method is necessary when not all the sharp graphs can be
        determined from the construction.

        INPUT:

         - arguments are integers specifying the indices of the graphs to set as sharp.
        """
        self.state("add_sharp_graphs", "yes")

        num_graphs = len(self._graphs)

        for arg in args:
            si = int(arg)
            if not si in range(num_graphs):
                raise ValueError
            if not si in self._sharp_graphs:
                self._sharp_graphs.append(si)
                self._sharp_graph_densities.append(Integer(0))
            else:
                sys.stdout.write("Warning: graph %d is already marked as sharp.\n" % si)

    def change_solution_bases(self, use_blocks=True):
        r"""
        Transforms the solution's Q matrices, so that they are (hopefully) positive definite. A
        construction should have been set previously, and this will be used to determine forced
        zero eigenvectors.

        This method is called from ``make_exact`` by default, and so is not normally explicitly
        invoked.

        INPUT:

         - ``use_blocks`` - Boolean (default: True). Specifies whether to apply an additional
           change of basis so that the matrices have a block structure with two blocks. This uses
           the invariant anti-invariant idea of Razborov.
        """

        if self.state("compute_flag_bases") != "yes":
            self.compute_flag_bases(use_blocks)

        self.state("transform_solution", "yes")

        num_types = len(self._types)

        sys.stdout.write("Transforming matrices")

        self._sdp_Qdash_matrices = []

        for ti in range(num_types):

            B = self._flag_bases[ti]
            if B.nrows() > 0:
                row_div = B.subdivisions()[0]
                M = B * self._sdp_Q_matrices[ti] * B.T
                M.subdivide(row_div, row_div)
                # zero out bits that should be zero. Note the copy() seems to be needed.
                M = block_diagonal_matrix([copy(M.subdivision(i,i)) for i in range(len(row_div) + 1)])
                M.set_immutable()
                self._sdp_Qdash_matrices.append(M)
            else:
                self._sdp_Qdash_matrices.append(matrix(self._approximate_field, 0, 0))
            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

    def compute_block_bases(self):
        r"""
        Computes a basis for the solution's Q matrices, that will give them a block structure with
        two blocks. This uses the invariant anti-invariant idea of Razborov. This method is not
        normally explicitly invoked. By default, ``change_problem_bases`` and
        ``change_solution_bases`` call it.
        """
        self.state("compute_block_bases", "yes")

        self._block_bases = []

        for ti in range(len(self._types)):

            B = self._flag_cls.flag_basis(self._types[ti], self._flags[ti])
            num_blocks, block_sizes, block_offsets = block_structure(B)
            sys.stdout.write("Type %d (%d flags) blocks: %s \n" % (ti, len(self._flags[ti]), block_sizes))
            self._block_bases.append(B)

    def compute_flag_bases(self, use_blocks=True, keep_rows=False, use_smaller=False):
        r"""
        Computes a basis for the solution's Q matrices, using the construction to determine forced
        zero eigenvectors. This method is used by ``change_problem_bases`` and
        ``change_solution_bases``, and would not usually be invoked directly.
        """
        self.state("compute_flag_bases", "yes")

        num_types = len(self._types)

        if use_blocks and self.state("compute_block_bases") != "yes":
            self.compute_block_bases()

        self._flag_bases = []

        sys.stdout.write("Creating bases")
        sys.stdout.flush()

        for ti in range(num_types):

            if use_blocks:
                num_blocks, block_sizes, block_offsets = block_structure(self._block_bases[ti])
            else:
                num_blocks, block_sizes, block_offsets = 1, [len(self._flags[ti])], [0]

            BS = []

            for bi in range(num_blocks):

                Z = self._zero_eigenvectors[ti]

                if use_blocks:
                    B = (self._block_bases[ti].subdivision(bi, 0) * Z.T).T
                else:
                    B = Z

                B = B.echelon_form()

                nzev = B.rank()
                B = B[:nzev, :]

                if nzev == 0:
                    B = identity_matrix(QQ, block_sizes[bi], sparse=True)

                elif nzev == block_sizes[bi]:
                    pass

                else:
                    B = B.stack(B.right_kernel().basis_matrix())

                if use_blocks:
                    B = B * self._block_bases[ti].subdivision(bi, 0)

                if not keep_rows:
                    B = B[nzev:, :]  # delete rows corresponding to zero eigenvectors

                if B.nrows() > 0:
                    BS.append(B)

            M = block_matrix([[B] for B in BS], subdivide=True)

            if M.nrows() == 0:
                M = matrix(self._field, 0, len(self._flags[ti]), sparse=True)
            else:
                M = safe_gram_schmidt(M)

            M.set_immutable()
            self._flag_bases.append(M)

            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

        self._inverse_flag_bases = []

        for ti in range(num_types):

            M = copy(self._flag_bases[ti])

            if use_smaller:
                MT = M.T
                MT.set_immutable()
                self._inverse_flag_bases.append(MT)
                for j in range(M.nrows()):
                    M[j, :] /= sum([x ** 2 for x in M.row(j)])
                M.set_immutable()
                self._flag_bases[ti] = M

            else:
                for j in range(M.nrows()):
                    M[j, :] /= sum([x ** 2 for x in M.row(j)])
                MT = M.T
                MT.set_immutable()
                self._inverse_flag_bases.append(MT)


    def compute_products(self):
        r"""
        Computes the products of the flags. This method is by default called from
        ``generate_flags``, and so would normally not need to be invoked directly.
        """
        self.state("compute_products", "yes")

        num_types = len(self._types)
        graph_block = make_graph_block(self._graphs, self._n)
        self._product_densities_arrays = []

        sys.stdout.write("Computing products")

        for ti in range(num_types):

            tg = self._types[ti]
            s = tg.n
            m = (self._n + s) / 2

            flags_block = make_graph_block(self._flags[ti], m)
            rarray = self._flag_cls.flag_products(graph_block, tg, flags_block, None)
            self._product_densities_arrays.append(rarray)

            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

    def _set_block_matrix_structure(self):

        self.state("set_block_matrix_structure", "yes")

        self._block_matrix_structure = []

        for ti in self._active_types:

            num_blocks, block_sizes, block_offsets = 1, [len(self._flags[ti])], [0]

            # Remove zero-sized blocks
            bi = 0
            while bi < num_blocks:
                if block_sizes[bi] == 0:
                    num_blocks -= 1
                    del block_sizes[bi]
                    del block_offsets[bi]
                else:
                    bi += 1

            for bi in range(num_blocks):
                self._block_matrix_structure.append((ti, block_sizes[bi], block_offsets[bi]))

    def _get_block_matrix_structure(self, ti):

        num_blocks = 0
        block_indices = []
        block_offsets = []
        block_sizes = []

        for bi in range(len(self._block_matrix_structure)):
            b = self._block_matrix_structure[bi]
            if b[0] == ti:
                num_blocks += 1
                block_indices.append(bi)
                block_sizes.append(b[1])
                block_offsets.append(b[2])

        return num_blocks, block_sizes, block_offsets, block_indices

    def solve_sdp(self, show_output=False, solver="csdp",
        force_sharp_graphs=False, force_zero_eigenvectors=False,
        check_solution=True, tolerance=1e-5, show_sorted=False, show_all=False,
        use_initial_point=False, import_solution_file=None):
        r"""
        Solves a semi-definite program to get a bound on the problem.

        INPUT:

         - ``show_output`` - Boolean (default: False). Whether to display output from the SDP
           solver.

         - ``solver`` - String (default: "csdp"). The SDP solver command to use. This can be one
           of the following:

            - "csdp" (Default) : the CSDP solver.
            - "sdpa" : the SDPA solver.
            - "sdpa_dd" : the double precision variant of SDPA
            - "sdpa_qd" : the quadruple precision variant of SDPA
            - "dsdp" : the DSDP solver.

            Note that the SDP solver must be present on the system; and it should be in a
            directory listed in PATH. The name of the solver should be "csdp", "sdpa", "sdpa_dd",
            "sdpa_qd" or "dsdp".

         - ``force_sharp_graphs`` - Boolean (default: False). If True, then the SDP is set up so
           that graphs that are supposed to be sharp are not given any "slack". Generally, this
           option is not particularly useful. It can sometimes improve the "quality" of a solution.

         - ``check_solution`` - Boolean (default: True). Whether to run ``check_solution`` to see
           if the bound appears to be tight; i.e. whether it is sufficiently close to the density
           given by the extremal construction.

         - ``tolerance`` - Number (default: 0.00001). This argument is passed to
           ``check_solution``. If a graph has a coefficient whose absolute difference with the
           bound is less than ``tolerance``, then it is considered to be sharp.

         - ``show_sorted`` - Boolean (default: False). This argument is passed to
           ``check_solution``. Whether to sort the sharp graphs according to their coefficients.
           If False, then they will be displayed in order of the graph indices.

          - ``show_all`` - Boolean (default: False). This argument is passed to
            ``check_solution``. If True, then the coefficients of all the graphs will be displayed
            instead of just the sharp graphs. In this case, the graphs that appear to be sharp are
            annotated with an "S" and the graphs that are forced to be sharp by the construction
            are annotated with a "C". (If these sets are not identical, then there is a problem.)

          - ``use_initial_point`` - Boolean (default: False). Whether to write an initial point
            file for the SDP solver. The initial point file is not used unless the solver is CSDP.
            Using an initial point can speed up the computation, but occasionally causes problems;
            so this option is False by default.

          - ``import_solution_file`` - Filename or None (default: None). If not None, then the SDP
            solver will not be run; instead the output file from a previous run of an SDP solver
            will be read. Care should be taken to ensure that the file being imported is for
            exactly the same problem, as minimal sanity-checking is done.
        """

        if import_solution_file is None:

            if self.state("write_sdp_input_file") != "yes":
                self.write_sdp_input_file(force_sharp_graphs=force_sharp_graphs,
                                          force_zero_eigenvectors=force_zero_eigenvectors)
            if use_initial_point and self.state("write_sdp_initial_point_file") != "yes":
                self.write_sdp_initial_point_file()
            self._run_sdp_solver(show_output=show_output, solver=solver,
                                 use_initial_point=use_initial_point)

        else:

            self._sdp_output_filename = import_solution_file
            # pretend we have run the solver!
            self.state("run_sdp_solver", "force_yes")

        self._read_sdp_output_file()

        if check_solution:
            self.check_solution(tolerance=tolerance, show_sorted=show_sorted, show_all=show_all)

    # TODO: add option for forcing sharps

    def write_sdp_input_file(self, force_sharp_graphs=False, force_zero_eigenvectors=False):
        r"""
        Writes an input file for the SDP solver, specifying the SDP to be solved. This method is
        by default called by ``solve_sdp``.

        INPUT:

         - ``force_sharp_graphs`` - Boolean (default: False). If True, then the SDP is set up so
           that graphs that are supposed to be sharp are not given any "slack". Generally, this
           option is not particularly useful. It can sometimes improve the "quality" of a solution.
        """
        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_active_densities = len(self._active_densities)
        num_density_coeff_blocks = len(self._density_coeff_blocks)

        if num_active_densities < 1:
            raise NotImplementedError("there must be at least one active density.")

        if num_density_coeff_blocks < 1:
            raise NotImplementedError("there must be at least one density coefficient block.")

        if self.state("set_block_matrix_structure") != "yes":
            self._set_block_matrix_structure()
        total_num_blocks = len(self._block_matrix_structure)

        if force_zero_eigenvectors:
            num_extra_matrices = sum(self._zero_eigenvectors[ti].nrows() for ti in self._active_types)
        else:
            num_extra_matrices = 0

        self.state("write_sdp_input_file", "yes")

        self._sdp_input_filename = os.path.join(unicode(SAGE_TMP), "sdp.dat-s")

        sys.stdout.write("Writing SDP input file...\n")

        with open(self._sdp_input_filename, "w") as f:

            f.write("%d\n" % (num_graphs + num_density_coeff_blocks + num_extra_matrices,))
            f.write("%d\n" % (total_num_blocks + 3 + (1 if force_zero_eigenvectors else 0),))

            f.write("1 ")
            for b in self._block_matrix_structure:
                f.write("%d " % b[1])

            f.write("-%d -%d" % (num_graphs, num_active_densities))
            if force_zero_eigenvectors:
                f.write(" -%d" % num_extra_matrices)
            f.write("\n")

            f.write("0.0 " * num_graphs)
            f.write("1.0 " * num_density_coeff_blocks)
            f.write("0.0 " * num_extra_matrices)
            f.write("\n")

            if not self._minimize:
                f.write("0 1 1 1 -1.0\n")
            else:
                f.write("0 1 1 1 1.0\n")

            if force_zero_eigenvectors:
                for mi in range(num_extra_matrices):
                    f.write("0 %d %d %d %s\n" % (total_num_blocks + 4, mi + 1, mi + 1, "1.0" if self._minimize else "-1.0"))

            for i in range(num_graphs):
                if not self._minimize:
                    f.write("%d 1 1 1 -1.0\n" % (i + 1,))
                else:
                    f.write("%d 1 1 1 1.0\n" % (i + 1,))
                if not (force_sharp_graphs and i in self._sharp_graphs):
                    f.write("%d %d %d %d 1.0\n" % (i + 1, total_num_blocks + 2, i + 1, i + 1))

            for i in range(num_graphs):
                for j in range(num_active_densities):
                    d = self._densities[self._active_densities[j]][i]
                    if d != 0:
                        if self._minimize:
                            d *= -1
                        f.write("%d %d %d %d %s\n" % (i + 1, total_num_blocks + 3, j + 1, j + 1, d.n(digits=64)))

            for i in range(num_density_coeff_blocks):
                for di in self._density_coeff_blocks[i]:
                    if di in self._active_densities:
                        j = self._active_densities.index(di)
                        f.write("%d %d %d %d 1.0\n" % (num_graphs + i + 1, total_num_blocks + 3, j + 1, j + 1))

            for ti in self._active_types:

                num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                for row in self._product_densities_arrays[ti]:
                    gi = row[0]
                    j = row[1]
                    k = row[2]
                    bi = num_blocks - 1
                    if bi > 0:
                        while block_offsets[bi] > j:
                            bi -= 1
                        j -= block_offsets[bi]
                        k -= block_offsets[bi]
                    value = Integer(row[3]) / Integer(row[4])
                    f.write("%d %d %d %d %s\n" %
                            (gi + 1, block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

            # TODO: get working with blocks, inactive types
            if force_zero_eigenvectors:
                mi = 0
                for ti in self._active_types:
                    nf = len(self._flags[ti])
                    for zi in range(self._zero_eigenvectors[ti].nrows()):
                        for j in range(nf):
                            for k in range(j, nf):
                                value = self._zero_eigenvectors[ti][zi, j] * self._zero_eigenvectors[ti][zi, k]
                                if value != 0:
                                    f.write("%d %d %d %d %s\n" %
                                            (num_graphs + num_density_coeff_blocks + mi + 1, ti + 2, j + 1, k + 1, value.n(digits=64)))
                        f.write("%d %d %d %d -1.0\n" % (num_graphs + num_density_coeff_blocks + mi + 1, total_num_blocks + 4, mi + 1, mi + 1))
                        mi += 1

    # TODO: handle no sharp graphs

    def write_sdp_initial_point_file(self, small_change=1/Integer(10)):
        r"""
        Writes an initial point file for the SDP solver. The zero matrix gives a feasible primal
        point. If a construction has been set, then this will be used to generate a feasible dual
        point, otherwise a zero matrix will be used for this as well. Using this method can reduce
        the time it takes the SDP solver to find a solution - typically by a third.

        INPUT:

         - ``small_change`` - Number (default: 1/10). Both primal and dual points must be
           perturbed by adding a small positive amount to the leading diagonals, in order that the
           matrices have no zero eigenvalues. Smaller values are not necessarily better here. The
           optimal value seems to depend on the problem and solver used.
        """

        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_active_densities = len(self._active_densities)

        self._sdp_initial_point_filename = os.path.join(unicode(SAGE_TMP), "sdp.ini-s")

        if self.state("set_block_matrix_structure") != "yes":
            self._set_block_matrix_structure()
        total_num_blocks = len(self._block_matrix_structure)

        self.state("write_sdp_initial_point_file", "yes")

        sys.stdout.write("Writing SDP initial point file...\n")

        with open(self._sdp_initial_point_filename, "w") as f:

            if self.state("set_construction") == "yes":

                for gi in range(num_graphs):
                    if gi in self._sharp_graphs:
                        si = self._sharp_graphs.index(gi)
                        f.write("%s " % self._sharp_graph_densities[si].n(digits=64))
                    else:
                        f.write("0.0 ")

                if not self._minimize:
                    f.write("%s\n" % (-self._target_bound).n(digits=64))
                else:
                    f.write("%s\n" % self._target_bound.n(digits=64))

                f.write("1 1 1 1 %s\n" % small_change.n(digits=64))

                for ti in range(num_types):

                    nf = len(self._flags[ti])
                    z_matrix = matrix(self._field, nf, nf)

                    for row in self._product_densities_arrays[ti]:
                        gi = row[0]
                        if not gi in self._sharp_graphs:
                            continue
                        si = self._sharp_graphs.index(gi)
                        j = row[1]
                        k = row[2]
                        value = Integer(row[3]) / Integer(row[4])
                        z_matrix[j, k] += value * self._sharp_graph_densities[si]

                    for j in range(nf):
                        z_matrix[j, j] += small_change

                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            for k in range(j, block_sizes[bi]):
                                value = z_matrix[block_offsets[bi] + j, block_offsets[bi] + k]
                                if value != 0:
                                    f.write("1 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

                for gi in range(num_graphs):
                    if gi in self._sharp_graphs:
                        si = self._sharp_graphs.index(gi)
                        value = self._sharp_graph_densities[si]
                    else:
                        value = 0
                    if value <= 0:
                        value = small_change
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

                for j in range(num_active_densities):
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, small_change.n(digits=64)))

            else:

                for gi in range(num_graphs + 1):
                    f.write("0.0 ")
                f.write("\n")
                f.write("1 1 1 1 %s\n" % small_change.n(digits=64))
                for ti in range(num_types):
                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)
                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            f.write("1 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, j + 1, small_change.n(digits=64)))
                for gi in range(num_graphs):
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, small_change.n(digits=64)))
                for j in range(num_active_densities):
                    f.write("1 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, small_change.n(digits=64)))

            # TODO: make this an exact Q check.
            if self.state("check_exact") == "yes":

                f.write("2 1 1 1 %s\n" % self._bound.n(digits=64))

                for ti in range(num_types):

                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            for k in range(j, block_sizes[bi]):
                                value = self._exact_Q_matrices[ti][block_offsets[bi] + j, block_offsets[bi] + k]
                                if j == k:
                                    value += small_change
                                if value != 0:
                                    f.write("2 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, k + 1, value.n(digits=64)))

                for gi in range(num_graphs):
                    value = self._bound - self._bounds[gi]
                    if value <= 0:
                        value = small_change
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

                for j in range(num_active_densities):
                    value = self._exact_density_coeffs[self._active_densities[j]]
                    if value <= 0:
                        value = small_change
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, value.n(digits=64)))

            else:

                densities = [sum([self._densities[j][gi] for j in self._active_densities])
                             / num_active_densities for gi in range(num_graphs)]

                bound = min(densities) if self._minimize else max(densities)

                value = bound
                if value <= 0:
                    value = small_change
                f.write("2 1 1 1 %s\n" % value.n(digits=64))

                num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                for ti in range(num_types):

                    num_blocks, block_sizes, block_offsets, block_indices = self._get_block_matrix_structure(ti)

                    for bi in range(num_blocks):
                        for j in range(block_sizes[bi]):
                            f.write("2 %d %d %d %s\n" % (block_indices[bi] + 2, j + 1, j + 1, small_change.n(digits=64)))

                for gi in range(num_graphs):
                    value = (bound - densities[gi]) * (-1 if self._minimize else 1)
                    if value <= 0:
                        value = small_change
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 2, gi + 1, gi + 1, value.n(digits=64)))

                value = Integer(1) / num_active_densities
                for j in range(num_active_densities):
                    f.write("2 %d %d %d %s\n" % (total_num_blocks + 3, j + 1, j + 1, value.n(digits=64)))

    # TODO: report error if problem infeasible

    def _run_sdp_solver(self, show_output=False, solver="csdp", use_initial_point=False):

        self.state("run_sdp_solver", "yes")

        previous_directory = os.getcwd()
        os.chdir(unicode(SAGE_TMP))

        if solver == "csdp":
            cmd = "%s %s sdp.out" % (cdsp_cmd, self._sdp_input_filename)

            if use_initial_point and self.state("write_sdp_initial_point_file") == "yes":
                cmd += " %s" % self._sdp_initial_point_filename

        elif solver == "dsdp":
            cmd = "%s %s -gaptol 1e-18 -print 1 -save sdp.out" % (dsdp_cmd, self._sdp_input_filename)

        elif solver == "sdpa":
            cmd = "%s -ds %s -o sdpa.out" % (sdpa_cmd, self._sdp_input_filename)

        elif solver == "sdpa_dd":
            cmd = "%s -ds %s -o sdpa.out" % (sdpa_dd_cmd, self._sdp_input_filename)

        elif solver == "sdpa_qd":
            cmd = "%s -ds %s -o sdpa.out" % (sdpa_qd_cmd, self._sdp_input_filename)

        else:
            raise ValueError("unknown solver.")

        sys.stdout.write("Running SDP solver...\n")

        # For maximization problems, the objective value returned by the SDP solver
        # must be negated.
        obj_value_factor = 1.0 if self._minimize else -1.0

        p = pexpect.spawn(cmd, timeout=60 * 60 * 24 * 7)
        obj_val = None
        self._sdp_solver_output = ""
        while True:
            if p.eof():
                break
            try:
                p.expect("\r\n")
                line = p.before.strip() + "\n"
                self._sdp_solver_output += line

                if show_output:
                    sys.stdout.write(line)

                if "Primal objective value:" in line:  # CSDP
                    obj_val = self._approximate_field(line.split()[-1]) * obj_value_factor
                elif "objValPrimal" in line:  # SDPA
                    obj_val = self._approximate_field(line.split()[-1]) * obj_value_factor
                elif "DSDP Solution" in line:  # DSDP: seems to print absolute value
                    obj_val = self._approximate_field(line.split()[-1])

            except pexpect.EOF:
                break

        p.close()
        self._sdp_solver_returncode = p.exitstatus

        sys.stdout.write("Returncode is %d. Objective value is %s.\n" % (
            self._sdp_solver_returncode, obj_val))

        # TODO: if program is infeasible, a returncode of 1 is given,
        # and output contains "infeasible"

        if "sdpa" in solver:

            with open("sdpa.out", "r") as inf:
                with open("sdp.out", "w") as f:

                    found, diagonal = False, False
                    t, row, col = 0, 1, 1

                    for line in inf:
                        if line[:6] == "yMat =":
                            break
                    else:
                        raise ValueError

                    for line in inf:

                        if line == "}":
                            break
                        elif line[:3] == "{ {":
                            t += 1
                            row = 1
                            diagonal = False
                        elif line[:2] == "{+" or line[:2] == "{-":
                            t += 1
                            row = 1
                            diagonal = True

                        line = line.replace("{", "")
                        line = line.replace("}", "")
                        col = 1
                        for a in line.split(","):
                            try:
                                v = a.strip()
                                vf = float(v)  # only done to see if we get ValueError
                                if diagonal:
                                    f.write("2 %d %d %d %s\n" % (t, row, col, v))
                                    row += 1
                                elif row <= col:
                                    f.write("2 %d %d %d %s\n" % (t, row, col, v))
                                col += 1
                            except ValueError:
                                pass

                        if col > 1:  # at least one number found...
                            row += 1

        self._sdp_output_filename = os.path.join(unicode(SAGE_TMP), "sdp.out")
        os.chdir(previous_directory)

    # TODO: read in dual solution

    def _read_sdp_output_file(self):

        self.state("read_solution", "yes")

        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_densities = len(self._densities)

        if self.state("set_block_matrix_structure") != "yes":
            self._set_block_matrix_structure()
        num_blocks = len(self._block_matrix_structure)

        with open(self._sdp_output_filename, "r") as f:

            self._sdp_Q_matrices = [matrix(self._approximate_field,
                len(self._flags[ti]), len(self._flags[ti])) for ti in range(num_types)]

            self._sdp_density_coeffs = [self._approximate_field(0) for i in range(num_densities)]

            for line in f:
                numbers = line.split()
                if numbers[0] != "2":
                    continue
                bi = int(numbers[1]) - 2
                if bi == num_blocks + 1:
                    j = int(numbers[2]) - 1
                    di = self._active_densities[j]
                    self._sdp_density_coeffs[di] = self._approximate_field(numbers[4])
                    continue
                if bi < 0 or bi >= num_blocks:
                    continue
                j = int(numbers[2]) - 1
                k = int(numbers[3]) - 1
                ti, size, offset = self._block_matrix_structure[bi]
                j += offset
                k += offset

                self._sdp_Q_matrices[ti][j, k] = self._approximate_field(numbers[4])
                self._sdp_Q_matrices[ti][k, j] = self._sdp_Q_matrices[ti][j, k]

        for ti in range(num_types):
            self._sdp_Q_matrices[ti].set_immutable()

    def check_solution(self, tolerance=1e-5, show_sorted=False, show_all=False):
        r"""
        Checks the approximate floating point bound given by the SDP solver, and determines which
        graphs appear to be sharp. The apparently sharp graphs are compared against the graphs
        that the construction forces to be sharp.

        INPUT:

         - ``tolerance`` - Number (default: 0.00001) If a graph has a coefficient whose absolute
           difference with the bound is less than ``tolerance``, then it is considered to be sharp.

         - ``show_sorted`` - Boolean (default: False). Whether to sort the sharp graphs according
           to their coefficients. If False, then they will be displayed in order of the graph
           indices.

          - ``show_all`` - Boolean (default: False). If True, then the coefficients of all the
            graphs will be displayed instead of just the sharp graphs. In this case, the graphs
            that appear to be sharp are annotated with an "S" and the graphs that are forced to be
            sharp by the construction are annotated with a "C". (If these sets are not identical,
            then there is a problem.)
        """
        self.state("check_solution", "yes")

        num_types = len(self._types)
        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        sys.stdout.write("Checking numerical bound...\n")

        fbounds = [sum([self._densities[j][i] * self._sdp_density_coeffs[j] for j in range(num_densities)]) for i in range(num_graphs)]

        for ti in self._active_types:
            for row in self._product_densities_arrays[ti]:
                gi, j, k, numer, denom = row
                d = Integer(numer) / Integer(denom)
                value = self._sdp_Q_matrices[ti][j, k]
                if j != k:
                    d *= 2
                if not self._minimize:
                    fbounds[gi] += d * value
                else:
                    fbounds[gi] -= d * value

        if not self._minimize:
            bound = max(fbounds)
        else:
            bound = min(fbounds)

        self._sdp_bounds = fbounds

        if self.state("set_construction") == "yes":

            if abs(bound - self._approximate_field(self._target_bound)) < tolerance:
                sys.stdout.write("Bound of %s appears to have been met.\n" % self._target_bound)
            else:
                sys.stdout.write("Warning: bound of %s appears to have not been met.\n" % self._target_bound)
                return
            sharp_graphs = self._sharp_graphs
        else:
            sharp_graphs = []  # set dummy sharp_graphs

        apparently_sharp_graphs = [gi for gi in range(num_graphs) if abs(fbounds[gi] - bound) < tolerance]

        if show_sorted or show_all:

            if not self._minimize:
                sorted_indices = sorted(range(num_graphs), key=lambda i: -fbounds[i])
            else:
                sorted_indices = sorted(range(num_graphs), key=lambda i: fbounds[i])

            for gi in sorted_indices:
                if gi in apparently_sharp_graphs:
                    sys.stdout.write("S")
                elif not show_all:
                    break
                else:
                    sys.stdout.write(" ")
                if gi in self._sharp_graphs:
                    sys.stdout.write("C")
                else:
                    sys.stdout.write(" ")
                sys.stdout.write(" %s : graph %d (%s) " % (fbounds[gi], gi, self._graphs[gi]))
                sys.stdout.write("\n")

        if self.state("set_construction") != "yes":
            return

        if not (show_sorted or show_all):
            sys.stdout.write("The following %d graphs appear to be sharp:\n" % len(apparently_sharp_graphs))
            for gi in apparently_sharp_graphs:
                sys.stdout.write("%.12f : graph %d (%s)\n" % (fbounds[gi], gi, self._graphs[gi]))

        extra_sharp_graphs = [gi for gi in apparently_sharp_graphs if not gi in self._sharp_graphs]
        missing_sharp_graphs = [gi for gi in self._sharp_graphs if not gi in apparently_sharp_graphs]

        if len(extra_sharp_graphs) > 0:
            sys.stdout.write("Warning: additional sharp graphs: %s\n" % (extra_sharp_graphs,))

        for gi in missing_sharp_graphs:
            sys.stdout.write("Warning: graph %d (%s) does not appear to be sharp.\n" % (gi, self._graphs[gi]))

    def import_solution(self, directory, complement=False):
        r"""
        Imports a solution found by Flagmatic 1.0 or 1.5.

        INPUT:

         - ``directory`` - the Flagmatic output directory, which must contain a flags.py file.

         - ``complement`` - Boolean (default: False). If True, then the solution will be assumed
           to be for the complementary problem. For example, if we are trying to minimize the
           density of k-cliques, the complementary problem is to minimize the density of
           independent sets of size k.
        """
        self.state("write_sdp_input_file", "yes")
        self.state("run_sdp_solver", "yes")

        num_graphs = len(self._graphs)
        num_types = len(self._types)
        num_densities = len(self._densities)

        if num_densities > 1:
            raise NotImplementedError

        sys.path.insert(0, directory)
        dont_write_bytecode = sys.dont_write_bytecode
        sys.dont_write_bytecode = True

        try:
            import flags
        except ImportError:
            sys.stdout.write("Cannot find flags.py in directory provided.\n")
            return

        # TODO: check admissible graphs, target bound, (sharps?).

        if flags.n != self._n:
            raise ValueError

        if self._flag_cls().oriented:
            if not ("oriented %d-graph" % self._flag_cls().r) in flags.description:
                raise ValueError
        else:
            if not ("%d-graph" % self._flag_cls().r) in flags.description:
                raise ValueError

        if flags.num_H != num_graphs:
            raise ValueError

        if flags.num_types != num_types:
            raise ValueError

        type_translations = []
        flag_translations = []

        for ti in range(num_types):
            tg = self._flag_cls(flags.types[ti])
            if complement:
                tg = tg.complement(True)
            for tj in range(num_types):
                if tg.is_labelled_isomorphic(self._types[tj]):
                    type_translations.append(tj)
                    break
            else:
                raise ValueError

            num_flags = len(self._flags[tj])
            ftr = []

            for fi in range(num_flags):
                fg = self._flag_cls(flags.flags[ti][fi])
                fg.t = tg.n
                if complement:
                    fg = fg.complement(True)
                for fj in range(num_flags):
                    if fg.is_labelled_isomorphic(self._flags[tj][fj]):
                        ftr.append(fj)
                        break
                else:
                    raise ValueError("solution has a flag that is not present.")

            flag_translations.append(ftr)

        self._sdp_Q_matrices = [matrix(self._approximate_field, len(self._flags[ti]),
                                len(self._flags[ti])) for ti in range(num_types)]

        try:
            f = open(directory + "/" + flags.out_filename, "r")
        except IOError:
            try:
                f = gzip.open(directory + "/" + flags.out_filename + ".gz", "rb")
            except IOError:
                print "Could not open %s or %s.gz" % (flags.out_filename, flags.out_filename)
                return

        for line in f:
            numbers = line.split()
            if numbers[0] != "2":
                continue
            ti = int(numbers[1]) - 2
            if ti >= 0 and ti < num_types:
                tj = type_translations[ti]
                if tj in self._active_types:
                    j = flag_translations[ti][int(numbers[2]) - 1]
                    k = flag_translations[ti][int(numbers[3]) - 1]
                    self._sdp_Q_matrices[tj][j, k] = numbers[4]
                    self._sdp_Q_matrices[tj][k, j] = self._sdp_Q_matrices[tj][j, k]

        f.close()

        self._sdp_density_coeffs = [1.0]

        for ti in range(num_types):
            self._sdp_Q_matrices[ti].set_immutable()

        sys.path.remove(directory)
        sys.dont_write_bytecode = dont_write_bytecode

    def make_exact(self, denominator=1024, meet_target_bound=True,
                   protect=None, use_densities=True, use_blocks=True, rank=None, show_changes=False,
                   check_exact_bound=True, diagonalize=True):
        r"""
        Makes an exact bound for the problem using the approximate floating point bound
        found by the SDP solver.

        INPUT:

         - ``denominator`` - Integer (default: 1024). The denominator to use when
            rounding. Higher numbers will cause the solution to be perturbed less, but can
            cause the resulting certificates to be larger.

         - ``meet_target_bound`` - Boolean (default: True). Determines whether the
           solution should be coerced to meet the target bound. If this is False, then a
           simpler method of rounding will be used. If there is no target bound, then this
           option will be disregarded, and the simpler method of rounding will be used.

         - ``protect`` - Array of integers or None (default: None). If an array of
            integers is given, then the entries of the matrices for the types with these
            indices will not be adjusted. (Using this option is not recommended.)

          - ``use_densities`` - Boolean (default: True). Whether to adjust the density
            coefficients to meet the target bound. This option will be disregarded if
            there is only one density.

          - ``use_blocks`` - Boolean (default: True). When computing the new basis for
             the solution's Q matrices, determines whether to give them a block structure
             with two blocks, using the invariant anti-invariant idea of Razborov.

          - ``rank`` - Integer or None (default: None). When computing the DR matrix,
             stop after ``rank`` columns have been found. This can save time in the
             case that the rank of the DR matrix is known (for example, from a previous
             run).

          - ``show_changes`` - Boolean (default: False). When meeting the target bound,
             display the changes being made to the matrix entries and the density
             coefficients.

          - ``check_exact_bound`` - Boolean (default: True). Whether to check the bound
             afterwards.

          - ``diagonalize`` - Boolean (default: True). Whether to diagonalize the Q
             matrices afterwards. If ``meet_target_bound`` is False, the Q matrices are
             always diagonalized.
        """

        if meet_target_bound and self.state("set_construction") != "yes":
            meet_target_bound = False
            sys.stdout.write("No target bound to meet.\n")

        if meet_target_bound:
            self.change_solution_bases(use_blocks=use_blocks)
            num_sharps = len(self._sharp_graphs)
        else:
            self._sdp_Qdash_matrices = self._sdp_Q_matrices
            # if non-tight, we don't have a separate diagonalization step
            self._exact_diagonal_matrices = []
            self._exact_r_matrices = []
            num_sharps = 0

        self.state("make_exact", "yes")

        num_types = len(self._types)
        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        if protect is None:
            protect = []

        q_sizes = [self._sdp_Qdash_matrices[ti].nrows() for ti in range(num_types)]

        def rationalize(f):
            return Integer(round(f * denominator)) / denominator

        sys.stdout.write("Rounding matrices")

        self._exact_Qdash_matrices = []

        for ti in range(num_types):

            if meet_target_bound:

                M = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
                for j in range(q_sizes[ti]):
                    for k in range(j, q_sizes[ti]):
                        value = rationalize(self._sdp_Qdash_matrices[ti][j, k])
                        if value != 0:
                            M[j, k] = value
                            M[k, j] = value

            else:

                try:
                    LF = numpy.linalg.cholesky(self._sdp_Qdash_matrices[ti])
                    # TODO: Consider using this:
                    # LF = self._sdp_Qdash_matrices[ti].cholesky_decomposition()
                except numpy.linalg.linalg.LinAlgError:
                # except ValueError:
                    sys.stdout.write("Could not compute Cholesky decomposition for type %d.\n" % ti)
                    return
                L = matrix(QQ, q_sizes[ti], q_sizes[ti], sparse=True)
                for j in range(q_sizes[ti]):
                    for k in range(j + 1):  # only lower triangle
                        L[j, k] = rationalize(LF[j, k])
                L.set_immutable()
                M = L * L.T
                if not meet_target_bound:
                    D = identity_matrix(QQ, q_sizes[ti], sparse=True)
                    D.set_immutable()
                    self._exact_diagonal_matrices.append(D)
                    self._exact_r_matrices.append(L)

            row_div = self._sdp_Qdash_matrices[ti].subdivisions()[0]
            M.subdivide(row_div, row_div)
            self._exact_Qdash_matrices.append(matrix(self._field, M))
            sys.stdout.write(".")
            sys.stdout.flush()

        sys.stdout.write("\n")

        self._exact_density_coeffs = [rationalize(self._sdp_density_coeffs[di]) for di in range(num_densities)]

        # TODO: Make density coeffs in each block sum to 1. Ever necessary when >1 density?
        for db in self._density_coeff_blocks:
            if len(db) == 1:
                self._exact_density_coeffs[db[0]] = Integer(1)

        # Force all inactive densities to be 0 (they should be anyway).
        for j in range(num_densities):
            if not j in self._active_densities:
                self._exact_density_coeffs[j] = Integer(0)

        if meet_target_bound:

            self.state("meet_target_bound", "yes")

            triples = [(ti, j, k) for ti in self._active_types for j in range(q_sizes[ti])
                       for k in range(j, q_sizes[ti])]

            num_triples = len(triples)
            triples.sort()
            triple_to_index = dict((triples[i], i) for i in range(num_triples))

            R = matrix(self._field, num_sharps, num_triples, sparse=True)

            sys.stdout.write("Constructing R matrix")

            # TODO: only use triples that correspond to middle blocks.

            for ti in self._active_types:

                Ds = [matrix(QQ, len(self._flags[ti]), len(self._flags[ti]))
                      for si in range(num_sharps)]

                for row in self._product_densities_arrays[ti]:
                    gi = row[0]
                    if not gi in self._sharp_graphs:
                        continue
                    si = self._sharp_graphs.index(gi)
                    j = row[1]
                    k = row[2]
                    value = Integer(row[3]) / Integer(row[4])
                    Ds[si][j, k] = value
                    Ds[si][k, j] = value

                if self.state("transform_solution") == "yes":
                    B = self._inverse_flag_bases[ti]
                    for si in range(num_sharps):
                        Ds[si] = B.T * Ds[si] * B

                for si in range(num_sharps):
                    for j in range(q_sizes[ti]):
                        for k in range(j, q_sizes[ti]):
                            trip = (ti, j, k)
                            value = Ds[si][j, k]
                            if j != k:
                                value *= 2
                            if self._minimize:
                                value *= -1
                            R[si, triple_to_index[trip]] = value

                sys.stdout.write(".")
                sys.stdout.flush()
            sys.stdout.write("\n")

            density_cols_to_use = []
            DR = matrix(self._field, num_sharps, 0)  # sparsity harms performance too much here
            EDR = DR.T

            sys.stdout.write("Constructing DR matrix")

            # Only if there is more than one density
            if num_densities > 1 and use_densities:

                for j in self._active_densities:

                    if not rank is None and DR.ncols() == rank:
                        break

                    new_col = matrix(QQ, [[self._densities[j][gi]] for gi in self._sharp_graphs])
                    if new_col.is_zero():
                        continue
                    EDR = EDR.stack(new_col.T)
                    EDR.echelonize()
                    if EDR[-1, :].is_zero():
                        EDR = EDR[:-1, :]
                        sys.stdout.write("~")
                        sys.stdout.flush()
                        continue

                    DR = DR.augment(new_col)
                    density_cols_to_use.append(j)
                    sys.stdout.write(".")
                    sys.stdout.flush()

                sys.stdout.write("\n")
                sys.stdout.write("DR matrix (density part) has rank %d.\n" % DR.ncols())

            col_norms = {}
            for i in range(num_triples):
                n = sum(x**2 for x in R.column(i))
                if n != 0:
                    col_norms[i] = n

            # Use columns with greatest non-zero norm - change minus to plus to
            # use smallest columns (not sure which is best, or maybe middle?)

            cols_in_order = sorted(col_norms.keys(), key = lambda i : -col_norms[i])
            cols_to_use = []

            for i in cols_in_order:

                if not rank is None and DR.ncols() == rank:
                    break

                ti, j, k = triples[i]
                if ti in protect:  # don't use protected types
                    continue
                new_col = R[:, i: i + 1]
                if new_col.is_zero():
                    continue
                EDR = EDR.stack(new_col.T)
                EDR.echelonize()
                if EDR[-1, :].is_zero():
                    EDR = EDR[:-1, :]
                    sys.stdout.write("~")
                    sys.stdout.flush()
                    continue

                DR = DR.augment(new_col)
                cols_to_use.append(i)
                sys.stdout.write(".")
                sys.stdout.flush()

            sys.stdout.write("\n")
            sys.stdout.write("DR matrix has rank %d.\n" % DR.ncols())

            T = matrix(self._field, num_sharps, 1)

            for si in range(num_sharps):

                gi = self._sharp_graphs[si]
                T[si, 0] = self._target_bound

                for j in range(num_densities):
                    if not j in density_cols_to_use:
                        T[si, 0] -= self._exact_density_coeffs[j] * self._densities[j][gi]

                for i in range(num_triples):
                    if not i in cols_to_use:
                        ti, j, k = triples[i]
                        T[si, 0] -= self._exact_Qdash_matrices[ti][j, k] * R[si, i]

            FDR = matrix(self._field, DR)
            try:
                X = FDR.solve_right(T)
            except ValueError:
                if rank is None:
                    raise ValueError("could not meet bound.")
                else:
                    raise ValueError("could not meet bound (try increasing the value of ``rank``).")

            RX = matrix(self._approximate_field, X.nrows(), 1)

            for i in range(len(density_cols_to_use)):
                di = density_cols_to_use[i]
                RX[i, 0] = self._exact_density_coeffs[di]
                self._exact_density_coeffs[di] = X[i, 0]

            for i in range(len(density_cols_to_use), X.nrows()):
                ti, j, k = triples[cols_to_use[i - len(density_cols_to_use)]]
                RX[i, 0] = self._sdp_Qdash_matrices[ti][j, k]
                self._exact_Qdash_matrices[ti][j, k] = X[i, 0]
                self._exact_Qdash_matrices[ti][k, j] = X[i, 0]

            if show_changes:
                for i in range(X.nrows()):
                    sys.stdout.write("%.11s -> %.11s " % (RX[i,0], RDF(X[i,0])))
                    if i < len(density_cols_to_use):
                        sys.stdout.write("(density %d)\n" % density_cols_to_use[i])
                    else:
                        sys.stdout.write("(matrix %d, entry [%d, %d])\n" % triples[cols_to_use[i - len(density_cols_to_use)]])

        for ti in range(num_types):
            self._exact_Qdash_matrices[ti].set_immutable()

        if check_exact_bound:
            self.check_exact_bound()

    def check_exact_bound(self, diagonalize=True):
        r"""
        Usually called by ``make_exact``. If the solution was transformed, then computes
        the Q matrices from the Q' matrices. If the solution was adjusted to meet the
        target bound, the Q matrices are checked (numerically) for negative eigenvalues.
        In all cases the bound is checked.

        If ``diagonalize`` is set to True, then ``diagonalize`` will be called at the
        end.
        """
        num_types = len(self._types)
        num_graphs = len(self._graphs)
        num_densities = len(self._densities)

        if num_densities > 1:
            negative_densities = [j for j in range(num_densities) if self._exact_density_coeffs[j] < 0]
            if len(negative_densities) > 0:
                sys.stdout.write("Warning! Densities %s have negative coefficients, so the bound is not valid.\n" % negative_densities)
                return
            sys.stdout.write("All density coefficients are non-negative.\n")

        # If we didn't try to meet the target bound, then the method of rounding is_exact
        # guaranteed to produce positive-semidefinite matrices.
        if self.state("meet_target_bound") == "yes":
            negative_types = []
            very_small_types = []
            for ti in self._active_types:
                if self._exact_Qdash_matrices[ti].nrows() > 0:
                    eigvals = sorted(numpy.linalg.eigvalsh(self._exact_Qdash_matrices[ti]))
                    if eigvals[0] < 0.0:
                        negative_types.append(ti)
                    elif eigvals[0] < 1e-6:
                        very_small_types.append(ti)

            if len(negative_types) > 0:
                sys.stdout.write("Warning! Types %s have negative eigenvalues, so the bound is not valid.\n" % negative_types)
                return
            sys.stdout.write("All eigenvalues appear to be positive.\n")
            if len(very_small_types) > 0:
                sys.stdout.write("Types %s have very small eigenvalues (but this is probably OK).\n" % very_small_types)

        self.state("check_exact", "yes")

        self._exact_Q_matrices = []

        if self.state("transform_solution") == "yes":
            for ti in range(num_types):
                B = self._inverse_flag_bases[ti]
                M = B * self._exact_Qdash_matrices[ti] * B.T
                self._exact_Q_matrices.append(M)
        else:
            self._exact_Q_matrices = self._exact_Qdash_matrices

        bounds = [sum([self._densities[j][i] * self._exact_density_coeffs[j]
                  for j in range(num_densities)]) for i in range(num_graphs)]

        for ti in self._active_types:
            for row in self._product_densities_arrays[ti]:
                gi, j, k, numer, denom = row
                d = Integer(numer) / Integer(denom)
                value = self._exact_Q_matrices[ti][j, k]
                if j != k:
                    value *= 2
                if not self._minimize:
                    bounds[gi] += d * value
                else:
                    bounds[gi] -= d * value

        if self._field == QQ:
            if not self._minimize:
                bound = max(bounds)
            else:
                bound = min(bounds)
        else:
            # Sorting doesn't currently work for number fields with embeddings, so use float approximation.
            # TODO: Check if Sage 5.0 fixes this.
            if not self._minimize:
                bound = max(bounds, key=lambda x: float(x))
            else:
                bound = min(bounds, key=lambda x: float(x))

        self._bounds = bounds
        self._bound = bound

        if self.state("meet_target_bound") != "yes":
            sys.stdout.write("Bound is %s (%s).\n" % (bound, bound.n(digits=10)))
            return

        if self._field == QQ:
            if not self._minimize:
                violators = [gi for gi in range(num_graphs) if bounds[gi] > self._target_bound]
            else:
                violators = [gi for gi in range(num_graphs) if bounds[gi] < self._target_bound]
        else:
            if not self._minimize:
                violators = [gi for gi in range(num_graphs) if float(bounds[gi]) > float(self._target_bound)]
            else:
                violators = [gi for gi in range(num_graphs) if float(bounds[gi]) < float(self._target_bound)]

        sys.stdout.write("Bound of %s attained by:\n" % self._target_bound)
        unexpectedly_sharp = []
        for gi in range(num_graphs):
            if bounds[gi] == self._target_bound:
                sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))
                if not gi in self._sharp_graphs:
                    unexpectedly_sharp.append(gi)

        if len(unexpectedly_sharp) > 0:
            sys.stdout.write("Warning: the following graphs unexpectedly attain the bound: %s\n"
                             % ", ".join([str(gi) for gi in unexpectedly_sharp]))

        if len(violators) > 0:
            sys.stdout.write("Bound violated by:")
            for gi in violators:
                sys.stdout.write("%s : graph %d (%s)\n" % (bounds[gi], gi, self._graphs[gi]))

        if diagonalize:
            self.diagonalize()

    def diagonalize(self):
        r"""
        For each matrix Q, produces a matrix R and a diagonal matrix M such that
        Q = R * M * R.T, where R.T denotes the transpose of R. Usually called from
        ``make_exact``. Note that if the solution has not been adjusted to meet a target
        bound, a simpler method of rounding is performed, and diagonalization is done
        at the same time.
        """

        self.state("diagonalize", "yes")

        self._exact_diagonal_matrices = []
        self._exact_r_matrices = []

        sys.stdout.write("Diagonalizing")

        for ti in range(len(self._types)):
            R, M = LDLdecomposition(self._exact_Qdash_matrices[ti])
            self._exact_diagonal_matrices.append(M)
            if self.state("transform_solution") == "yes":
                R = self._inverse_flag_bases[ti] * R
            self._exact_r_matrices.append(R)
            sys.stdout.write(".")
            sys.stdout.flush()
        sys.stdout.write("\n")

        # Q can now be computed as Q = R * M * R.T

        sys.stdout.write("Verifying")

        for ti in range(len(self._types)):
            Q = self._exact_r_matrices[ti] * self._exact_diagonal_matrices[ti] * self._exact_r_matrices[ti].T
            if Q != self._exact_Q_matrices[ti]:
                raise ValueError  # TODO: choose appropriate error
            sys.stdout.write(".")
            sys.stdout.flush()
        sys.stdout.write("\n")

    def describe(self):
        r"""
        Returns a human-readable description of the problem. This is used by
        ``make_certificate``.
        """

        description = self._flag_cls.description() + "; "
        description += "minimize " if self._minimize else "maximize "
        description += self._describe_density()

        forbidden = []
        for g in self._forbidden_graphs:
            forbidden.append(str(g))
        for g in self._forbidden_induced_graphs:
            forbidden.append("induced %s" % g)
        for pair in self._forbidden_edge_numbers:
            forbidden.append("induced %d.%d" % pair)
        if len(forbidden) > 0:
            description += "; forbid " + ", ".join(forbidden)

        return description

    def _describe_density(self):

        if len(self._density_graphs) == 0:
            return ""
        elif len(self._density_coeff_blocks) == 1 and len(self._density_coeff_blocks[0]) == 1:
            density = self._density_graphs[0]
            if len(density) == 1 and density[0][1] == 1:
                return "%s density" % density[0][0]
            else:
                return "density expression %s" % self._flag_cls.format_combination(density)
        else:
            return "combination of quantum graphs"

    def _augment_certificate(self, data):
        pass

    def write_certificate(self, filename):
        r"""
        Writes a certificate in a (hopefully) clear format, in JSON, to the file specified
        by ``filename``. For more information about the contents of the certificates, see
        the User's Guide.
        """

        self.state("write_certificate", "yes")

        def upper_triangular_matrix_to_list(M):
            return [list(M.row(i))[i:] for i in range(M.nrows())]

        def matrix_to_list(M):
            return [list(M.row(i)) for i in range(M.nrows())]

        if self.state("meet_target_bound") != "yes" or self.state("diagonalize") == "yes":
            qdash_matrices = self._exact_diagonal_matrices
            r_matrices = self._exact_r_matrices
        else:
            qdash_matrices = self._exact_Qdash_matrices
            r_matrices = self._inverse_flag_bases

        data = {
            "description": self.describe(),
            "bound": self._bound,
            "order_of_admissible_graphs": self._n,
            "number_of_admissible_graphs": len(self._graphs),
            "admissible_graphs": self._graphs,
            "number_of_types": len(self._types),
            "types": self._types,
            "numbers_of_flags": [len(L) for L in self._flags],
            "flags": self._flags,
            "qdash_matrices": [upper_triangular_matrix_to_list(M) for M in qdash_matrices],
            "r_matrices": [matrix_to_list(M) for M in r_matrices],
        }

        if len(self._density_graphs) == 1:
            data["admissible_graph_densities"] = self._densities[0]

        # Allow subclasses to add more things
        self._augment_certificate(data)

        def default_handler(O):
            # Only output an int if it is less than 2^53 (to be valid JSON).
            if O in ZZ and O < 9007199254740992:
                return int(Integer(O))
            return repr(O)

        try:
            with open(filename, "w") as f:
                json.dump(data, f, indent=4, default=default_handler)
            sys.stdout.write("Written certificate.\n")

        except IOError:
            sys.stdout.write("Cannot open file for writing.\n")


def ThreeGraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type 3-graph problem. For help
    with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(ThreeGraphFlag, order, **kwargs)


def GraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type graph problem. For help
    with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(GraphFlag, order, **kwargs)


def OrientedGraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type oriented graph problem. For
    help with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(OrientedGraphFlag, order, **kwargs)


def TwoMultigraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type 2-multigraph problem. For
    help with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(TwoMultigraphFlag, order, **kwargs)


def ThreeMultigraphProblem(order=None, **kwargs):
    r"""
    Returns a Problem object, that will represent a Turán-type 3-multigraph problem. For
    help with Problem objects, enter

    sage: help(Problem)
    """
    return Problem(ThreeMultigraphFlag, order, **kwargs)
