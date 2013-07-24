"""

Copyright (c) 2011, E. R. Vaughan. All rights reserved.

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

import fractions
import getopt
import itertools
import json
import re
import sys

HELP_TEXT = """
Usage: inspect_certificate.py CERTIFICATE OPTIONS
Possible options:
--help                       Display this message.
--admissible-graphs          Display the admissible graphs.
--flags                      Display the types and flags.
--r-matrices                 Display the R matrices.
--qdash-matrices             Display the Q' matrices.
--q-matrices                 Compute and display the Q matrices.
--pair-densities             Compute and display the flag pair densities.
--verify-bound               Verify the bound.
--sharp-graphs               Display the admissible graphs that are sharp.
--flag-algebra-coefficients  Display each admissible graph's flag algebra coefficient.
"""

try:
    from sage.all import *
    using_sage = True
except ImportError:
    using_sage = False

should_print_help = False

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], "", ["help", "admissible-graphs",
                                                      "flags", "r-matrices", "qdash-matrices",
                                                      "pair-densities", "q-matrices", "verify-bound",
                                                      "sharp-graphs", "flag-algebra-coefficients"])

except getopt.GetoptError:
    should_print_help = True
    opts, args = ((), ())

if len(args) != 1:
    should_print_help = True

action = ""

for o, a in opts:
    if o == "--help":
        should_print_help = True
    elif o == "--admissible-graphs":
        action = "print admissible graphs"
    elif o == "--flags":
        action = "print flags"
    elif o == "--r-matrices":
        action = "print r_matrices"
    elif o == "--qdash-matrices":
        action = "print qdash_matrices"
    elif o == "--pair-densities":
        action = "print pair densities"
    elif o == "--q-matrices":
        action = "print q_matrices"
    elif o == "--verify-bound":
        action = "verify bound"
    elif o == "--sharp-graphs":
        action = "print sharp graphs"
    elif o == "--flag-algebra-coefficients":
        action = "print flag algebra coefficients"

if should_print_help:
    print HELP_TEXT
    sys.exit(0)

certificate_filename = args[0]

try:
    if certificate_filename[-3:] == ".gz":
        import gzip
        certf = gzip.open(certificate_filename)
    elif certificate_filename[-4:] == ".bz2":
        import bz2
        certf = bz2.BZ2File(certificate_filename)
    else:
        certf = open(certificate_filename)
except IOError:
    sys.stdout.write("Could not open certificate.\n")
    sys.exit(1)

certificate = json.load(certf)

print 'Problem is "{}".'.format(certificate["description"])
print 'Claimed bound is {}.'.format(certificate["bound"])

minimize = "minimize" in certificate["description"]

if using_sage:
    x = polygen(QQ)
else:
    x = None

if "field" in certificate.keys():

    if certificate["field"] != "RationalField()":

        print 'Field used is "{}".'.format(certificate["field"])

        if not using_sage:
            print "This script must be run using Sage's python, as it does not use the rational field."
            sys.exit(1)

        try:
            K = sage_eval(certificate["field"], locals={'x': x})
        except AttributeError:
            K = RationalField()
        x = K.gen()

if action == "":
    sys.exit(0)

if "3-graph" in certificate["description"]:
    edge_size = 3
elif "2-graph" in certificate["description"]:
    edge_size = 2
else:
    print "Unsupported graph kind."
    sys.exit(1)

oriented = "oriented" in certificate["description"]


class Flag(object):

    edge_size = edge_size
    oriented = oriented

    __slots__ = ("n", "t", "edges")

    def __init__(self, s=None):

        if not s:
            self.n = 0
            self.t = 0
            self.edges = ()
            return

        m = re.match(r'(\d+):(\d*)(?:\((\d+)\)|)', s)
        if not m:
            raise ValueError
        n = int(m.group(1))
        t = int(m.group(3)) if m.group(3) else 0
        if not (0 <= n <= 9 and 0 <= t <= n):
            raise ValueError

        edges = []
        for i in range(0, len(m.group(2)), self.edge_size):
            edges.append(tuple(map(int, m.group(2)[i:i + self.edge_size])))

        self.n = n
        self.t = t
        self.edges = tuple(edges)

    def __repr__(self):
        s = "{}:{}".format(self.n, "".join(str(x) for e in self.edges for x in e))
        if self.t > 0:
            s += "({})".format(self.t)
        return s

    def __eq__(self, other):
        return self.n == other.n and self.t == other.t and self.edges == other.edges

    def __ne__(self, other):
        return not self.__eq__(other)

    def induced_subgraph(self, S):
        if not all(0 <= x <= self.n for x in S):
            raise ValueError
        good_edges = [e for e in self.edges if all(x in S for x in e)]
        p = [0 for _ in range(self.n + 1)]
        for i, x in enumerate(S):
            p[x] = i + 1
        h = Flag()
        h.n = len(S)
        h.t = 0
        if self.oriented:
            h.edges = tuple(sorted([tuple([p[x] for x in e]) for e in good_edges]))
        else:
            h.edges = tuple(sorted([tuple(sorted([p[x] for x in e])) for e in good_edges]))
        return h

    def minimal_isomorph(self):
        min_edges = self.edges
        for p in itertools.permutations(range(self.t + 1, self.n + 1), self.n - self.t):
            pt = tuple(range(1, self.t + 1)) + p
            if self.oriented:
                edges = tuple(sorted([tuple([pt[x - 1] for x in e]) for e in self.edges]))
            else:
                edges = tuple(sorted([tuple(sorted([pt[x - 1] for x in e])) for e in self.edges]))
            if edges < min_edges:
                min_edges = edges
        h = Flag()
        h.n = self.n
        h.t = self.t
        h.edges = min_edges
        return h

    def induced_flag(self, tv, ov):
        type_vertices = list(tv)
        other_vertices = list(ov)
        edges = []
        for e in self.edges:
            if any(not (x in tv or x in ov) for x in e):
                continue
            edge = []
            for x in e:
                if x in tv:
                    edge.append(1 + type_vertices.index(x))
                else:
                    edge.append(len(tv) + 1 + other_vertices.index(x))
            edges.append(tuple(edge))
        h = Flag()
        h.n = len(tv) + len(ov)
        h.t = len(tv)
        h.edges = tuple(edges)
        return h


admissible_graphs = [Flag(s) for s in certificate["admissible_graphs"]]
types = [Flag(s) for s in certificate["types"]]
flags = [[Flag(s) for s in f] for f in certificate["flags"]]

if action == "print admissible graphs":

    print "There are {} admissible graphs:".format(len(admissible_graphs))
    for i, g in enumerate(admissible_graphs):
        print "{}. {}".format(i + 1, g)
    sys.exit(0)

if action == "print flags":
    for ti, _type in enumerate(types):
        nf = len(flags[ti])
        print "Type {} ({}) has {} flags:".format(ti + 1, _type, nf)
        for i in range(nf):
            print "   {}. {}".format(i + 1, flags[ti][i])
    sys.exit(0)


def stringify(a):
    if a is None:
        return "Not used."
    if isinstance(a, list):
        return map(stringify, a)
    s = str(a)
    try:
        n = int(s)
        return n
    except ValueError:
        return s


if action == "print r_matrices":

    for ti, _type in enumerate(types):
        print "R matrix for type {} ({}):".format(ti + 1, _type)
        print "{}".format(stringify(certificate["r_matrices"][ti]))
    sys.exit(0)

if action == "print qdash_matrices":

    for ti, _type in enumerate(types):
        print "Q' matrix for type {} ({}):".format(ti + 1, _type)
        print "{}".format(stringify(certificate["qdash_matrices"][ti]))
    sys.exit(0)


print "Computing Q matrices..."

Qs = []

for ti, _type in enumerate(types):

    if certificate["qdash_matrices"][ti] is None:
        Qs.append(None)
        continue

    if using_sage:
        QD = [[sage_eval(str(s), locals={'x': x}) for s in row] for row in certificate["qdash_matrices"][ti]]
    else:
        QD = [[fractions.Fraction(s) for s in row] for row in certificate["qdash_matrices"][ti]]

    if certificate["r_matrices"][ti] is None:
        Qs.append(QD)
        continue

    if using_sage:
        R = [[sage_eval(str(s), locals={'x': x}) for s in row] for row in certificate["r_matrices"][ti]]
    else:
        R = [[fractions.Fraction(s) for s in row] for row in certificate["r_matrices"][ti]]

    nq = len(QD)
    nf = len(flags[ti])

    Q = [[0 for j in range(i, nf)] for i in range(nf)]
    for l in range(nq):
        for k in range(l, nq):
            qlk = QD[l][k - l]
            if qlk != 0:
                for i in range(nf):
                    for j in range(i, nf):
                        Q[i][j - i] += R[i][l] * R[j][k] * qlk
                        if k != l:
                            Q[i][j - i] += R[i][k] * R[j][l] * qlk
    Qs.append(Q)


if action == "print q_matrices":

    for ti, _type in enumerate(types):
        print "Q matrix for type {} ({}):".format(ti + 1, _type)
        print "{}".format(stringify(Qs[ti]))
    sys.exit(0)

print "Computing pair densities..."

pair_densities = {}
n = certificate["order_of_admissible_graphs"]

for ti, _type in enumerate(types):
    nf = len(flags[ti])
    s = _type.n
    if len(set(flag.n for flag in flags[ti])) != 1:
        raise ValueError("Flags for a given type must all have the same order.")
    m = flags[ti][0].n
    for g in admissible_graphs:
        pairs_found = [[0 for j in range(k, nf)] for k in range(nf)]
        total = 0
        for p in itertools.permutations(range(1, n + 1)):
            if p[s] > p[m]:
                continue
            is_good_permutation = True
            for c in range(s, n):
                if c in (s, m, 2 * m - s):
                    continue
                if p[c - 1] > p[c]:
                    is_good_permutation = False
                    break
            if not is_good_permutation:
                continue
            total += 1

            it = g.induced_subgraph(p[:s])
            if it != _type:
                continue

            f1 = g.induced_flag(p[:s], p[s: m]).minimal_isomorph()
            f2 = g.induced_flag(p[:s], p[m: 2 * m - s]).minimal_isomorph()

            if1 = flags[ti].index(f1)
            if2 = flags[ti].index(f2)

            if if1 <= if2:
                pairs_found[if1][if2 - if1] += 1
            else:
                pairs_found[if2][if1 - if2] += 1

        for k in range(nf):
            for j in range(k, nf):
                pf = pairs_found[k][j - k]
                if pf > 0:
                    if j == k:
                        pf *= 2
                    if using_sage:
                        pair_densities[(ti, admissible_graphs.index(g), k, j)] = Integer(pf) / (total * 2)
                    else:
                        pair_densities[(ti, admissible_graphs.index(g), k, j)] = fractions.Fraction(pf, total * 2)

if action == "print pair densities":

    for i, g in enumerate(admissible_graphs):
        print "Pair densities for admissible graph {} ({}):".format(i + 1, g)

        for ti, _type in enumerate(types):
            print "   Non-zero densities for type {} ({}):".format(ti + 1, _type)

            for key, d in pair_densities.items():
                if key[:2] == (ti, i):
                    print "      Flags {} and {} ({} and {}): {}".format(key[2] + 1, key[3] + 1,
                                                                         flags[ti][key[2]], flags[ti][key[3]], d)
    sys.exit(0)

print "Computing bound..."

if certificate["admissible_graph_densities"] and isinstance(certificate["admissible_graph_densities"][0], list):
    graph_densities = certificate["admissible_graph_densities"]
    density_coefficients = certificate["density_coefficients"]
else:
    graph_densities = [certificate["admissible_graph_densities"]]
    density_coefficients = [1]

bounds = [0 for _ in certificate["admissible_graphs"]]
for di, sdc in enumerate(density_coefficients):
    if using_sage:
        dc = sage_eval(str(sdc), locals={'x': x})
    else:
        dc = fractions.Fraction(sdc)
    for i, s in enumerate(graph_densities[di]):
        if using_sage:
            bounds[i] += dc * sage_eval(str(s), locals={'x': x})
        else:
            bounds[i] += dc * fractions.Fraction(s)

for key in pair_densities.keys():
    t, i, j, k = key
    if Qs[t] is None:  # check that type is used
        continue
    d = pair_densities[key]
    v = Qs[t][j][k - j]
    if minimize:
        v *= -1
    if j == k:
        bounds[i] += d * v
    else:
        bounds[i] += d * v * 2

bound = min(bounds) if minimize else max(bounds)

print "Bound is {}.".format(bound)

if action == "print sharp graphs":

    print "Sharp graphs are:"
    for i, g in enumerate(admissible_graphs):
        if bounds[i] == bound:
            print "{}. {}".format(i + 1, g)
    sys.exit(0)

if action == "print flag algebra coefficients":

    print "There are {} admissible graphs:".format(len(admissible_graphs))
    for i, g in enumerate(admissible_graphs):
        print "{}. ({}) : {}".format(i + 1, g, bounds[i])
    sys.exit(0)
