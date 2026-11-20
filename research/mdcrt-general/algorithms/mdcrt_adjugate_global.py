"""Solve general MD-CRT using one row HNF of the scaled adjugate MLC."""

from math import gcd, lcm
from operator import index
from time import perf_counter_ns

from mdcrt_adjugate_dual import Q, Z, adjugate_rows, column, eye, integer_matrix, rectangle, rows, vector
from mlc_mdcrt_exact.bidirectional import prepare


class AdjugateReduction:
    """Fixed scalar constraints shared by the two global HNF routes."""

    def __init__(self, matrices):
        matrices = [tuple(tuple(index(value) for value in row) for row in m) for m in matrices]
        if not matrices or not matrices[0]:
            raise ValueError('at least one nonempty matrix modulus is required')
        n = len(matrices[0])
        if any(len(m) != n or any(len(row) != n for row in m) for m in matrices):
            raise ValueError('matrix moduli must be square and have the same dimension')
        factor = gcd(*(value for m in matrices for row in m for value in row))
        if not factor:
            raise ValueError('matrix modulus must be nonsingular')
        self.n, self.count, self.factor = n, len(matrices), factor
        self.blocks = [adjugate_rows(Z([[value // factor for value in row] for row in m]))
                       for m in matrices]
        self.coefficients = [row for coefficients, _ in self.blocks for row in coefficients]
        self.moduli = [d for _, moduli in self.blocks for d in moduli]

    def rhs(self, residues):
        residues = [tuple(index(value) for value in r) for r in residues]
        if len(residues) != self.count or any(len(r) != self.n for r in residues):
            raise ValueError('one residue vector of the correct dimension is required per modulus')
        offset = tuple(value % self.factor for value in residues[0])
        rhs = []
        for residue, (coefficients, moduli) in zip(residues, self.blocks):
            reduced = []
            for value, base in zip(residue, offset):
                quotient, remainder = divmod(value - base, self.factor)
                if remainder:
                    return None
                reduced.append(quotient)
            for row, modulus in zip(coefficients, moduli):
                rhs.append(sum(a * s for a, s in zip(row, reduced)) % modulus)
        return offset, rhs

    def lift(self, offset, solution):
        return tuple(base + self.factor * value for base, value in zip(offset, solution))


class AuxiliaryHNF:
    """Prepare the fixed augmented integer equations of a reduced MD-CRT."""

    def __init__(self, reduction):
        self.reduction = reduction
        n, m = reduction.n, len(reduction.moduli)
        self.solver = None
        basis = eye(n)
        if m:
            self.solver = prepare({'kind': 'mlc', 'A': reduction.coefficients,
                                   'moduli': reduction.moduli}, 'hnf_tri')
            basis = Z(self.solver.B)
        self.B = tuple(map(tuple, rows(basis * reduction.factor)))
        self.work_shape = (m + n, m) if m else None

    def solve(self, residues):
        mapped = self.reduction.rhs(residues)
        if mapped is None:
            return None
        offset, rhs = mapped
        if self.solver is None:
            return offset
        solution = self.solver.solve(rhs)
        return None if solution is None else self.reduction.lift(offset, solution)


class CongruenceHNF:
    """Cache the fixed part of the row lattice and the first complete period basis."""

    def __init__(self, reduction):
        self.reduction = reduction
        self.common = lcm(*reduction.moduli)
        self.scales = [self.common // d for d in reduction.moduli]
        self.fixed = [[scale * a for a in row] for scale, row in zip(self.scales, reduction.coefficients)]
        self.identity = rows(eye(reduction.n + 1) * self.common)
        self.work_shape = (reduction.n + 1 + len(self.scales), reduction.n + 1) if self.scales else None
        self.inverse = None
        self.hnf_ns = 0
        self.B = None if self.scales else tuple(map(tuple, rows(eye(reduction.n) * reduction.factor)))

    def solve(self, residues):
        mapped = self.reduction.rhs(residues)
        if mapped is None:
            return None
        offset, rhs = mapped
        if not rhs:
            return offset
        n = self.reduction.n
        generators = self.identity + [row + [scale * b]
                                      for row, scale, b in zip(self.fixed, self.scales, rhs)]
        matrix = Z(generators)
        started = perf_counter_ns()
        hermite = matrix.hnf()
        self.hnf_ns += perf_counter_ns() - started
        if hermite[n, n] != self.common:
            return None
        basis = rectangle(hermite, 0, n, 0, n)
        if self.inverse is None:
            # For compatible RHS values this is the canonical HNF of the
            # projection onto the fixed first n coordinates.
            self.inverse = Q(basis).inv()
            self.B = tuple(map(tuple, rows(integer_matrix(
                self.inverse * (self.reduction.factor * self.common)))))
        rhs_column = column(hermite[i, n] for i in range(n))
        solution = vector(integer_matrix(self.inverse * Q(rhs_column)))
        return self.reduction.lift(offset, solution)


def solve_mdcrt_global(matrices, residues):
    """Return (particular solution, column basis) or None for inconsistent data."""
    solver = CongruenceHNF(AdjugateReduction(matrices))
    solution = solver.solve(residues)
    return None if solution is None else (solution, solver.B)
