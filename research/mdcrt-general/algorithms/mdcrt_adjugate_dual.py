"""Exact MD-CRT through adjugate congruences and a dual lattice."""

from math import gcd, lcm
from operator import index
from mlc_mdcrt_exact.bidirectional import Q, Z, column, eye, rows, vector


def integer_matrix(matrix):
    entries = matrix.entries()
    if any(value.denominator != 1 for value in entries):
        raise ArithmeticError('expected an exact integer matrix')
    return Z(matrix.nrows(), matrix.ncols(), [int(value) for value in entries])


def rectangle(matrix, first_row, row_count, first_column, column_count):
    return Z([[matrix[i, j] for j in range(first_column, first_column + column_count)]
              for i in range(first_row, first_row + row_count)])


def reduce_matrix(matrix, modulus):
    return Z(matrix.nrows(), matrix.ncols(), [value % modulus for value in matrix.entries()])


def adjugate_rows(matrix):
    """Return an equivalent list of scalar congruences for a single block."""
    determinant = int(matrix.det())
    if determinant == 0:
        raise ValueError('matrix modulus must be nonsingular')
    if abs(determinant) == 1:
        return [], []
    adjugate = integer_matrix(matrix.inv() * determinant)
    coefficients, moduli = [], []
    for row in rows(adjugate):
        factor = gcd(determinant, *row)
        modulus = abs(determinant) // factor
        if modulus == 1:
            continue
        reduced = [(value // factor) % modulus for value in row]
        if modulus == abs(determinant):
            # The resulting map onto Z/modulus has the same kernel index
            # as the original block, so this one row defines the block.
            return [reduced], [modulus]
        coefficients.append(reduced)
        moduli.append(modulus)
    return coefficients, moduli


class AdjugateDualMDCRT:
    """Prepare a system of full-rank square integer matrix congruences.

    B is a complete integer column basis of the homogeneous solutions.
    solve(residues) returns one integer solution or None if incompatible.
    """

    def __init__(self, matrices):
        matrices = tuple(tuple(tuple(index(value) for value in row) for row in m) for m in matrices)
        if not matrices or not matrices[0]:
            raise ValueError('at least one nonempty matrix modulus is required')
        n = len(matrices[0])
        if any(len(m) != n or any(len(row) != n for row in m) for m in matrices):
            raise ValueError('matrix moduli must be square and have the same dimension')
        factor = gcd(*(v for m in matrices for row in m for v in row))
        if factor == 0:
            raise ValueError('matrix modulus must be nonsingular')

        self.dimension, self.count, self.factor = n, len(matrices), factor
        self.blocks = [adjugate_rows(Z([[v // factor for v in row] for row in m])) for m in matrices]
        self.row_count = sum(len(moduli) for _, moduli in self.blocks)
        self.modulus = lcm(*(d for _, moduli in self.blocks for d in moduli))
        identity = eye(n)
        basis = identity
        previous = 1
        stages = []

        for coefficients, moduli in self.blocks:
            if not moduli:
                stages.append(None)
                continue
            current = lcm(previous, *moduli)
            scaled = [[(current // d) * v for v in row] for row, d in zip(coefficients, moduli)]
            generator = Z(rows(basis * (current // previous)) + scaled)
            hermite, transform = generator.hnf(transform=True)
            basis = rectangle(hermite, 0, n, 0, n)
            old = rectangle(transform, 0, n, 0, n)
            new = rectangle(transform, 0, n, n, len(moduli))
            contribution = new * Z(scaled) * (self.modulus // current)
            stages.append((reduce_matrix(old, self.modulus), reduce_matrix(contribution, self.modulus)))
            previous = current

        # Compile the fixed RHS lift by composing the block maps backwards.
        carry = identity
        lifts = []
        for stage in reversed(stages):
            if stage is None:
                lifts.append([[0] * n for _ in range(n)])
                continue
            old, contribution = stage
            lifts.append(rows(reduce_matrix(carry * contribution, self.modulus)))
            carry = reduce_matrix(carry * old, self.modulus)
        self.lift = Z([[value for block in reversed(lifts) for value in block[i]] for i in range(n)])

        self.dual_basis = integer_matrix(Q(basis).inv() * self.modulus)
        self.B = tuple(tuple(row) for row in rows(self.dual_basis * factor))
        self.check_rows = Z([row for coefficients, _ in self.blocks for row in coefficients]) if self.row_count else None
        self.check_moduli = tuple(d for _, moduli in self.blocks for d in moduli)

    def solve(self, residues):
        residues = tuple(tuple(index(value) for value in r) for r in residues)
        if len(residues) != self.count or any(len(r) != self.dimension for r in residues):
            raise ValueError('one residue vector of the correct dimension is required per modulus')
        g = self.factor
        if g == 1:
            offset = (0,) * self.dimension
            reduced = residues
        else:
            offset = tuple(value % g for value in residues[0])
            reduced = []
            for residue in residues:
                quotients = []
                for value, base in zip(residue, offset):
                    quotient, remainder = divmod(value - base, g)
                    if remainder:
                        return None
                    quotients.append(quotient)
                reduced.append(quotients)
        if not self.row_count:
            return offset

        lifted = self.lift * column(value for residue in reduced for value in residue)
        scaled = self.dual_basis * column(value % self.modulus for value in vector(lifted))
        y = []
        for value in vector(scaled):
            quotient, remainder = divmod(value, self.modulus)
            if remainder:
                return None
            y.append(quotient)

        # The lift solves a basis of the row lattice. Check all original scalar
        # congruences to reject inconsistent assignments of their RHS values.
        computed = vector(self.check_rows * column(y))
        cursor = 0
        for residue, (coefficients, moduli) in zip(reduced, self.blocks):
            for row, modulus in zip(coefficients, moduli):
                rhs = sum(a * v for a, v in zip(row, residue))
                if (computed[cursor] - rhs) % modulus:
                    return None
                cursor += 1
        return tuple(base + g * value for base, value in zip(offset, y))
