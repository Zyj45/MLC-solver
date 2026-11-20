"""Experimental exact adapters for the frozen bidirectional comparison.

All lattices are column lattices. No benchmark certificates enter prepare().
The production APIs and historical experiments are deliberately untouched.
"""
from dataclasses import dataclass
from math import gcd, lcm
from pathlib import Path
from time import perf_counter_ns

from flint import fmpz_mat as Z, fmpq_mat as Q, ctx
ctx.threads = 1

MLC_METHODS = ('hnf_tri', 'hnf_weight', 'row_md_tri', 'row_md_weight')
MD_METHODS = ('global_tri', 'global_weight', 'cascade_tri', 'cascade_weight',
              'adj_raw_tri', 'adj_raw_weight', 'adj_reduced_tri', 'adj_reduced_weight')


def bits(a):
    return max((abs(int(v)).bit_length() for v in a.entries()), default=0)


def rows(a):
    return [[int(a[i, j]) for j in range(a.ncols())] for i in range(a.nrows())]


def eye(n):
    return Z([[int(i == j) for j in range(n)] for i in range(n)])


def block_rows(a, start, count):
    return Z([[a[i, j] for j in range(a.ncols())] for i in range(start, start + count)])


def rational_join(a, b):
    return Q([[a[i, j] for j in range(a.ncols())] +
              [b[i, j] for j in range(b.ncols())] for i in range(a.nrows())])


def decomposition(c, projection, kernel=True, observe=False):
    """Return G, first projected rows of Z_part/Z_kernel: C Z=[G|0]."""
    h, t = c.transpose().hnf(transform=True)
    rank, width = c.nrows(), c.ncols()
    g = Z([[h[j, i] for j in range(rank)] for i in range(rank)])
    p = Z([[t[j, i] for j in range(rank)] for i in range(projection)])
    k = (Z([[t[j, i] for j in range(rank, width)] for i in range(projection)])
         if kernel else None)
    return g, p, k, bits(t) if observe else 0


def triangular(g, rhs):
    values = []
    for i, b in enumerate(rhs):
        q, rem = divmod(int(b) - sum(row * v for row, v in zip(g[i][:i], values)), g[i][i])
        if rem:
            return None
        values.append(q)
    return values


def column(values):
    return Z([[int(v)] for v in values])


def vector(a):
    return tuple(int(a[i, 0]) for i in range(a.nrows()))


@dataclass
class Compiled:
    numerator: object
    denominator: tuple

    @classmethod
    def make(cls, a):
        ds = tuple(lcm(*(int(a[i, j].denominator) for j in range(a.ncols())))
                   for i in range(a.nrows()))
        nums = Z([[int(a[i, j] * ds[i]) for j in range(a.ncols())]
                  for i in range(a.nrows())])
        return cls(nums, ds)

    def apply(self, rhs):
        v = vector(self.numerator * column(rhs))
        answer = []
        for a, d in zip(v, self.denominator):
            q, rem = divmod(a, d)
            if rem:
                return None
            answer.append(q)
        return tuple(answer)

    def max_bits(self):
        return max(bits(self.numerator), max((d.bit_length() for d in self.denominator), default=0))


@dataclass
class Prepared:
    kind: str
    method: str
    task: str
    original: dict
    B: tuple | None
    data: object
    checker: object
    metadata: dict

    def solve(self, rhs):
        if self.kind == 'mlc':
            flat = tuple(int(v) for v in rhs)
        else:
            flat = tuple(int(v) for r in rhs for v in r)
        mode, payload, mapping = self.data
        mapped = vector(mapping * column(flat)) if mapping is not None else flat
        if mode == 'hnf':
            g, part = payload
            y = triangular(g, mapped)
            answer = None if y is None else vector(part * column(y))
        elif mode == 'weight':
            answer = payload.apply(mapped)
            if answer is not None and not self.check(answer, flat):
                answer = None
        else:
            singles, stages = payload
            if singles is not None:
                ps = []
                for value, (g, beta) in zip(flat, singles):
                    q, rem = divmod(value, g)
                    if rem:
                        return None
                    ps.append(tuple(q * b for b in beta))
            else:
                n = len(self.original['matrices'][0])
                ps = [flat[i:i+n] for i in range(0, len(flat), n)]
            answer = ps[0]
            for p, (g, lift) in zip(ps[1:], stages):
                y = triangular(g, [a - b for a, b in zip(p, answer)])
                if y is None:
                    return None
                delta = vector(lift * column(y))
                answer = tuple(a + b for a, b in zip(answer, delta))
        return answer

    def check(self, x, flat):
        if self.kind == 'mlc':
            ax = vector(self.checker * column(x))
            return all((a - b) % m == 0 for a, b, m in zip(ax, flat, self.original['moduli']))
        n = len(x)
        return all(check.apply([a - b for a, b in zip(x, flat[i*n:(i+1)*n])]) is not None
                   for i, check in enumerate(self.checker))


def augmented(a, mods):
    r = len(mods)
    return Z([list(row) + [-mods[i] if i == j else 0 for j in range(r)]
              for i, row in enumerate(a)])


def md_global(matrices):
    n, k = len(matrices[0]), len(matrices)
    out = []
    for t, m in enumerate(matrices):
        for i in range(n):
            row = [int(i == j) for j in range(n)] + [0] * (n*k)
            row[n+t*n:n+(t+1)*n] = [-v for v in m[i]]
            out.append(row)
    return Z(out)


def prepare(problem, method, task='full', *, profile=False):
    start = perf_counter_ns()
    watch = bits if profile else lambda a: 0
    kind = problem['kind']
    if task not in ('particular', 'full'):
        raise ValueError('invalid task')
    if method not in (MLC_METHODS if kind == 'mlc' else MD_METHODS + ('formula16',)):
        raise ValueError('invalid method')
    weight = method.endswith('weight') or method == 'formula16'
    n = len(problem['A'][0]) if kind == 'mlc' else len(problem['matrices'][0])
    meta = {'method': method, 'task': task, 'native_threads': 1, 'work_shapes': [],
            'transform_checkpoint_bits': 0, 'intermediate_basis_checkpoint_bits': 0,
            'weight_bits': 0, 'weight_compile_ns': 0, 'adjugate_bits': 0,
            'rows_before': 0, 'rows_after': 0, 'internal_peak_observed': False,
            'diagnostic_profile': profile}
    if kind == 'mlc':
        a = [[int(v) for v in row] for row in problem['A']]
        mods = [int(v) for v in problem['moduli']]
        if not a or len(mods) != len(a) or any(len(row) != n for row in a) or any(m < 1 for m in mods):
            raise ValueError('invalid MLC')
        checker = Z(a)
        meta['input_bits'] = max(watch(checker), max(m.bit_length() for m in mods)) if profile else None
    else:
        matrices = [Z(m) for m in problem['matrices']]
        if not matrices or any(m.nrows() != n or m.ncols() != n or not m.det() for m in matrices):
            raise ValueError('invalid MD-CRT')
        meta['input_bits'] = max(watch(m) for m in matrices) if profile else None
        meta['determinant_bits'] = [abs(int(m.det())).bit_length() for m in matrices] if profile else None
        checker = None
    meta['validation_ns'] = perf_counter_ns() - start
    conv_start = perf_counter_ns()
    mapping = None
    singles = None
    md_matrices = None
    source_vectors = None
    if kind == 'mlc' and method.startswith('row_md'):
        singles, md_matrices, source_vectors = [], [], []
        for row, m in zip(a, mods):
            g, beta, bmat, tb = decomposition(Z([row + [-m]]), n, observe=profile)
            gi = int(g[0, 0])
            singles.append((gi, vector(beta)))
            source_vectors.append(Q(beta) / gi)
            md_matrices.append(bmat)
            meta['transform_checkpoint_bits'] = max(meta['transform_checkpoint_bits'], tb)
        meta['rows_before'] = meta['rows_after'] = len(a)
    elif kind == 'md' and method.startswith('adj_'):
        a, mods, maps = [], [], []
        for t, m in enumerate(matrices):
            det = int(m.det())
            adj = Q(m.inv()) * det
            if profile:meta['adjugate_bits'] = max(meta['adjugate_bits'], watch(Z(rows(adj))))
            for i in range(n):
                coeff = [int(adj[i, j]) for j in range(n)]
                modulus = abs(det)
                if method.startswith('adj_reduced'):
                    d = gcd(modulus, *coeff)
                    coeff, modulus = [v // d for v in coeff], modulus // d
                    if modulus == 1:
                        continue
                a.append(coeff)
                mods.append(modulus)
                maps.append([0]*(t*n) + coeff + [0]*((len(matrices)-t-1)*n))
        meta['rows_before'], meta['rows_after'] = len(matrices)*n, len(a)
        if not a:
            a, mods, maps = [[0]*n], [1], [[0]*(len(matrices)*n)]
        mapping = Z(maps)
    elif kind == 'md' and method.startswith('cascade'):
        md_matrices = matrices
    meta['conversion_ns'] = perf_counter_ns() - conv_start
    prep_start = perf_counter_ns()
    bfinal = None
    if method == 'formula16':
        ps = [int(m[0,0]) for m in matrices]
        if any(p<=1 or m != eye(n)*p for p,m in zip(ps,matrices)) or any(gcd(p,q)!=1 for i,p in enumerate(ps) for q in ps[i+1:]):
            raise ValueError('formula16 control requires scalar pairwise-coprime matrix moduli')
        total = 1
        for p in ps: total *= p
        ws = [total//p * pow(total//p,-1,p) for p in ps]
        payload = Compiled.make(Q([[ws[b] if i==j else 0 for b in range(len(ps)) for j in range(n)] for i in range(n)]))
        bfinal, mode = eye(n)*total, 'weight'
        meta['weight_bits'] = payload.max_bits() if profile else None
        meta['work_coefficient_bits'] = 0
    elif md_matrices is None:
        c = (md_global(problem['matrices']) if kind == 'md' and method.startswith('global')
             else augmented(a, mods))
        meta['work_shapes'].append([c.nrows(), c.ncols()])
        meta['work_coefficient_bits'] = watch(c)
        g, part, bfinal, tb = decomposition(c, n, task == 'full', observe=profile)
        meta['transform_checkpoint_bits'] = max(meta['transform_checkpoint_bits'], tb)
        if weight:
            t0 = perf_counter_ns()
            # Solve G^T W^T=P^T rather than materializing a full inverse.
            w = Q(g.transpose()).solve(Q(part.transpose())).transpose()
            payload = Compiled.make(w)
            meta['weight_compile_ns'] += perf_counter_ns() - t0
            meta['weight_bits'] = payload.max_bits() if profile else None
            mode = 'weight'
        else:
            payload, mode = (rows(g), part), 'hnf'
    else:
        bfinal, stages = md_matrices[0], []
        w = source_vectors[0] if source_vectors is not None else Q(eye(n))
        for i, m in enumerate(md_matrices[1:], 1):
            c = Z([br + [-v for v in mr] for br, mr in zip(rows(bfinal), rows(m))])
            meta['work_shapes'].append([n, 2*n])
            meta['work_coefficient_bits'] = max(meta.get('work_coefficient_bits', 0), watch(c))
            g, part, kernel, tb = decomposition(c, n, observe=profile)
            lift = bfinal * part
            stages.append((rows(g), lift))
            if weight:
                t0 = perf_counter_ns()
                k = Q(g.transpose()).solve(Q(lift.transpose())).transpose()
                new = k * source_vectors[i] if source_vectors is not None else k
                w = rational_join((Q(eye(n)) - k) * w, new)
                meta['weight_compile_ns'] += perf_counter_ns() - t0
            bfinal = bfinal * kernel
            meta['transform_checkpoint_bits'] = max(meta['transform_checkpoint_bits'], tb)
            meta['intermediate_basis_checkpoint_bits'] = max(meta['intermediate_basis_checkpoint_bits'], watch(bfinal))
        if weight:
            t0 = perf_counter_ns()
            payload = Compiled.make(w)
            meta['weight_compile_ns'] += perf_counter_ns() - t0
            meta['weight_bits'] = payload.max_bits() if profile else None
            mode = 'weight'
        else:
            payload, mode = (singles, stages), 'cascade'
    if kind == 'md' and weight:
        checker = tuple(Compiled.make(m.inv()) for m in matrices)
    meta['preparation_including_weight_ns'] = perf_counter_ns() - prep_start
    output_start = perf_counter_ns()
    basis = tuple(tuple(row) for row in rows(bfinal)) if task == 'full' else None
    meta['basis_bits'] = watch(bfinal) if bfinal is not None and profile else None
    meta['output_setup_ns'] = perf_counter_ns() - output_start
    meta['prepare_total_ns'] = perf_counter_ns() - start
    # The stage time includes weight compilation; totals must not add it twice.
    meta['preparation_excluding_weight_ns'] = meta['preparation_including_weight_ns'] - meta['weight_compile_ns']
    return Prepared(kind, method, task, problem, basis, (mode, payload, mapping), checker, meta)


def verify_solution(problem, rhs, answer):
    if answer is None:
        return
    x = column(answer)
    if problem['kind'] == 'mlc':
        ax = vector(Z(problem['A']) * x)
        assert all((v-b) % m == 0 for v, b, m in zip(ax, rhs, problem['moduli']))
    else:
        for m, r in zip(problem['matrices'], rhs):
            z = Z(m).inv() * column([a-b for a, b in zip(answer, r)])
            assert all(v.denominator == 1 for v in z.entries())


def oracle(problem):
    """Independent full augmented-SNF oracle; small correctness tests only."""
    from sympy import ZZ
    from sympy.polys.matrices import DomainMatrix
    from sympy.polys.matrices.normalforms import smith_normal_decomp
    c = augmented(problem['A'], problem['moduli']) if problem['kind'] == 'mlc' else md_global(problem['matrices'])
    s, u, v = smith_normal_decomp(DomainMatrix.from_list(rows(c), ZZ))
    s, u, v = [m.to_Matrix() for m in (s, u, v)]
    r = c.nrows()
    n = len(problem['A'][0]) if problem['kind'] == 'mlc' else len(problem['matrices'][0])
    b = Z([[int(v[i, j]) for j in range(r, c.ncols())] for i in range(n)])
    def accepts(rhs):
        flat = list(rhs) if problem['kind'] == 'mlc' else [a for rr in rhs for a in rr]
        return all(sum(int(u[i,j])*flat[j] for j in range(r)) % int(s[i,i]) == 0 for i in range(r))
    return b, accepts


def equal_lattices(a, b):
    a, b = Z(a), Z(b)
    for x in (a.inv()*b, b.inv()*a):
        assert all(v.denominator == 1 for v in x.entries())
    return True
