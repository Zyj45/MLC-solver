"""Deterministic original-input generators and oracle-only certificates."""
from hashlib import sha256
from math import gcd, lcm
import json
import random
from .bidirectional import (Z, Q, eye, rows, vector, column, decomposition,
                            augmented, md_global, Compiled, bits)

BASE_SEED = 20261105
MLC_FAMILIES = ('dense_shared', 'repeated_nested', 'sparse_rankdef')
MD_FAMILIES = ('dense_shared', 'cyclic_shared', 'nested', 'diagonal_coprime', 'rebased')


def seeded(*parts):
    return random.Random(int.from_bytes(sha256('|'.join(map(str,(BASE_SEED,)+parts)).encode()).digest()[:8], 'big'))


def primes(w, count, rng):
    from flint import fmpz
    low, high = 1 << (w-1), 1 << w
    p, found = rng.randrange(low, high) | 1, []
    while len(found) < count:
        if p >= high:
            p = low | 1
        if fmpz(p).is_prime() and p not in found:
            found.append(p)
        p += 2
    return found


def remainders(matrices, x):
    out = []
    for m in matrices:
        inv = Z(m).inv()
        q = [int(v.numerator)//int(v.denominator) for v in (inv*column(x)).entries()]
        shift = vector(Z(m)*column(q))
        out.append([a-b for a,b in zip(x,shift)])
    return out


def make_matrices(family, n, k, w, seed):
    base_family = 'dense_shared' if family == 'rebased' else family
    rng = seeded('md',base_family,n,k,w,seed)
    ms = []
    if family == 'diagonal_coprime':
        return [rows(eye(n)*p) for p in primes(w,k,rng)]
    if family == 'nested':
        u = rows(eye(n))
        for i in range(n-1):
            u[i][i+1] = rng.choice((-1,1))
        for t in range(k):
            ds = [1 << (1 + (t*(w-2))//max(1,k-1)) for _ in range(n)]
            ms.append(rows(Z(u)*Z([[ds[i] if i==j else 0 for j in range(n)] for i in range(n)])))
        return ms
    for t in range(k):
        if base_family == 'dense_shared':
            bound = max(2, ((1<<w)-1)//2)
            while True:
                m = Z([[2*rng.randint(-bound,bound) for j in range(n)] for i in range(n)])
                if m.det():
                    break
        else:
            # One nonunit Smith factor before the shared factor 2.
            u, v = rows(eye(n)), rows(eye(n))
            for j in range(n-1):
                u[j][j+1] = rng.choice((-1,1))
                v[j+1][j] = rng.choice((-1,1))
            ds = [1]*(n-1) + [(1 << max(1,w-3)) + 1]
            m = 2*Z(u)*Z([[ds[i] if i==j else 0 for j in range(n)] for i in range(n)])*Z(v)
        ms.append(rows(m))
    if family == 'rebased':
        # Same lattices and numeric observations; deliberate extra 8-bit basis shear.
        rr = seeded('rebase',n,k,w,seed)
        out = []
        for m in ms:
            v = rows(eye(n))
            i,j = rr.sample(range(n),2)
            v[i][j] = rr.choice((-256,256))
            out.append(rows(Z(m)*Z(v)))
        ms = out
    return ms


def make_case(kind, family, n, size, w, seed, q=1000):
    rng = seeded(kind,family,n,size,w,seed)
    limit = (1<<w)-1
    if kind == 'mlc':
        r = size
        a = [[rng.randint(-limit,limit) for _ in range(n)] for _ in range(r)]
        if family == 'dense_shared':
            mods = [6*rng.randint(max(1,(1<<(w-1))//6),max(1,limit//6)) for _ in range(r)]
        else:
            mods = [1 << (w-1-(i % max(1,w//2))) for i in range(r)]
            if family == 'repeated_nested':
                for i in range(1,r,3):
                    a[i] = a[i-1][:]
            else:
                d = max(1,n//2)
                a = [[v if j < d and rng.random() < .5 else 0 for j,v in enumerate(row)] for row in a]
        problem = {'kind':kind,'A':a,'moduli':mods}
        c = augmented(a,mods)
    else:
        matrices = make_matrices(family,n,size,w,seed)
        problem = {'kind':kind,'matrices':matrices}
        c = md_global(matrices)
    # Independent original embedding only; never passed into prepare().
    g, part, basis, tb = decomposition(c,n)
    invg = Compiled.make(g.inv())
    bad = next((i for i,d in enumerate(invg.denominator) if d>1),None)
    certificate = None
    perturb = None
    if bad is not None:
        d = invg.denominator[bad]
        num = rows(invg.numerator)[bad]
        perturb = next(j for j,v in enumerate(num) if v % d)
        assert all(int(v) % d == 0 for v in (Z([num])*c).entries())
        certificate = {'row':num,'modulus':d,'perturb_coordinate':perturb}
    panels = {'compatible':[], 'incompatible':[]}
    witness = []
    rhs_rng = seeded('rhs',kind, 'dense_shared' if family=='rebased' else family,n,size,w,seed)
    residue_matrices = (make_matrices('dense_shared',n,size,w,seed) if family=='rebased' else problem.get('matrices'))
    # Compile exact floor maps once; input generation must not invert per query.
    residue_maps = []
    if kind == 'md':
        for m in residue_matrices:
            residue_maps.append((Z(m),Compiled.make(Z(m).inv())))
    for _ in range(q):
        x = [rhs_rng.randint(-limit,limit) for _ in range(n)]
        if kind == 'mlc':
            rhs = [v % m for v,m in zip(vector(Z(a)*column(x)),mods)]
        else:
            rhs = []
            for m, inv in residue_maps:
                numerators = vector(inv.numerator*column(x))
                quotient = [v//d for v,d in zip(numerators,inv.denominator)]
                shift = vector(m*column(quotient))
                rhs.append([a-b for a,b in zip(x,shift)])
        panels['compatible'].append(rhs)
        witness.append(x)
        if certificate:
            flat = list(rhs) if kind=='mlc' else [v for row in rhs for v in row]
            flat[perturb] += 1
            assert sum(a*b for a,b in zip(certificate['row'],flat)) % certificate['modulus']
            invalid = flat if kind=='mlc' else [flat[i:i+n] for i in range(0,len(flat),n)]
            panels['incompatible'].append(invalid)
    case_id = f'{kind}_{family}_n{n}_s{size}_w{w}_seed{seed}'
    return {'id':case_id,'kind':kind,'family':family,'n':n,'size':size,'width':w,'seed':seed,
            'problem':problem,'panels':panels,'oracle_basis':rows(basis),'witness':witness,
            'incompatible_certificate':certificate,
            'input_bits':bits(c),'native_shape':[c.nrows(),c.ncols()],
            'source_seed':BASE_SEED,'representation_pair':('dense_shared' if family=='rebased' else None)}


def descriptors():
    """Frozen breadth-first range: skeleton, bit extensions, shape extensions."""
    groups = []
    for phase in ('skeleton','bits','shape'):
        for n in (4,8,16,32):
            kinds = []
            for kind, families in (('mlc',MLC_FAMILIES),('md',MD_FAMILIES)):
                sizes = [n] if kind=='mlc' else [2]
                widths = [16]
                if phase=='bits': widths=[8,32]
                if phase=='shape': sizes=[n//2,2*n] if kind=='mlc' else [4]
                kinds.append([(kind,f,n,s,w,0) for f in families for s in sizes for w in widths])
            for i in range(max(map(len,kinds))):
                for cases in kinds:
                    if i<len(cases): groups.append(cases[i])
    return groups
