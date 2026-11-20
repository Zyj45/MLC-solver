"""Compare exact full-solution MD-CRT routes on identical compatible queries."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import subprocess
import sys
from time import perf_counter_ns

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'algorithms'))

from mdcrt_adjugate_dual import AdjugateDualMDCRT
from mlc_mdcrt_exact.bidirectional import Z, column, prepare, vector
from mlc_mdcrt_exact.bidirectional_inputs import make_matrices


METHODS = ('global_weight', 'cascade_weight', 'adjugate_dual')
GRIDS = (('cyclic_shared', 32, (2, 8, 16)), ('dense_shared', 8, (2, 8)),
         ('dense_primitive', 16, (2, 8)))
QUERIES = (1, 10, 100, 1000)


def trial_case(family, n, k, seed, method):
    base_family = 'dense_shared' if family == 'dense_primitive' else family
    matrices = make_matrices(base_family, n, 16, 16, seed)[:k]
    if family == 'dense_primitive':
        matrices = [[[value // 2 for value in row] for row in matrix] for matrix in matrices]
    native = [Z(matrix) for matrix in matrices]
    inverses = [matrix.inv() for matrix in native]
    rng = random.Random(20261111 + seed)
    queries = []
    for _ in range(QUERIES[-1]):
        x = column([rng.randint(-65535, 65535) for _ in range(n)])
        residues = []
        for matrix, inverse in zip(native, inverses):
            coordinates = inverse * x
            quotient = column([int(v.numerator) // int(v.denominator) for v in coordinates.entries()])
            residues.append(vector(x - matrix * quotient))
        queries.append(residues)

    def make_solver():
        if method == 'adjugate_dual':
            return AdjugateDualMDCRT(matrices)
        return prepare({'kind': 'md', 'matrices': matrices}, method)

    assert make_solver().solve(queries[0]) is not None
    measurements, results = [], []
    for repetition in range(1, 4):
        start = perf_counter_ns()
        solver = make_solver()
        setup = perf_counter_ns() - start
        totals, answers = {}, []
        start = perf_counter_ns()
        for q, residues in enumerate(queries, 1):
            answers.append(solver.solve(residues))
            if q in QUERIES:
                totals[f'total_q{q}_ns'] = setup + perf_counter_ns() - start
        measurements.append({'trial': repetition, 'prepare_ns': setup, **totals})
        results.append((solver.B, answers))

    # Validation happens outside the measurement window. Compare canonical
    # column lattices across routes and check all returned solutions directly.
    digests = set()
    for basis, answers in results:
        basis = Z(basis)
        assert basis.det()
        canonical = basis.transpose().hnf()
        digests.add(hashlib.sha256('|'.join(hex(int(v)) for v in canonical.entries()).encode()).hexdigest())
        for inverse in inverses:
            assert all(v.denominator == 1 for v in (inverse * basis).entries())
        for residues, answer in zip(queries, answers):
            assert answer is not None
            for residue, inverse in zip(residues, inverses):
                diff = column([x - r for x, r in zip(answer, residue)])
                assert all(value.denominator == 1 for value in (inverse * diff).entries())
    assert len(digests) == 1
    return {'trials': measurements, 'basis_hash': digests.pop(),
            'shared_factor': getattr(solver, 'factor', ''),
            'scalar_rows': getattr(solver, 'row_count', ''),
            'common_modulus_bits': getattr(solver, 'modulus', 0).bit_length() or ''}


def benchmark(path, timeout):
    fields = ['family', 'n', 'k', 'seed', 'method', 'status', 'trial', 'prepare_ns',
              *(f'total_q{q}_ns' for q in QUERIES), 'basis_hash', 'shared_factor',
              'scalar_rows', 'common_modulus_bits']
    with path.open('x', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        for family, n, counts in GRIDS:
            for k in counts:
                for seed in (1, 2, 3):
                    hashes = set()
                    offset = seed % len(METHODS)
                    for method in METHODS[offset:] + METHODS[:offset]:
                        base = {'family': family, 'n': n, 'k': k, 'seed': seed, 'method': method}
                        command = [sys.executable, '-B', __file__, '--case', family, str(n), str(k), str(seed), method]
                        try:
                            process = subprocess.run(command, capture_output=True, text=True, timeout=timeout, check=True)
                        except subprocess.TimeoutExpired:
                            writer.writerow({**base, 'status': 'worker_timeout'})
                            print(f'Timeout {family} n={n} k={k} seed={seed} {method}', flush=True)
                        except subprocess.CalledProcessError as error:
                            raise RuntimeError(error.stderr) from error
                        else:
                            result = json.loads(process.stdout)
                            hashes.add(result['basis_hash'])
                            for trial in result.pop('trials'):
                                writer.writerow({**base, 'status': 'verified', **result, **trial})
                        file.flush()
                    assert len(hashes) <= 1, (family, n, k, seed, hashes)
                print(f'Finished {family} n={n} k={k}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--case', nargs=5)
    group.add_argument('--out', type=Path)
    parser.add_argument('--timeout', type=float, default=45)
    args = parser.parse_args()
    if args.case:
        family, n, k, seed, method = args.case
        print(json.dumps(trial_case(family, int(n), int(k), int(seed), method)))
    else:
        benchmark(args.out, args.timeout)
