"""Three-route benchmark for unstructured, compatible MD-CRT inputs."""

import argparse
import csv
from hashlib import sha256
from itertools import product
import json
from math import lcm
from pathlib import Path
import platform
from random import Random
import subprocess
import sys
from time import perf_counter_ns

import flint

HERE = Path(__file__).resolve().parent
ALGORITHMS = HERE.parents[1] / 'algorithms'
sys.path.insert(0, str(ALGORITHMS))

from mdcrt_adjugate_global import AdjugateReduction, AuxiliaryHNF, CongruenceHNF
from mlc_mdcrt_exact.bidirectional import Compiled, Z, column, decomposition, rows, triangular, vector

flint.ctx.threads = 1
# Exact period bases in the large cases can exceed Python's default JSON digit limit.
sys.set_int_max_str_digits(100_000)
SIZES = (16, 32, 64)
QUERIES = (1, 100)
METHODS = ('cascade_hnf', 'adj_aux', 'adj_congruence')
BASE_SEED = 20261125
ENTRY_LIMIT = 127
WITNESS_LIMIT = 32767
WORKER_TIMEOUT = 7200
LABELS = dict(zip(METHODS, ('矩阵级联', '伴随＋辅助变量 HNF', '伴随＋同余行格 HNF')))


def stable_json(value):
    return json.dumps(value, separators=(',', ':'), ensure_ascii=False).encode('utf-8')


def digest(value):
    return sha256(stable_json(value)).hexdigest()


def integer_rows_digest(matrix):
    checksum = sha256()
    for row in matrix:
        checksum.update(len(row).to_bytes(8, 'little'))
        for item in row:
            item = int(item)
            data = abs(item).to_bytes((abs(item).bit_length() + 7) // 8, 'little')
            checksum.update(bytes([item < 0]))
            checksum.update(len(data).to_bytes(8, 'little'))
            checksum.update(data)
    return checksum.hexdigest()


def max_bits(values):
    return max((abs(int(value)).bit_length() for value in values), default=0)


def make_input(n, k):
    rng = Random(f'{BASE_SEED}:matrix:{n}')
    matrices = []
    determinants = []
    rejected = 0
    while len(matrices) < k:
        matrix = [[rng.randint(-ENTRY_LIMIT, ENTRY_LIMIT) for _ in range(n)] for _ in range(n)]
        det = int(Z(matrix).det())
        if det:
            matrices.append(matrix)
            determinants.append(abs(det))
        else:
            rejected += 1

    rng = Random(f'{BASE_SEED}:witness:{n}')
    witnesses = [tuple(rng.randint(-WITNESS_LIMIT, WITNESS_LIMIT) for _ in range(n))
                 for _ in range(QUERIES[-1])]
    native = list(map(Z, matrices))
    columns = list(map(column, witnesses))
    queries = [[] for _ in witnesses]
    duals = []
    for matrix, determinant in zip(native, determinants):
        inverse = matrix.inv()
        inverse_map = Compiled.make(inverse)
        # Rows of inv(M) generate the dual of the column lattice M Z^n.
        duals.append((Z([[int(v * determinant) for v in row]
                         for row in inverse.tolist()]), determinant))
        for witness, query in zip(columns, queries):
            coordinates = vector(inverse_map.numerator * witness)
            quotient = column(a // d for a, d in zip(coordinates, inverse_map.denominator))
            query.append(vector(witness - matrix * quotient))
    return matrices, queries, witnesses, duals, determinants, rejected


def reference_basis(duals):
    """Independent homogeneous construction via balanced sums of dual lattices."""
    layer = list(duals)
    while len(layer) > 1:
        merged = []
        for i in range(0, len(layer) - 1, 2):
            left, d_left = layer[i]
            right, d_right = layer[i + 1]
            common = lcm(d_left, d_right)
            generators = rows(left * (common // d_left)) + rows(right * (common // d_right))
            hermite = Z(generators).hnf()
            merged.append((Z(rows(hermite)[:left.nrows()]), common))
        layer = merged + (layer[-1:] if len(layer) % 2 else [])
    basis, denominator = layer[0]
    return basis.inv() * denominator


def require_integer(matrix, error):
    if any(value.denominator != 1 for value in matrix.entries()):
        raise ArithmeticError(error)


class MatrixCascade:
    """Cascade the same local merge rules while normalizing intermediate column bases."""

    def __init__(self, matrices):
        native = list(map(Z, matrices))
        n = native[0].nrows()
        basis = native[0].transpose().hnf().transpose()
        self.stages = []
        for matrix in native[1:]:
            augmented = Z([br + [-value for value in mr]
                           for br, mr in zip(rows(basis), rows(matrix))])
            triangle, lift, kernel, _ = decomposition(augmented, n)
            lift = basis * lift
            basis = (basis * kernel).transpose().hnf().transpose()
            self.stages.append((rows(triangle), lift, rows(basis)))
        self.B = tuple(map(tuple, rows(basis)))

    def solve(self, residues):
        answer = residues[0]
        for residue, (triangle, lift, basis) in zip(residues[1:], self.stages):
            point = triangular(triangle, [a - b for a, b in zip(residue, answer)])
            if point is None:
                return None
            pending = [a + b for a, b in zip(answer, vector(lift * column(point)))]
            # Subtract integer combinations of columns of the lower-triangular
            # period basis, keeping an equivalent representative for the next stage.
            for j, row in enumerate(basis):
                quotient, pending[j] = divmod(pending[j], row[j])
                if quotient:
                    for i in range(j + 1, len(pending)):
                        pending[i] -= quotient * basis[i][j]
            answer = tuple(pending)
        return answer


def build(method, matrices):
    start = perf_counter_ns()
    if method == 'cascade_hnf':
        solver = MatrixCascade(matrices)
        end = perf_counter_ns()
        return solver, end - start, 0, {'merge_steps': len(matrices) - 1}
    reduction = AdjugateReduction(matrices)
    converted = perf_counter_ns()
    solver = (AuxiliaryHNF if method == 'adj_aux' else CongruenceHNF)(reduction)
    end = perf_counter_ns()
    return solver, end - start, converted - start, {
        'retained_rows': len(reduction.moduli),
        'rows_before': reduction.n * reduction.count,
        'single_row_blocks': sum(len(moduli) == 1 for _, moduli in reduction.blocks),
        'global_factor': reduction.factor,
        'reduced_modulus_bits': max_bits(reduction.moduli),
        'work_shape': solver.work_shape,
        'unified_modulus_bits': solver.common.bit_length() if method == 'adj_congruence' else None,
    }


def append(path, record):
    with path.open('a', encoding='utf-8') as log:
        log.write(json.dumps(record, ensure_ascii=False) + '\n')


def generate(n, k, directory):
    start = perf_counter_ns()
    matrices, queries, witnesses, duals, determinants, rejected = make_input(n, k)
    generated = perf_counter_ns()
    reference = reference_basis(duals)
    require_integer(reference, 'noninteger reference basis')
    basis = rows(reference)
    # The original moduli independently confirm every generated witness/residue pair.
    for i, (numerators, denominator) in enumerate(duals):
        differences = Z([tuple(a - b for a, b in zip(witness, query[i]))
                         for witness, query in zip(witnesses, queries)]).transpose()
        if any(v % denominator for v in (numerators * differences).entries()):
            raise ArithmeticError('witness does not satisfy the original blocks')
    verified = perf_counter_ns()
    data = {'matrices': matrices, 'queries': queries, 'witnesses': witnesses, 'reference': basis}
    (directory / 'input.json').write_bytes(stable_json(data))
    info = {'n': n, 'k': k, 'seed': BASE_SEED, 'input_sha256': digest(data),
            'reference_sha256': digest(basis), 'singular_rejections': rejected,
            'matrix_entry_bits': max_bits(v for m in matrices for row in m for v in row),
            'determinant_max_bits': max_bits(determinants),
            'residue_entry_bits': max_bits(v for query in queries for residue in query for v in residue),
            'reference_basis_bits': max_bits(v for row in basis for v in row),
            'generation_ns': generated - start, 'reference_and_witness_check_ns': verified - generated}
    (directory / 'input_metadata.json').write_bytes(stable_json(info))


def trial(n, k, method, directory):
    data = json.loads((directory / 'input.json').read_bytes())
    info = json.loads((directory / 'input_metadata.json').read_bytes())
    if (info['n'], info['k']) != (n, k) or digest(data) != info['input_sha256']:
        raise ValueError('input record does not match its metadata')
    reference = Z(data['reference'])
    reference_inverse = reference.inv()
    log = directory / f'{method}.jsonl'
    append(log, {'phase': 'preparing'})
    solver, setup, conversion, details = build(method, data['matrices'])
    append(log, {'phase': 'prepared', 'setup_ns': setup, **details})
    query_ns = 0
    check_ns = 0
    answers = []
    for q, (rhs, witness) in enumerate(zip(data['queries'], data['witnesses']), 1):
        started = perf_counter_ns()
        answer = solver.solve(rhs)
        query_ns += perf_counter_ns() - started
        checked = perf_counter_ns()
        if answer is None:
            raise ArithmeticError(f'query {q} reported inconsistent despite certified witness')
        if q == 1:
            # Equality of column lattices needs both inclusions.
            actual = Z(solver.B)
            require_integer(reference_inverse * actual, 'basis too large')
            require_integer(actual.inv() * reference, 'basis too small')
        require_integer(reference_inverse * column(a - b for a, b in zip(answer, witness)),
                        f'query {q} differs from witness outside the period lattice')
        check_ns += perf_counter_ns() - checked
        answers.append(answer)
        if q in QUERIES:
            record = dict(n=n, k=k, method=method, q=q, setup_ns=setup, conversion_ns=conversion,
                          query_ns=query_ns, query_hnf_ns=getattr(solver, 'hnf_ns', None),
                          total_ns=setup + query_ns, verification_ns=check_ns,
                          answer_sha256=integer_rows_digest(answers), input_sha256=info['input_sha256'],
                          reference_sha256=info['reference_sha256'],
                          basis_sha256=integer_rows_digest(solver.B), **details)
            append(log, {'phase': 'verified', **record})


def run_worker(directory, name, arguments, timeout):
    start = perf_counter_ns()
    print(f'{directory.name}: {name}', flush=True)
    with (directory / f'{name}.log').open('wb') as log:
        try:
            proc = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()), *arguments],
                                  stdout=log, stderr=log, timeout=timeout)
            status = 'ok' if proc.returncode == 0 else 'error'
        except subprocess.TimeoutExpired:
            status = 'timeout'
    elapsed = perf_counter_ns() - start
    append(directory.parent / 'jobs.jsonl',
           {'case': directory.name, 'worker': name, 'status': status, 'wall_ns': elapsed})
    print(f'{directory.name}: {name} {status}, {elapsed / 1e9:.2f} s', flush=True)
    if status == 'error':
        raise RuntimeError(f'{directory.name}: {name} failed, see {name}.log')
    return status


def write_csv(path, records):
    fields = list(dict.fromkeys(key for record in records for key in record))
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields)
        writer.writeheader()
        writer.writerows(records)


def export(directory, sizes):
    measurements, inputs = [], []
    for n, k in product(sizes, repeat=2):
        case = directory / f'n{n}_k{k}'
        info = case / 'input_metadata.json'
        if info.exists():
            inputs.append(json.loads(info.read_bytes()))
        for method in METHODS:
            log = case / f'{method}.jsonl'
            if log.exists():
                for line in log.read_bytes().splitlines():
                    event = json.loads(line)
                    if event.pop('phase') == 'verified':
                        measurements.append(event)
    write_csv(directory / 'measurements.csv', measurements)
    write_csv(directory / 'inputs.csv', inputs)
    manifest = json.loads((directory / 'manifest.json').read_bytes())
    index = {(r['n'], r['k'], r['q'], r['method']): r for r in measurements}
    report = ['# 随机可解 MD-CRT：三路线实测', '',
              '两档时间均包含一次预处理和该档的全部查询。每次求解返回一个特解，共享完整周期格基。',
              '每组同一随机输入按三个独立进程串行运行，每条路线计时一次。单位为毫秒。', '',
              f"实测后端为 python-flint {manifest['python_flint']}，FLINT {manifest['flint']}，单线程。", '']
    for q in QUERIES:
        report += [f'## Q = {q}', '',
                   '| n | k | 矩阵级联 | 伴随＋辅助变量 HNF | 伴随＋同余行格 HNF | 最快 |',
                   '|---:|---:|---:|---:|---:|---|']
        for n, k in product(sizes, repeat=2):
            group = [index.get((n, k, q, method)) for method in METHODS]
            cells = ['—' if row is None else f"{row['total_ns']/1e6:,.3f}" for row in group]
            fastest = LABELS[min(group, key=lambda r: r['total_ns'])['method']] if all(group) else '缺少配对'
            report.append(f'| {n} | {k} | ' + ' | '.join(cells) + f' | {fastest} |')
        report.append('')
    report += ['## 核验覆盖', '', f'实测并通过精确核验的计时点共 {len(measurements)} 个，'
               f'目标 {len(sizes)**2 * len(QUERIES) * len(METHODS)} 个。']
    (directory / 'REPORT.md').write_text('\n'.join(report) + '\n', encoding='utf-8')


def run(args):
    root = args.directory.resolve()
    root.mkdir(parents=True, exist_ok=False)
    sizes = (4,) if args.smoke else SIZES
    sources = [Path(__file__).resolve(), ALGORITHMS / 'mdcrt_adjugate_global.py',
               ALGORITHMS / 'mdcrt_adjugate_dual.py',
               ALGORITHMS / 'mlc_mdcrt_exact' / 'bidirectional.py']
    manifest = {'sizes': sizes, 'queries': QUERIES, 'methods': METHODS, 'base_seed': BASE_SEED,
                'entry_limit': ENTRY_LIMIT, 'witness_limit': WITNESS_LIMIT,
                'repetitions': 1, 'timeout_s': args.timeout, 'threads': flint.ctx.threads,
                'python': platform.python_version(), 'platform': platform.platform(),
                'processor': platform.processor(), 'python_flint': flint.__version__,
                'flint': flint.__FLINT_VERSION__, 'sha256_by_source_file':
                {path.name: sha256(path.read_bytes()).hexdigest() for path in sources}}
    (root / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    for order, (n, k) in enumerate(product(sizes, repeat=2)):
        case = root / f'n{n}_k{k}'
        case.mkdir()
        dimensions = ['--n', str(n), '--k', str(k)]
        if run_worker(case, 'generate', ['generate', str(case), *dimensions], args.timeout) == 'ok':
            methods = METHODS[order % len(METHODS):] + METHODS[:order % len(METHODS)]
            for method in methods:
                run_worker(case, method, ['trial', str(case), *dimensions, '--method', method], args.timeout)
                export(root, sizes)
    export(root, sizes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    execute = sub.add_parser('run')
    execute.add_argument('directory', type=Path)
    execute.add_argument('--timeout', type=int, default=WORKER_TIMEOUT)
    execute.add_argument('--smoke', action='store_true')
    for name in ('generate', 'trial'):
        worker = sub.add_parser(name)
        worker.add_argument('directory', type=Path)
        worker.add_argument('--n', type=int, required=True)
        worker.add_argument('--k', type=int, required=True)
        if name == 'trial':
            worker.add_argument('--method', choices=METHODS, required=True)
    args = parser.parse_args()
    if args.command == 'run':
        run(args)
    elif args.command == 'generate':
        generate(args.n, args.k, args.directory)
    else:
        trial(args.n, args.k, args.method, args.directory)


if __name__ == '__main__':
    main()
