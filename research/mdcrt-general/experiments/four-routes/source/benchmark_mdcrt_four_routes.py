"""Four-route exact MD-CRT experiment, with independently verified checkpoints."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
from itertools import product
import json
from math import gcd, lcm
from pathlib import Path
import pickle
import platform
from random import Random
from statistics import median
import subprocess
import sys
from time import monotonic, perf_counter_ns
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / 'algorithms'))

from mdcrt_adjugate_global import AdjugateReduction, AuxiliaryHNF, CongruenceHNF
from mlc_mdcrt_exact.bidirectional import Compiled, Z, column, prepare, rows, vector
from mlc_mdcrt_exact.bidirectional_inputs import make_matrices
import flint

flint.ctx.threads = 1
METHODS = ('global_tri', 'cascade_tri', 'adj_aux', 'adj_congruence')
FAMILIES = ('cyclic_shared', 'dense_shared', 'dense_primitive', 'diagonal_coprime')
SIZES = (8, 16, 32)
QS = (1, 10, 100, 1000)
BASE_SEED = 20261105
NAMES = dict(zip(METHODS, ('原生堆叠', '原生级联', '伴随＋辅助变量', '伴随＋同余行格')))
LABELS = ('循环共享因子', '稠密共享因子', '稠密无公共标量因子', '对角互素')
DESCRIPTIONS = (
    '每个模矩阵带公共因子 2。提取 2 后，各块的整数商群为循环群，测试单块压缩的收益。',
    '用有符号的偶数元素生成非奇异稠密块，测试公共因子与块内约简叠加的效果。',
    '将稠密偶数块除以 2，并核验最短前缀的全局元素 gcd 为 1。',
    '每块采用不同的 16 位素数乘单位矩阵，保留各坐标的独立约束，测试难以删行的情况。',
)


def dump(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    with gzip.open(temporary, 'wb', compresslevel=1) as stream:
        pickle.dump(value, stream, protocol=5)
    temporary.replace(path)
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    # Only read local files produced by this experiment.
    with gzip.open(path, 'rb') as stream:
        return pickle.load(stream)


def emit(path, value):
    with path.open('a', encoding='utf-8') as stream:
        stream.write(json.dumps(value, ensure_ascii=False) + '\n')


def events(path):
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_bytes().splitlines(keepends=True)
            if line.endswith(b'\n')]


def bits(values):
    return max((abs(int(value)).bit_length() for value in values), default=0)


def inverse_maps(matrices):
    return [Compiled.make(Z(matrix).inv()) for matrix in matrices]


def reference_basis(inverses, n):
    """Homogeneous oracle from all inverse rows, with no adjugate row selection."""
    common = lcm(*(d for inverse in inverses for d in inverse.denominator))
    generators = [[common if i == j else 0 for j in range(n)] for i in range(n)]
    for inverse in inverses:
        for row, d in zip(rows(inverse.numerator), inverse.denominator):
            if d > 1:
                generators.append([(value % d) * (common // d) for value in row])
    hnf = Z(generators).hnf()
    basis = Z([[hnf[i, j] for j in range(n)] for i in range(n)]).inv() * common
    assert all(value.denominator == 1 for value in basis.entries())
    return rows(basis)


def generate(directory, family, n, seed):
    source = 'dense_shared' if family == 'dense_primitive' else family
    matrices = make_matrices(source, n, 32, 16, seed)
    if family == 'dense_primitive':
        matrices = [[[value // 2 for value in row] for row in matrix] for matrix in matrices]
        assert gcd(*(v for matrix in matrices[:8] for row in matrix for v in row)) == 1
    inverses = inverse_maps(matrices)
    native = list(map(Z, matrices))
    rng = Random(f'{BASE_SEED}|{family}|{n}|{seed}')
    witnesses = [tuple(rng.randint(-65535, 65535) for _ in range(n)) for _ in range(QS[-1])]
    queries = []
    for witness in witnesses:
        x = column(witness)
        residues = []
        for matrix, inverse in zip(native, inverses):
            coordinates = vector(inverse.numerator * x)
            quotient = column(value // d for value, d in zip(coordinates, inverse.denominator))
            residues.append(vector(x - matrix * quotient))
        queries.append(tuple(residues))
    refs = {k: reference_basis(inverses[:k], n) for k in SIZES}
    payload = {'matrices': matrices, 'queries': queries, 'witnesses': witnesses, 'references': refs,
               'determinant_bits': [abs(int(m.det())).bit_length() for m in native]}
    digest = dump(directory / 'input.pkl.gz', payload)
    emit(directory / 'events.jsonl', {'event': 'generated', 'sha256': digest,
                                      'input_bits': bits(v for m in matrices for row in m for v in row),
                                      'determinant_bits': payload['determinant_bits']})


def build(method, matrices):
    start = perf_counter_ns()
    if method.startswith('adj_'):
        reduction = AdjugateReduction(matrices)
        converted = perf_counter_ns()
        solver = (AuxiliaryHNF if method == 'adj_aux' else CongruenceHNF)(reduction)
        end = perf_counter_ns()
        meta = {'conversion_ns': converted - start, 'decomposition_setup_ns': end - converted,
                'common_factor': reduction.factor, 'scalar_rows': len(reduction.moduli),
                'rows_before': reduction.n * reduction.count,
                'modulus_bits': bits(reduction.moduli),
                'coefficient_bits': bits(a for row in reduction.coefficients for a in row),
                'unified_modulus_bits': solver.common.bit_length() if method == 'adj_congruence' else None,
                'work_shapes': [solver.work_shape] if solver.work_shape is not None else []}
    else:
        solver = prepare({'kind': 'md', 'matrices': matrices}, method)
        end = perf_counter_ns()
        meta = {'conversion_ns': 0, 'validation_ns': solver.metadata['validation_ns'],
                'decomposition_setup_ns': solver.metadata['preparation_excluding_weight_ns'],
                'work_shapes': [list(reversed(shape)) for shape in solver.metadata['work_shapes']]}
    return solver, end - start, meta


def check_basis(basis, reference):
    actual, expected = Z(basis), Z(reference)
    for source, destination in ((actual, expected), (expected, actual)):
        assert all(v.denominator == 1 for v in (source.inv() * destination).entries())


def check_solution(solution, rhs, inverses):
    assert solution is not None
    for residue, inverse in zip(rhs, inverses):
        assert inverse.apply([x - r for x, r in zip(solution, residue)]) is not None


def trial(directory, input_path, k, method, target):
    log = directory / 'events.jsonl'
    data = load(input_path)
    matrices = data['matrices'][:k]
    queries = [batch[:k] for batch in data['queries'][:target]]
    inverses = inverse_maps(matrices)
    emit(log, {'event': 'preparation'})
    solver, preparation, meta = build(method, matrices)
    emit(log, {'event': 'prepared', 'prepare_ns': preparation, **meta,
               'input_bits': bits(v for m in matrices for row in m for v in row),
               'determinant_bits': max(data['determinant_bits'][:k])})
    query_ns = 0
    answers = []
    first = 0
    for q, rhs in enumerate(queries, 1):
        started = perf_counter_ns()
        solution = solver.solve(rhs)
        query_ns += perf_counter_ns() - started
        answers.append(solution)
        if q not in QS:
            continue
        emit(log, {'event': 'verify', 'q': q})
        if q == 1:
            basis = solver.B
            check_basis(basis, data['references'][k])
            basis_hash = hashlib.sha256(pickle.dumps(data['references'][k], protocol=5)).hexdigest()
        for answer, observation in zip(answers, queries[first:q]):
            check_solution(answer, observation, inverses)
        emit(log, {'event': 'evidence', 'q': q})
        digest = dump(directory / f'q{q}.pkl.gz',
                      {'basis': basis if q == 1 else None, 'answers': answers, 'start': first})
        emit(log, {'event': 'checkpoint', 'q': q, 'prepare_ns': preparation,
                   'query_ns': query_ns, 'total_ns': preparation + query_ns,
                   'query_hnf_ns': getattr(solver, 'hnf_ns', None),
                   'basis_hash': basis_hash, 'evidence_sha256': digest,
                   'basis_bits': bits(v for row in basis for v in row)})
        emit(log, {'event': 'query'})
        first, answers = q, []
    emit(log, {'event': 'complete'})


def self_test():
    base = [[[2, 0], [2, 6]], [[2, 0], [4, 10]]]
    cases = [
        (base, [(1, 5), (1, 1)], True),
        (base, [(1, 5), (0, 1)], False),
        ([[[4]], [[4]], [[6]]], [(0,), (2,), (0,)], False),
        ([[[2, 0], [0, -3]]], [(1, 2)], True),
        ([[[2, 0], [0, 2]], [[0, 2], [-2, 0]]], [(1, 2), (3, 4)], True),
    ]
    for matrices, query, expected in cases:
        inverses = inverse_maps(matrices)
        reference = reference_basis(inverses, len(matrices[0]))
        for method in METHODS:
            solver, _, _ = build(method, matrices)
            answer = solver.solve(query)
            assert (answer is not None) == expected, (method, matrices, query)
            if expected:
                check_basis(solver.B, reference)
                check_solution(answer, query, inverses)
                next_query = [tuple(value + 4 for value in r) for r in query]
                check_solution(solver.solve(next_query), next_query, inverses)
                check_basis(solver.B, reference)
    assert len(AdjugateReduction([[[2, 0], [0, -3]]]).moduli) == 2
    print(f'{len(cases)*len(METHODS)} route/edge-case checks passed')


def serial_worker(directory, arguments, cap, cutoff):
    directory.mkdir(parents=True, exist_ok=True)
    started = monotonic()
    with (directory / 'worker.log').open('wb') as log:
        try:
            process = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()), *arguments],
                                     stdout=log, stderr=log,
                                     timeout=max(.001, min(cap, cutoff - started)))
        except subprocess.TimeoutExpired:
            status = 'budget_cutoff' if monotonic() >= cutoff else 'timeout'
        else:
            status = 'ok' if process.returncode == 0 else 'error'
    return status, monotonic() - started


def grid(smoke):
    if smoke:
        return [(FAMILIES[0], 8, 8), (FAMILIES[-1], 8, 8)]
    return [(family, n, k) for n, k in product(SIZES, repeat=2) for family in FAMILIES]


def schedule(smoke, started, seconds):
    configurations = grid(smoke)
    end = started + seconds
    yield 'screen', 0, 1, 1, configurations, (end if smoke else started + seconds * 5 / 12)
    yield 'batch', 0, 1, (10 if smoke else 1000), configurations, (end if smoke else started + seconds * .75)
    if not smoke:
        for case in configurations:
            for rep in (1, 2, 3):
                for seed in (1, 2):
                    yield 'confirm', seed, rep, 1000, [case], end


def snapshot(directory, seconds, cap, smoke):
    import mlc_mdcrt_exact.bidirectional as backend
    import mlc_mdcrt_exact.bidirectional_inputs as inputs
    local = [Path(__file__), Path(__file__).with_name('mdcrt_adjugate_global.py'),
             Path(__file__).with_name('mdcrt_adjugate_dual.py')]
    sources = [(p, 'local') for p in local] + [(Path(p), 'backend') for p in
               (backend.__file__, inputs.__file__)]
    hashes = {}
    with zipfile.ZipFile(directory / 'sources.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
        for source, group in sources:
            name = f'{group}/{source.name}'
            archive.write(source, name)
            hashes[name] = hashlib.sha256(source.read_bytes()).hexdigest()
    manifest = {'seconds': seconds, 'trial_seconds': cap, 'smoke': smoke, 'bits': 16,
                'seed_base': BASE_SEED, 'grid': grid(smoke), 'methods': METHODS, 'queries': QS,
                'python': sys.version, 'flint': flint.__version__, 'platform': platform.platform(),
                'threads': 1, 'source_sha256': hashes, 'executable': sys.executable}
    (directory / 'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    (directory / 'PROTOCOL.md').write_text(
        '# 冻结实验协议\n\n'
        '四条路线使用同一批原始模矩阵和余数，全部返回特解和完整格基。'
        '两条伴随路线分别执行相同的逐块转换并计时。\n\n'
        '种子 0 先做单次检查，再做批量序列。种子 1、2 各三次。'
        '阶段分别截止于预算的 5/12、3/4 和终点，前一阶段提前结束时立即进入下一阶段。'
        '每个工作进程串行运行，硬限包含进程启动、计时、核验和存证。'
        '进程墙钟与算法计时分开。达到限时前完成且核验的检查点会保留。\n\n'
        '每组输入生成 32 个模矩阵并以相同的余数前缀比较 k=8、16、32。'
        '用有界整数见证生成余数。输入位长参数为 16，实际位长单独记录。\n\n'
        '基参考来自未压缩的原始逆矩阵行。逐一代回原同余，'
        '基与参考做双向整数坐标检查，失败立即暂停。'
        '原始数据通过本地 gzip/pickle 存证。报告再从事件重建。\n', encoding='utf-8')


def execute(directory, seconds, cap, smoke):
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=False)
    snapshot(directory, seconds, cap, smoke)
    started = monotonic()
    journal = directory / 'jobs.jsonl'
    generation = {}
    halted = False
    for phase, seed, rep, target, cases, cutoff in schedule(smoke, started, seconds):
        if halted or monotonic() >= started + seconds:
            break
        print(f'Phase {phase} seed={seed} repeat={rep}', flush=True)
        for family, n, k in cases:
            if monotonic() >= cutoff or halted:
                break
            key = (family, n, seed)
            location = directory / 'inputs' / f'{family}_n{n}_s{seed}'
            if key not in generation or generation[key] == 'budget_cutoff':
                arguments = ['--generate', str(location), family, str(n), str(seed)]
                status, wall = serial_worker(location, arguments, cap, cutoff)
                generation[key] = status
                answer = next((e for e in events(location / 'events.jsonl') if e['event']=='generated'), {})
                if status == 'ok' and not answer:
                    status = generation[key] = 'error'
                emit(journal, {'op':'generate', 'family':family, 'n':n, 'seed':seed,
                               'status':status, 'wall_s':wall, **answer})
                if status == 'error':
                    halted = True
            if generation[key] != 'ok':
                continue
            shift = (seed + rep + n.bit_length() + k.bit_length() + FAMILIES.index(family)) % len(METHODS)
            for method in METHODS[shift:] + METHODS[:shift]:
                if monotonic() >= cutoff:
                    break
                name = f'{phase}_{family}_n{n}_k{k}_s{seed}_r{rep}_{method}'
                job = directory / 'jobs' / name
                args = ['--trial', str(job), str(location / 'input.pkl.gz'), str(k), method, str(target)]
                status, wall = serial_worker(job, args, cap, cutoff)
                entries = events(job / 'events.jsonl')
                last = entries[-1]['event'] if entries else None
                if status == 'ok' and last != 'complete':
                    status = 'error'
                if status != 'ok' and last in ('verify', 'evidence'):
                    status = f'{last}_{status}'
                verified = [p['q'] for p in entries if p['event']=='checkpoint']
                emit(journal, {'op':'trial', 'phase':phase, 'family':family, 'n':n, 'k':k, 'seed':seed,
                               'rep':rep, 'method':method, 'target':target,'status':status,'wall_s':wall,
                               'path':str(job.relative_to(directory)), 'verified_q':verified})
                print(f'{family} n={n} k={k} {method}: {status} {verified}', flush=True)
                if status.endswith('error'):
                    halted = True
                    break
    (directory / 'completion.json').write_text(json.dumps(
        {'wall_seconds': monotonic() - started, 'halted': halted}, indent=2), encoding='utf-8')
    report(directory)
    if halted:
        raise SystemExit('Experiment halted on a worker error. See its worker.log.')


def write_csv(path, records, fields):
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def report(directory):
    manifest = json.loads((directory / 'manifest.json').read_text(encoding='utf-8'))
    completion = directory / 'completion.json'
    stop = json.loads(completion.read_text(encoding='utf-8')) if completion.exists() else {}
    trials = [job for job in events(directory / 'jobs.jsonl') if job['op']=='trial']
    raw = []
    for job in trials:
        entries = events(directory / job['path'] / 'events.jsonl')
        meta = next((e for e in entries if e['event']=='prepared'), {})
        raw.extend({**job, **meta, **point} for point in entries if point['event']=='checkpoint')
    dims = ('phase', 'family', 'n', 'k', 'seed', 'rep', 'method')
    write_csv(directory / 'measurements.csv', raw, [*dims, 'q', 'status', 'total_ns', 'prepare_ns',
              'query_ns', 'query_hnf_ns', 'validation_ns', 'conversion_ns', 'decomposition_setup_ns',
              'scalar_rows', 'rows_before', 'common_factor', 'coefficient_bits', 'modulus_bits',
              'unified_modulus_bits', 'input_bits', 'determinant_bits', 'work_shapes', 'basis_bits',
              'basis_hash', 'evidence_sha256', 'wall_s', 'path'])
    seen = {tuple(j[d] for d in dims):j for j in trials}
    coverage = []
    for phase, seed, rep, target, cases, _ in schedule(manifest['smoke'], 0, manifest['seconds']):
        for family, n, k in cases:
            for method in METHODS:
                key = (phase, family, n, k, seed, rep, method)
                coverage.append(seen.get(key, dict(zip(dims, key)) | {'status':'not_run', 'target':target}))
    write_csv(directory / 'coverage.csv', coverage, [*dims,'target','status','verified_q','wall_s','path'])
    cells = {}
    for p in raw:
        if p['phase']=='batch' and p['q']==1:
            continue  # Screening Q=1 is independent from the seed-0 batch.
        key = tuple(p[d] for d in ('family','n','k','seed','q','method'))
        cells.setdefault(key, []).append(p)
    summaries = {}
    for key, values in cells.items():
        required = 1 if key[3] == 0 else 3
        if len(values) == required:
            assert len({v['basis_hash'] for v in values}) == 1
            summaries[key] = {'time':median(v['total_ns'] for v in values),
                              'basis_hash':values[0]['basis_hash'],
                              'm':values[0].get('scalar_rows')}
    pairs = []
    for key, result in sorted(summaries.items()):
        family,n,k,seed,q,method = key
        for baseline in METHODS:
            if baseline == method:
                continue
            base = summaries.get((*key[:-1], baseline))
            if base:
                assert base['basis_hash'] == result['basis_hash']
                pairs.append(dict(family=family,n=n,k=k,seed=seed,q=q,method=method,baseline=baseline,
                                  trial_count=1 if seed==0 else 3,candidate_ns=result['time'],
                                  baseline_ns=base['time'],speedup=base['time']/result['time']))
    write_csv(directory / 'paired.csv', pairs, ['family','n','k','seed','q','method','baseline',
                                               'trial_count','candidate_ns','baseline_ns','speedup'])
    wall = f'{stop["wall_seconds"]:.1f} 秒' if stop else '运行中'
    lines = ['# 四条 MD-CRT 路线的配对实验', '',
             f'正式墙钟 {wall}。已核验计时检查点 {len(raw)} 个。'
             f'方法试跑状态 {dict(Counter(row["status"] for row in coverage))}。', '',
             '四条路线从相同的原始模矩阵和精确余数开始，输出一个特解与完整周期格。'
             '原生堆叠、原生级联与伴随辅助变量路线复用分解。同余行格复用伴随归约与周期格，'
             '每个新右端仍做一次普通 HNF。总耗时包含准备、转换和查询。'
             '生成、核验、存证及进程启动占运行预算，单列于算法耗时之外。', '',
             '数字为毫秒。种子 0 的单次和批量结果是筛选数据。'
             '种子 1、2 各自取三次中位数，仅展示完成三次的单元。'
             '四条路线都有数据时才标最快。配对加速比使用同一输入、种子和 Q。', '']
    for family, label, description in zip(FAMILIES, LABELS, DESCRIPTIONS):
        lines += [f'## {label}', '', description, '',
                  '| n | k | Q | 种子 | 约简行数 | 原生堆叠 | 原生级联 | 伴随＋辅助变量 | 伴随＋同余行格 | 最快 |',
                  '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
        for kind,n,k in manifest['grid']:
            if kind != family:
                continue
            for q in QS:
                for seed in (0,1,2):
                    cells_here = [summaries.get((family,n,k,seed,q,m)) for m in METHODS]
                    if seed and not any(cells_here):
                        continue
                    times = [f'{p["time"]/1e6:.3f}' if p else '—' for p in cells_here]
                    scalar_rows = next((p['m'] for p in cells_here[2:] if p), '—')
                    winner = min(range(4), key=lambda i:cells_here[i]['time']) if all(cells_here) else None
                    lines.append(f'| {n} | {k} | {q} | {seed} | {scalar_rows} | '+' | '.join(times)+
                                 f' | {NAMES[METHODS[winner]] if winner is not None else "—"} |')
        lines.append('')
    lines += ['原始计时与矩阵尺寸、位长见 measurements.csv，任务覆盖见 coverage.csv，'
              '各路线间的配对见 paired.csv。超时记录不参与加速比。']
    (directory / 'REPORT.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser()
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument('--self-test', action='store_true')
    actions.add_argument('--out', type=Path)
    actions.add_argument('--generate', nargs=4)
    actions.add_argument('--trial', nargs=5)
    actions.add_argument('--report', type=Path)
    parser.add_argument('--seconds', type=float, default=3600)
    parser.add_argument('--cap', type=float, default=30)
    parser.add_argument('--smoke', action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
    elif args.generate:
        out, family, n, seed = args.generate
        generate(Path(out), family, int(n), int(seed))
    elif args.trial:
        out, data, k, method, target = args.trial
        trial(Path(out), Path(data), int(k), method, int(target))
    elif args.report:
        report(args.report)
    else:
        execute(args.out, args.seconds, args.cap, args.smoke)


if __name__ == '__main__':
    main()
