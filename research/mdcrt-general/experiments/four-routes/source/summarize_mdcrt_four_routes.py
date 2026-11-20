"""Summarize only complete, matched comparisons from the frozen experiment."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
from statistics import geometric_mean, median

METHODS = ('global_tri', 'cascade_tri', 'adj_aux', 'adj_congruence')
NAMES = ('原生堆叠', '原生级联', '伴随＋辅助变量', '伴随＋同余行格')
LABELS = dict(zip(METHODS, NAMES))


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def summary(directory):
    rows = read_csv(directory / 'paired.csv')
    coverage = read_csv(directory / 'coverage.csv')
    measurements = read_csv(directory / 'measurements.csv')
    completion = json.loads((directory / 'completion.json').read_text(encoding='utf-8'))
    cells = defaultdict(dict)
    for row in rows:
        if row['baseline'] != 'global_tri':
            continue
        key = (row['family'], *(int(row[field]) for field in ('n', 'k', 'seed', 'q')))
        cells[key][row['method']] = int(row['candidate_ns'])
        cells[key]['global_tri'] = int(row['baseline_ns'])
    cohorts = defaultdict(list)
    for key, methods in cells.items():
        if all(method in methods for method in METHODS):
            cohorts[key[3] != 0, key[4]].append((key, methods))

    lines = ['# 四路线 MD-CRT 实验结论', '',
             f'实测墙钟 {completion["wall_seconds"]:.2f} 秒，'
             f'{len(measurements)} 个检查点通过原式代回和完整格基核验。', '',
             '同一行只统计四条路线均有核验结果的输入。筛选表每个配对用种子 0 的一次试跑，'
             '确认表用独立种子 1、2 各自完成三次的中位数。', '',
             '| 阶段 | Q | 四路线齐全的配对 | 堆叠最快 | 级联最快 | 辅助变量最快 | 同余行格最快 |',
             '|---|---:|---:|---:|---:|---:|---:|']
    for (confirmed, q), cases in sorted(cohorts.items()):
        wins = Counter(min(METHODS, key=times.__getitem__) for _, times in cases)
        stage = '确认' if confirmed else '筛选'
        lines.append(f'| {stage} | {q} | {len(cases)} | ' +
                     ' | '.join(str(wins[method]) for method in METHODS) + ' |')

    lines += ['', '以下是三组方法之间的配对几何平均加速比。'
              '每行采用同一批四路线齐全的配对，系数大于 1 表示右侧路线更快。', '',
              '| 阶段 | Q | 堆叠 / 辅助变量 | 级联 / 辅助变量 | 堆叠 / 同余行格 | 级联 / 同余行格 | 辅助变量 / 同余行格 |',
              '|---|---:|---:|---:|---:|---:|---:|']
    comparisons = [(0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]
    for (confirmed, q), cases in sorted(cohorts.items()):
        factors = [geometric_mean(times[METHODS[a]] / times[METHODS[b]]
                                  for _, times in cases) for a, b in comparisons]
        stage = '确认' if confirmed else '筛选'
        lines.append(f'| {stage} | {q} | ' + ' | '.join(f'{factor:.3f}' for factor in factors) + ' |')

    lines += ['', '## 各结构的单次结果', '',
              '| 输入结构 | 原始约束数 | 约简行数中位数 | 堆叠 / 辅助变量 | 堆叠 / 同余行格 |',
              '|---|---:|---:|---:|---:|']
    groups = defaultdict(list)
    for key, times in cohorts[False, 1]:
        groups[key[0]].append((key, times))
    after = {(row['family'], int(row['n']), int(row['k'])):int(row['scalar_rows'])
             for row in measurements if row['phase']=='screen' and
             row['method']=='adj_aux' and row['scalar_rows']}
    for family, cases in sorted(groups.items()):
        before = [n*k for (_, n, k, _, _), _ in cases]
        reduced = [after[f,n,k] for (f,n,k,_,_), _ in cases]
        factors = [geometric_mean(times['global_tri']/times[m] for _,times in cases)
                   for m in METHODS[2:]]
        lines.append(f'| {family} | {min(before)}–{max(before)} | {median(reduced):g} | '
                     f'{factors[0]:.3f} | {factors[1]:.3f} |')

    lines += ['', '## 实际覆盖和解释', '']
    for phase in ('screen','batch','confirm'):
        counts = Counter(row['status'] for row in coverage if row['phase'] == phase)
        lines.append(f'- {phase} 任务状态：{dict(counts)}。')
    lines += [
        '- 辅助变量构造中带变换 HNF 的宽度由约简后保留的标量行数决定。'
        '同余行格构造的宽度固定为 n+1，仍需读入约简行并处理统一模数的位长。',
        '- 循环共享因子输入可把块约简到一条标量同余。对角互素块通常保留 n 条。'
        '稠密族各实例约简行数及统一模数位长见 measurements.csv。',
        '- 同余行格路线每个右端重做 HNF。大 Q 中的实际 HNF 时间可用 query_hnf_ns '
        '与 query_ns 对照，三条缓存路线的分解成本只在 prepare_ns 中发生一次。',
        '- 各 Q 的配对样本可能不同。限时点没有被替换成固定耗时，不能从不同覆盖的横向平均推断交叉点。',
        '- 完整的九种规模分表及每个种子的比较见 REPORT.md。'
    ]
    (directory / 'SUMMARY.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    print('\n'.join(lines[:25]))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', type=Path)
    summary(parser.parse_args().directory)
