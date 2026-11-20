"""Render the verified seed-0 MD-CRT timings as grouped 3D bars."""

import argparse
from base64 import b64encode
import csv
from html import escape
from math import log10
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager, pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


FAMILIES = (
    ("cyclic_shared", "循环共享因子", "模矩阵共享因子 2，约掉后每块的整数商群为循环群。"),
    ("dense_shared", "稠密共享因子", "非奇异稠密模矩阵由有符号偶数元素生成。"),
    ("dense_primitive", "稠密无公共标量因子", "稠密模矩阵的全部元素全局 gcd 为 1。"),
    ("diagonal_coprime", "对角互素", "模矩阵为互不相同的 16 位素数乘单位矩阵。"),
)
METHODS = (
    ("global_tri", "原生堆叠 HNF", "#3366A3"),
    ("cascade_tri", "原生矩阵级联", "#C36519"),
    ("adj_aux", "伴随＋辅助变量", "#087F71"),
    ("adj_congruence", "伴随＋同余行格", "#895294"),
)
SIZES = (8, 16, 32)
QUERIES = (1, 10, 100, 1000)
FLOOR_MS = 0.5
CEILING_MS = 40000
TICKS_MS = (0.5, 1, 10, 100, 1000, 10000)
MISSING = {
    "timeout": "试跑限时",
    "verify_timeout": "核验限时",
    "evidence_timeout": "存证限时",
    "budget_cutoff": "预算截断",
    "verify_budget_cutoff": "核验时预算截断",
    "not_run": "未运行",
}


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def experiment_data(directory):
    values = {}
    for point in read_csv(directory / "measurements.csv"):
        q = int(point["q"])
        phase = "screen" if q == 1 else "batch"
        if point["seed"] != "0" or point["phase"] != phase:
            continue
        key = (point["family"], int(point["n"]), int(point["k"]), q, point["method"])
        assert key not in values and point["evidence_sha256"], key
        values[key] = int(point["total_ns"]) / 1e6
    statuses = {
        (row["family"], int(row["n"]), int(row["k"]), row["phase"], row["method"]): row["status"]
        for row in read_csv(directory / "coverage.csv") if row["seed"] == "0"
    }
    expected = {(family, n, k, q, method) for family, *_ in FAMILIES for n in SIZES
                for k in SIZES for q in QUERIES for method, *_ in METHODS}
    assert values.keys() <= expected
    assert all(FLOOR_MS < value < CEILING_MS for value in values.values())
    return values, statuses


def height(milliseconds):
    return log10(milliseconds / FLOOR_MS)


def plot(path, family, label, q, values):
    fig = plt.figure(figsize=(9, 6.1), facecolor="#F7F9FC")
    ax = fig.add_subplot(projection="3d", computed_zorder=True, facecolor="#F7F9FC")
    ax.set_position([0.08, 0.145, 0.84, 0.73])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.90, 0.93, 0.97, 0.55))
        axis.line.set_color("#A3B0C1")
    legend = []
    valid = []
    missing = []
    for method_index, (method, name, color) in enumerate(METHODS):
        legend.append(Patch(facecolor=color, label=name))
        for x, n in enumerate(SIZES):
            for y, k in enumerate(SIZES):
                xpos = x + (method_index - 1.5) * 0.19
                value = values.get((family, n, k, q, method))
                if value is None:
                    missing.append((xpos, y, color))
                else:
                    valid.append((xpos, y, height(value), color))
    xs, ys, zs, colors = zip(*valid)
    ax.bar3d([x - 0.072 for x in xs], [y - 0.095 for y in ys], [0] * len(xs),
             0.144, 0.19, zs, color=colors, shade=True, edgecolor="#FFFFFF", linewidth=0.35)
    if missing:
        for x, y, color in missing:
            ax.plot([x], [y], [0.02], marker="x", markersize=7, markeredgewidth=1.8,
                    color=color, linestyle="none", zorder=30)
        legend.append(Line2D([], [], color="#526174", marker="x", linestyle="none",
                             label="无已核验计时"))
    ax.set(xlim=(-0.48, 2.48), ylim=(-0.42, 2.42), zlim=(0, height(CEILING_MS)),
           xticks=range(3), xticklabels=SIZES, yticks=range(3), yticklabels=SIZES)
    ax.yaxis.set_inverted(True)
    ax.set_zticks([height(t) for t in TICKS_MS])
    ax.set_zticklabels([f"{t:g}" for t in TICKS_MS])
    ax.set_xlabel("未知量维度 n", labelpad=10)
    ax.set_ylabel("模矩阵数 k", labelpad=12)
    ax.set_zlabel("总耗时 / ms（对数）", labelpad=10)
    ax.set_box_aspect((1.35, 1.08, 0.95), zoom=0.97)
    ax.view_init(elev=27, azim=-48)
    ax.tick_params(labelsize=9, colors="#39465A")
    fig.text(0.075, 0.943, f"{label} · Q = {q:,}", fontsize=20, weight="bold", color="#15283F")
    fig.text(0.076, 0.89, f"种子 0 · 一次准备 + {q:,} 次查询 · 已核验 {len(valid)}/36",
             fontsize=10, color="#46566D")
    fig.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 0.11), ncol=len(legend),
               frameon=False, fontsize=9.2, handlelength=1.2, columnspacing=1.3)
    fig.text(0.5, 0.055, "共同时间刻度 · 柱体从 0.5 ms 起绘制 · 顶端表示实测总耗时",
             ha="center", fontsize=9.6, color="#46566D")
    fig.savefig(path, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)
    return len(valid)


def timing_table(family, q, values, statuses):
    headers = ("n", "k", *(label for _, label, _ in METHODS))
    output = ["<table><thead><tr>" + "".join(f"<th>{h}</th>" for h in headers) + "</tr></thead><tbody>"]
    phase = "screen" if q == 1 else "batch"
    for n in SIZES:
        for k in SIZES:
            readings = [values.get((family, n, k, q, method)) for method, *_ in METHODS]
            fastest = min((v for v in readings if v is not None), default=None) if all(readings) else None
            fields = [f"<th>{n}</th>", f"<th>{k}</th>"]
            for value, (method, _, _) in zip(readings, METHODS):
                if value is None:
                    status = statuses.get((family, n, k, phase, method), "not_run")
                    fields.append(f'<td class="missing">{escape(MISSING.get(status, "缺少有效计时"))}</td>')
                elif value == fastest:
                    fields.append(f'<td class="fast">{value:,.3f}</td>')
                else:
                    fields.append(f"<td>{value:,.3f}</td>")
            output.append("<tr>" + "".join(fields) + "</tr>")
    return "".join(output) + "</tbody></table>"


def render(source, output):
    font = Path("C:/Windows/Fonts/msyh.ttc")
    if font.is_file():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=font).get_name()
    plt.rcParams["axes.unicode_minus"] = False
    values, statuses = experiment_data(source)
    output.mkdir(parents=True, exist_ok=True)
    navigation = " ".join(f'<a href="#{family}">{label}</a>' for family, label, _ in FAMILIES)
    sections = []
    for family, label, description in FAMILIES:
        cards = []
        for q in QUERIES:
            path = output / f"{family}_q{q}.png"
            verified = plot(path, family, label, q, values)
            data = b64encode(path.read_bytes()).decode("ascii")
            table = timing_table(family, q, values, statuses)
            cards.append(f"""<article>
<img src="data:image/png;base64,{data}" alt="{label}，Q={q}，四种方法的分组 3D 耗时柱状图，{verified}/36 个有效计时">
<details><summary>查看 Q = {q:,} 的精确时间（ms）和限时状态</summary>{table}</details>
</article>""")
        sections.append(f'<section id="{family}"><h2>{label}</h2><p>{description}</p>'
                        f'<div class="charts">{"".join(cards)}</div></section>')
    html = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MD-CRT 四路线 3D 性能对比</title>
<style>
:root {{color-scheme: light; font-family: "Microsoft YaHei", "Noto Sans SC", sans-serif; color:#16283f; background:#f0f3f8}}
body {{margin:0; padding:clamp(16px,4vw,64px); line-height:1.65}}
main {{max-width:1500px; margin:auto}}
h1 {{font-size:clamp(28px,4vw,46px); margin:0 0 .3em}}
h2 {{font-size:26px; margin:0}}
p {{color:#42536a; margin:.45em 0 1em}}
nav {{display:flex; gap:.55rem; flex-wrap:wrap; margin:1.15rem 0 1.8rem}}
nav a {{background:#fff; padding:.38rem .8rem; border:1px solid #d5deea; border-radius:6px; color:#17477a}}
section {{margin-bottom:2.85rem}}
.charts {{display:grid; grid-template-columns:repeat(auto-fit,minmax(min(100%,530px),1fr)); gap:16px}}
article {{background:white; border:1px solid #dbe3ed; border-radius:8px; overflow:hidden}}
img {{display:block; width:100%; height:auto}}
details {{padding:.65rem .8rem; overflow:auto; font-size:.87rem}}
summary {{cursor:pointer; color:#17477a; font-weight:600}}
table {{width:100%; border-collapse:collapse; margin-top:.6rem; font-variant-numeric:tabular-nums; white-space:nowrap}}
th,td {{text-align:right; padding:.26rem .4rem; border-bottom:1px solid #e4eaf2}}
thead th {{background:#eff3fa}} .missing {{color:#526174}} .fast {{font-weight:750; color:#006650}}
</style></head><body><main>
<header><h1>MD-CRT 四路线 3D 性能对比</h1>
<p>横轴 n 为未知量维度，纵深轴 k 为模矩阵数，柱高表示包含一次预处理的总耗时，单位毫秒。四种方法的颜色在所有图中保持一致。</p>
<p>统一使用覆盖九种规模的种子 0。Q = 1 取单次筛选，其余 Q 取批量检查点。独立种子用于原实验的确认统计。各图时间轴共用同一对数刻度，柱底为 0.5 ms，彩色 × 表示对应方法缺少该点的已核验计时。精确表中只给四条路线都有数据的行标出最快值。</p></header>
<nav aria-label="按输入结构跳转">{navigation}</nav>
{"".join(sections)}
</main></body></html>"""
    (output / "index.html").write_text(html, encoding="utf-8")
    print(f"Rendered {len(FAMILIES) * len(QUERIES)} charts with {len(values)} verified points.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    arguments = parser.parse_args()
    render(arguments.input, arguments.output)
