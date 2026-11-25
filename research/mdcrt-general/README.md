# 一般 MD-CRT 研究材料

本目录汇总一般矩阵中国剩余定理（MD-CRT）的伴随归约、全局 HNF 求解方法及已完成的性能实验。

## 内容

- theory/adjugate-hnf-guide.md 从问题定义开始，介绍伴随归约、逐块约简、两种全局 HNF 构造、正确性引理与数值算例。文档第 1.1 节也给出一般 MD-CRT 的有解充要条件。
- algorithms/mdcrt_adjugate_global.py 实现伴随归约后的辅助变量 HNF 与含右端同余行格 HNF 两条路线。
- algorithms/mdcrt_adjugate_dual.py 实现伴随同余和对偶格路线。algorithms/tests/ 保存对应验证用例。
- experiments/four-routes/ 保存四种路线的冻结协议、完整报告、汇总、原始计时表、配对表、覆盖表及 16 张 3D 图。网页报告位于 charts_3d/index.html。
- experiments/adjugate-dual/ 保存对偶格路线的补充实验结果、计时数据和运行脚本。
- [experiments/random-large/REPORT.md](experiments/random-large/REPORT.md) 比较随机可解输入的矩阵级联和两种伴随 HNF，覆盖 n、k 各取 16、32、64，查询次数为 1 和 100。输入生成和核验规则见[实验协议](experiments/random-large/PROTOCOL.md)。
- requirements.txt 记录绘图和精确整数运算依赖版本。

四路线实验比较原生堆叠 HNF、原生矩阵级联、伴随＋辅助变量和伴随＋同余行格。输入按循环共享因子、稠密共享因子、稠密无公共标量因子、对角互素分组。全部结果和未完成覆盖均保留在报告与 CSV 中。

## 实验数据范围

正式实验生成的输入快照约 32 MB、逐进程运行日志约 2.6 GB，均未放入本目录。计时记录、核验覆盖、确定性生成参数和冻结源码快照均已保留。正式报告说明了运行时限和未完成项，缺少计时的位置没有用估算值填充。

## 依赖与源码

algorithms/mlc_mdcrt_exact/ 包含正式四路线实验使用的精确算术适配器。代码使用 FLINT 的 Python 接口。当前环境版本见 experiments/four-routes/manifest.json。归档清理了本机专属的绝对可执行文件路径。运行补充单元测试：

~~~bash
python -m unittest discover -s algorithms/tests -p "test_*.py"
~~~

运行入口和完整参数见 experiments/four-routes/PROTOCOL.md 及相应脚本的 --help 输出。
