# 四路线正式实验

本目录对应一小时正式实验。PROTOCOL.md 记录冻结协议，REPORT.md 是逐结构结果，SUMMARY.md 汇总配对速度比较，measurements.csv、paired.csv 和 coverage.csv 保存计时与覆盖数据。manifest.json 记录环境和随机种子。

## 图表

在本地浏览器打开 charts_3d/index.html 可查看四类输入、每类四个查询数 Q 的分组 3D 柱状图。每张图以 n 和 k 为平面坐标，以总毫秒数为高度，四条路线用固定颜色区分。图表使用共同对数时间刻度，缺少核验计时的位置以叉号标出。网页内含各点精确计时表，PNG 可单独查看。

## 保留与省略

CSV 中保留每个已计时检查点和覆盖状态。原始 84 个压缩输入样本约 32 MB，逐进程日志约 2.6 MB，未加入此归档。source/sources_frozen.zip 保存实验源码冻结快照，便于校对 manifest.json 中列出的源码哈希。
