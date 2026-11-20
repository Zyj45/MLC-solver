# 四类 MD-CRT 输入的 3D 性能图

横轴为未知量维度 n，纵深轴为模矩阵数 k，柱高为总耗时。颜色依次表示原生堆叠 HNF、原生矩阵级联、伴随＋辅助变量、伴随＋同余行格。时间轴使用共同对数刻度。

## 循环共享因子

### Q = 1

![循环共享因子输入，Q=1](./cyclic_shared_q1.png)

### Q = 10

![循环共享因子输入，Q=10](./cyclic_shared_q10.png)

### Q = 100

![循环共享因子输入，Q=100](./cyclic_shared_q100.png)

### Q = 1000

![循环共享因子输入，Q=1000](./cyclic_shared_q1000.png)

## 稠密共享因子

### Q = 1

![稠密共享因子输入，Q=1](./dense_shared_q1.png)

### Q = 10

![稠密共享因子输入，Q=10](./dense_shared_q10.png)

### Q = 100

![稠密共享因子输入，Q=100](./dense_shared_q100.png)

### Q = 1000

![稠密共享因子输入，Q=1000](./dense_shared_q1000.png)

## 稠密无公共标量因子

### Q = 1

![稠密无公共标量因子输入，Q=1](./dense_primitive_q1.png)

### Q = 10

![稠密无公共标量因子输入，Q=10](./dense_primitive_q10.png)

### Q = 100

![稠密无公共标量因子输入，Q=100](./dense_primitive_q100.png)

### Q = 1000

![稠密无公共标量因子输入，Q=1000](./dense_primitive_q1000.png)

## 对角互素

### Q = 1

![对角互素输入，Q=1](./diagonal_coprime_q1.png)

### Q = 10

![对角互素输入，Q=10](./diagonal_coprime_q10.png)

### Q = 100

![对角互素输入，Q=100](./diagonal_coprime_q100.png)

### Q = 1000

![对角互素输入，Q=1000](./diagonal_coprime_q1000.png)
