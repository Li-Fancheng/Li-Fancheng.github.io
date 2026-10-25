# 同一周期 Vertex 网络

本算例连接 Note 03 的刚性、Note 04 的有限收缩路径和 Note 05 的受驱动响应。它是自洽的小网络计算，不用于拟合论文数据或估计普适临界点。

## 复现

```bash
python3 -m venv /tmp/vertex-unified
/tmp/vertex-unified/bin/pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /tmp/vertex-unified/bin/python unified_vertex.py
/tmp/vertex-unified/bin/python verify_unified_case.py
/tmp/vertex-unified/bin/python plot_unified_case.py
```

在仓库内运行时，给脚本加上 `notes/vertex_model/computation/` 前缀。解压包中直接运行计算和校验；绘图脚本在原仓库写入对应章节；单独解压时写入包内的 generated 目录。

默认主动计算为 24 对轨迹，`--repeats` 可更改数量。随机输入由种子 1000–1023 决定。实际使用的环境为 Python 3.13、NumPy 2.5.4、SciPy 1.18.1、Autograd 1.9.1 和 Matplotlib 3.11.2；不同优化库版本可能改变近零能量构型中的具体坐标。

## 定义与输出

- 盒边长 4，16 个细胞、32 个三价顶点，总面积 16。
- 初态为 `4×4` 方格中心加均匀位置扰动后生成的周期 Voronoi 多边形，种子 3。后续变量为交界顶点坐标。
- 能量为 `0.5 Σ[(A−1)²+(P−p0)²]`。每个目标周长独立使用同一初态和同一拓扑最小化。
- `J`、`H` 和残差加权的几何 Hessian `S` 分别求导，再检查分解；模量另外用重新最小化的能量中心差分核对。
- `unified_vertex_data.json` 包含制备、静态统计、最低模态、路径坐标和每对动力学轨迹。
- `unified_vertex_matrices.npz` 保存 `p0=3.75` 的顶点、`J,H,S`、含仿射剪切的扩展 Hessian 与平移投影矩阵。
- `verification.json` 保留求解误差和独立数值审计；没有用预设曲线替换未收敛结果。

## 计算边界

T1 追踪顶点对 `(2,3)`。给定旧边和新边的方向、逐步改变边向量，再松弛其余坐标；在零边长处切换邻接。因此这是受约束局部路径的势垒，不是全空间最小能量路径搜索。不同反应坐标、其他先发生的重排或制备可以给出不同势垒。

主动计算使用细胞对标签上的独立 OU 张力场；非接触标签也继续演化。新接触读取其标签的场，不继承旧边张力。边长阈值为 0.025，新边长为 0.0375；换接不要求总能量单调下降。有限阈值引起的能量跳变已记录，尚未建立阈值趋零的收敛。

驱动与 `±0.15` 阶跃同时开启，记录对称有限幅度应力差分。它不自动等于稳态、线性极限中的 `G(t)`，24 对轨迹也不足以确定无限时间平台。较小的被动阶跃是振幅检查；持续剪切流动、有限尺寸和多制备平均需另行计算。
