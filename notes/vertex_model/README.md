# Vertex Model：从细胞几何到组织力学

八篇主线和一篇拓展，入口为 `/notes/vertex_model/`。页面元信息在 `content/zh/notes/morphogenesis/vertex_model/`，保留显式 URL；正文由 `shared-note` 读取本目录。

| 部分 | 源文件 |
| --- | --- |
| I · 几何与刚性 | [01 几何、能量与力学平衡](foundations/body.md)；[02 线性响应与非仿射弹性](linear_response/body.md)；[03 几何不相容、预应力与刚性转变](rigidity/body.md) |
| II · 重排与动力学 | [04 T1 重排、能量势垒与塑性](rearrangements/body.md)；[05 组织流变与主动动力学](rheology_activity/body.md)；[06 涨落、相关与组织重塑](fluctuations/body.md) |
| III · 几何与生物推广 | [07 各向异性、曲面与三维组织](geometry_3d/body.md) |
| IV · 实验与机制判别 | [08 形态发生、实验检验与开放问题](morphogenesis/body.md) |
| 专题拓展 | [主动软物质的堆积与拓扑](active_packing/body.md) |

果蝇和斑马鱼是主要发育案例；培养上皮与线虫提供流变、运动停滞、边界、细胞挤出和骨架更新的比较。经典汇合 Vertex、非汇合界面、细胞中心接触图及连续模型分别标明自由度与假设。首页集中给出理论目录，生物案例和视频索引放在可展开区域，不重复展示文章卡片。

## 贯穿计算：同一周期网络

[computation/README.md](computation/README.md) 给出定义、运行方法与边界。16 细胞、32 顶点的周期网络使用同一组初态与拓扑，连接以下量：

- Note 03：面积与周长残差、完整 Hessian、Gram 项、预应力项、最低模态、仿射与非仿射剪切模量。
- Note 04：同一交界的受约束松弛收缩与展开路径，零边长处切换邻接。
- Note 05：明确的 OU 张力、耗散、T1 规则及对称有限阶跃应力差分，保留逐条轨迹与误差范围。

静态模量另用重新最小化核对；Hessian 另用力的方向差分核对；应力另用能量虚功核对。源代码、JSON 数据、矩阵、校验和复现包保存在 computation 目录。没有将受约束路径宣称为全局最小能量路径，也没有将有限时间、有限幅度响应宣称为无限时间线性模量。

Note 05 集中推导运动学分解、能量受力关系和重排闭合如何组成 Maxwell 方程，并比较分裂、移除与参考结构更新。单模态 OU 方差和频谱推导集中在 Note 06。长期平台统一记为 `G_eq`，并区分固定拓扑平台与瞬时模量。

## 其他计算图的复现

```bash
python3 -m venv /tmp/vertex-model-figures
/tmp/vertex-model-figures/bin/pip install -r notes/vertex_model/requirements.txt
/tmp/vertex-model-figures/bin/python notes/vertex_model/generate_figures.py
/tmp/vertex-model-figures/bin/python notes/vertex_model/generate_active_matter_figures.py
/tmp/vertex-model-figures/bin/python notes/vertex_model/generate_fluctuation_figures.py
```

脚本导出 SVG 和相应 300 dpi PNG，数学字体为 Computer Modern，文字优先 Arial。前两个脚本不使用随机采样；涨落脚本的 OU 样本使用种子 20261020，其余曲线由公式和积分计算。贯穿网络的计算另外依赖 SciPy 与 Autograd，依赖列表和命令见 computation。

## 论文图与视频

原图保持下载时的字节，DOI、图号、来源、日期和 SHA-256 见 [paper_figure_sources.json](paper_figure_sources.json)。作者稿版本在正文标明，权利与许可按原文保留。

```bash
python3 notes/vertex_model/download_paper_figures.py
```

脚本先核对本地校验值，仅获取缺失文件；不覆盖与记录不同的文件。每幅图的网页宽度与类型见 [data/vertex_figures.json](../../data/vertex_figures.json)，点击可看原图。

视频使用期刊的远程原媒体和 HTML5 控件，不自动播放、不在仓库中存储或转换。实验、分析叠加与模拟分别标明；图注保留原文视频号、标记和时间尺度。来源、尺寸与时长见 [paper_video_sources.json](paper_video_sources.json)。
