# 趋化：从随机游走到集体模式形成

西湖大学 2026 年秋季《定量生物学导论》课程学习笔记。

正文见[趋化笔记](lecture_notes.md)。从单细胞的随机游动开始，先计算细胞团怎样散开，再加入化学梯度和细胞分泌，建立 Keller–Segel 模型，判断微小密度起伏何时增长。最后回到实验，比较均匀背景信号对聚集的影响。技术附录 A 连接连续时间 Master equation、Fokker–Planck 与 Langevin 描述，附录 B 给出完整两场增长率，附录 C 推导化学响应核。

文中 18 张图片按阅读顺序统一编号为 Fig. 1–18（均在正文中）：4 张论文实验图，以及 14 张由代码生成的模型示意图、理论曲线和模拟图。图片文件名保留生成脚本中的标识，文中图号可按下表查找：

| 文中图号 | 图片文件（位于 `figures/`） |
|---|---|
| Fig. 1 | `00_ecoli_run_tumble_taute2015.png` |
| Fig. 2 | `01_random_walk_trajectories.png` |
| Fig. 3 | `A_walk_probability_bookkeeping.png` / `.svg` |
| Fig. 4 | `B_diffusion_boundary_flux.png` / `.svg` |
| Fig. 5 | `02_random_walk_to_diffusion.png` |
| Fig. 6 | `paper_grognot2021_fig1.png` |
| Fig. 7 | `03_biased_random_walk.png` |
| Fig. 8 | `C_gradient_to_jump_bias.png` / `.svg` |
| Fig. 9 | `04_chemotaxis_in_fixed_field.png` |
| Fig. 10 | `paper_laganenka2016_fig1.jpg` |
| Fig. 11 | `paper_laganenka2016_fig2.jpg` |
| Fig. 12 | `05_field_mediated_feedback.png` / `.svg` |
| Fig. 13 | `E_signal_range_and_superposition.png` / `.svg` |
| Fig. 14 | `D_cosine_mode_and_growth.png` / `.svg` |
| Fig. 15 | `08_optional_growth_rate.png` |
| Fig. 16 | `06_collective_1d_space_time.png` |
| Fig. 17 | `07_collective_1d_snapshots.png` |
| Fig. 18 | `09_optional_2d_patterns.png` |


## 文件布局

```text
chemotaxis/
├── lecture_notes.md                 # 笔记正文，图片为相对路径
├── README.md                        # 本说明
├── requirements.txt                 # NumPy 和 Matplotlib
├── figures/
│   ├── paper_grognot2021_fig1.png
│   ├── paper_laganenka2016_fig1.jpg
│   ├── paper_laganenka2016_fig2.jpg
│   ├── 00_ecoli_run_tumble_taute2015.png
│   ├── 01_random_walk_trajectories.png
│   ├── 02_random_walk_to_diffusion.png
│   ├── 03_biased_random_walk.png
│   ├── 04_chemotaxis_in_fixed_field.png
│   ├── 05_field_mediated_feedback.png
│   ├── 05_field_mediated_feedback.svg
│   ├── 06_collective_1d_space_time.png
│   ├── 07_collective_1d_snapshots.png
│   ├── 08_optional_growth_rate.png
│   ├── 09_optional_2d_patterns.png
│   ├── A_walk_probability_bookkeeping.png / .svg
│   ├── B_diffusion_boundary_flux.png / .svg
│   ├── C_gradient_to_jump_bias.png / .svg
│   ├── D_cosine_mode_and_growth.png / .svg
│   └── E_signal_range_and_superposition.png / .svg
└── code/
    ├── common.py
    ├── random_walk.py
    ├── fixed_field.py
    ├── feedback_diagram.py
    ├── teaching_diagrams.py
    ├── collective_1d.py
    ├── collective_2d.py
    └── make_all_figures.py
```

## 阅读笔记

打开 `lecture_notes.md`（VS Code 的 Markdown 预览、Typora、Obsidian 或支持数学公式的 Markdown 阅读器均可）。请保持 `figures/` 与正文的相对路径不变。

笔记使用标准 Markdown 加 `$...$`、`$$...$$` LaTeX 数学语法。图片保存在本地；正文中的实验视频链接需要访问网络。

## 安装及运行

进入 `notes/qbio_course/chemotaxis` 目录：

```bash
python -m pip install -r requirements.txt
python code/make_all_figures.py
```

若只想生成某一类图，可分别运行：

```bash
python code/random_walk.py
python code/fixed_field.py
python code/feedback_diagram.py
python code/teaching_diagrams.py
python code/collective_1d.py
python code/collective_2d.py
```

运行后图像会输出到 `figures/`。示意图使用英文标注，正文配中文图注。SVG 中的文字已转为路径，查看时不需要安装原字体。计算依赖 NumPy 与 Matplotlib，不需要 SciPy，也不需要网络连接。脚本随机种子固定，输出结果可复现。

## 参数比较

- `code/random_walk.py`：将 `p_right=0.58` 改为 `0.50`、`0.55`、`0.65`，比较平均位移与分布宽度。
- `code/fixed_field.py`：改变 `chi` 或 `D`，比较外加浓度峰下的积累程度。
- `code/collective_1d.py`：改 `low_chi`、`high_chi`（例如 0.25、0.45、0.7），研究化学反馈是否使扰动增长。物理参数在 `PARAMS` 中。
- `code/collective_2d.py`：二维图展示短时演化，不建议把很晚时刻的尖锐峰解释为有固定大小的稳定菌落。

## 数值方法与正确性检查

一维、二维细胞密度使用守恒的有限体积法（finite-volume method），趋化漂移采用迎风通量（upwind flux），化学物质采用显式反应–扩散更新。程序自动限制时间步长并检查非负性及细胞总数守恒。

采用**周期边界条件**（用于避免边界输入/输出干扰）。不同图片的色条可能使用不同范围；不能只通过颜色判断两组结果的绝对大小。

模型将细菌的直游与转向（run-and-tumble）压缩为有效运动偏置，未显式加入感知记忆、拥挤排斥或趋化饱和。尤其是最简二维 Keller–Segel 在后期可能出现分辨率敏感的高密度峰，所以二维图故意只截取短时间窗。

## 内容与引用

正文末尾给出了 Keller–Segel、Berg–Brown 和 Budrene–Berg 的经典论文 DOI。Fig. 1 摘自 Taute 等（2015），原文 Fig. 4a，DOI: https://doi.org/10.1038/ncomms9776，按 CC BY 4.0（https://creativecommons.org/licenses/by/4.0/）转载，仅裁出面板，未改动数据。来源 PDF：https://ir.amolf.nl/pub/6550/16036OA.pdf，第 5 页。另有 Fig. 6 引自 Grognot 与 Taute（2021），原文 Fig. 1，DOI: https://doi.org/10.1038/s42003-021-02190-2；Fig. 10、11 引自 Laganenka 等（2016），原文 Fig. 1、2，DOI: https://doi.org/10.1038/ncomms12984。这三张图片取自出版社提供的完整图片，未经修改，均按 CC BY 4.0 转载。其余 14 张图由本目录代码生成；生成脚本不会覆盖论文图。
