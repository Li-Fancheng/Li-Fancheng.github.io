# 趋化：从随机游走到集体模式形成

西湖大学 2026 年秋季《定量生物学导论》课程学习笔记。

正文见 [趋化笔记](lecture_notes.md)，从一维随机游走推导扩散和漂移，再加入吸引物分泌与化学场反馈，建立 Keller–Segel 类型的最简趋化模型。图像为可复现的理论图和数值模拟结果。

## 文件布局

```text
chemotaxis/
├── lecture_notes.md                 # 笔记正文，图片为相对路径
├── README.md                        # 本说明
├── requirements.txt                 # NumPy 和 Matplotlib
├── figures/
│   ├── 01_random_walk_trajectories.png
│   ├── 02_random_walk_to_diffusion.png
│   ├── 03_biased_random_walk.png
│   ├── 04_chemotaxis_in_fixed_field.png
│   ├── 05_field_mediated_feedback.png
│   ├── 05_field_mediated_feedback.svg
│   ├── 06_collective_1d_space_time.png
│   ├── 07_collective_1d_snapshots.png
│   ├── 08_optional_growth_rate.png
│   └── 09_optional_2d_patterns.png
└── code/
    ├── common.py
    ├── random_walk.py
    ├── fixed_field.py
    ├── feedback_diagram.py
    ├── collective_1d.py
    ├── collective_2d.py
    └── make_all_figures.py
```

## 阅读笔记

打开 `lecture_notes.md`（VS Code 的 Markdown 预览、Typora、Obsidian 或支持数学公式的 Markdown 阅读器均可）。请保持 `figures/` 与正文的相对路径不变。

笔记使用标准 Markdown 加 `$...$`、`$$...$$` LaTeX 数学语法。图片是本包代码生成的原创示意图／模拟图，不需要访问网络。

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
python code/collective_1d.py
python code/collective_2d.py
```

运行后图像会输出到 `figures/`。计算依赖 NumPy 与 Matplotlib，不需要 SciPy，也不需要网络连接。脚本随机种子固定，输出结果可复现。

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

正文末尾给出了 Keller–Segel、Berg–Brown 和 Budrene–Berg 的经典论文 DOI。材料中的所有插图均为本项目源代码生成；参考文献提供现象背景，不代表本包复刻了论文中的实验结果。
