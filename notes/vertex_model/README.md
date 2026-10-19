# Vertex Model：从细胞几何到组织力学

五部分、七篇讲义式理论笔记；网站入口为 `/notes/vertex_model/`。正文包含逐步推导、完整算例、论文读图、适用条件与技术附录，共 34 张论文原图和 25 张自制图。

| 部分 | 源文件 |
| --- | --- |
| I · 模型建立 | [01 几何、能量与力学平衡](foundations/body.md) |
| II · 刚性的几何起源 | [02 线性响应与非仿射弹性](linear_response/body.md)；[03 几何不相容、预应力与刚性转变](rigidity/body.md) |
| III · 从固体到流动 | [04 T1 重排、能量势垒与塑性](rearrangements/body.md)；[05 组织流变与主动动力学](rheology_activity/body.md) |
| IV · 几何推广 | [06 各向异性、曲面与三维组织](geometry_3d/body.md) |
| V · 生物应用与边界 | [07 形态发生、实验检验与开放问题](morphogenesis/body.md) |

页面元信息在 `content/zh/notes/vertex_model/`，正文通过现有 `shared-note` shortcode 渲染，图与正文保存在同一主题目录内。

## 斑马鱼与物理主线

保留七篇结构，在 01、03、05、06、07 中贯穿斑马鱼：自由表面与接触角、胚盘刚性渗流、主动液滴流变、MPZ–PSM 材料梯度、黏附与密度的独立干预，以及细胞外通道对 Nodal 运输的反馈。经典汇合 Vertex Model、非汇合界面模型、有效接触图与连续运输分别标明自由度与适用范围。

## 复现图示

需要 Python 3、NumPy 和 Matplotlib：

```bash
python3 -m venv /tmp/vertex-model-figures
/tmp/vertex-model-figures/bin/pip install -r notes/vertex_model/requirements.txt
/tmp/vertex-model-figures/bin/python notes/vertex_model/generate_figures.py
```

脚本导出 25 张 SVG 和相应 300 dpi PNG，数学字体为 Computer Modern，文字优先 Arial；无随机采样。包括正多边形几何、真实仿射六边形路径、线性模态、预应力弹簧、激活率分布、Maxwell 响应、自由主动粒子与参考长度可塑性，以及接触角力平衡、铰链约束矩阵、重排活性与方向偏置、孔隙中的溶质守恒等计算。正文明确标注模型性质与参数，不把教学模型当作完整组织仿真或实验数据。

## 论文原图与展示规范

原图文件保持下载时的原始字节，来源、DOI、原文图号、检索日期与 SHA-256 见 [paper_figure_sources.json](paper_figure_sources.json)。图注保留作者与来源；其中两张 SVG 使用正式论文对应的作者稿图像，并在正文明确标注版本。

需要检查或重新获取原图时运行：

```bash
python3 notes/vertex_model/download_paper_figures.py
```

脚本优先检查本地校验值，仅下载缺失文件；远端文件与记录不一致时不覆盖已有图。图片权利与许可按原文保留。

所有 59 张图的网页宽度与类型记录在仓库 [data/vertex_figures.json](../../data/vertex_figures.json)。宽幅多面板论文图与单一几何示意分别设置展示宽度，点击均可查看原图；使用与 qbio_course 相同的长文、公式滚动、目录和图片展示方式。

## 论文视频

六篇正文中嵌入 8 段原始视频：线虫皮层切割、果蝇顶点滑移、斑马鱼刚性簇图谱、A-PSM 细胞起伏、孔隙与信号响应、胚盘铺展和不同材料条件下的接触动态，以及一段主动泡沫模型动画。显微成像、成像叠加分析与模拟分别标明；每段附标记、原文视频编号、时间尺度、观看重点和物理联系。

使用原期刊媒体服务的远程地址与 HTML5 原生控件，无自动播放，不在仓库中存储或转换视频。`preload="metadata"` 用于加载尺寸、时长及浏览器可用的首帧；点击后播放原媒体。兼容性或网络问题时，可使用每段下方的原文和原始视频链接。

来源、DOI、原文视频号、显示尺寸与检索日期记录在 [paper_video_sources.json](paper_video_sources.json)。专题目录提供逐段跳转入口。视频沿用原始画面与播放内容，实验和模型的判断以原论文测量、参数和图例为准。
