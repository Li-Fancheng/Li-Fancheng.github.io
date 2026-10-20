# Vertex Model：从细胞几何到组织力学

六部分、九篇笔记；网站入口为 `/notes/vertex_model/`。

| 部分 | 源文件 |
| --- | --- |
| I · 模型建立 | [01 几何、能量与力学平衡](foundations/body.md) |
| II · 刚性的几何起源 | [02 线性响应与非仿射弹性](linear_response/body.md)；[03 几何不相容、预应力与刚性转变](rigidity/body.md) |
| III · 从固体到流动 | [04 T1 重排、能量势垒与塑性](rearrangements/body.md)；[05 组织流变与主动动力学](rheology_activity/body.md) |
| IV · 几何推广 | [06 各向异性、曲面与三维组织](geometry_3d/body.md) |
| V · 生物应用与边界 | [07 形态发生、实验检验与开放问题](morphogenesis/body.md) |
| VI · 主动软物质 | [08 主动软物质的堆积与拓扑](active_packing/body.md) |
| VI · 主动软物质 | [09 涨落、相关与组织重塑](fluctuations/body.md) |

页面元信息在 `content/zh/notes/vertex_model/`，正文通过现有 `shared-note` shortcode 渲染，图与正文保存在同一主题目录内。

## 生物案例与物理主线

以果蝇和斑马鱼为主要发育案例。果蝇贯穿翅上皮的排列与受力、蛹翅外部连接切割、胚带定向重排与玫瑰结、腹沟脉冲收缩、卵室基底收缩；斑马鱼保留胚盘与尾部的材料测量、接触角、液滴探针、孔隙与信号反馈。MDCK 绕障碍物流动、人支气管上皮培养与线虫用于比较本构、运动停滞、边界和骨架塑性。第 08 篇比较主动胶体、微管—马达体系与 MDCK 上皮，将堆积、取向缺陷和细胞挤出接回 Vertex Model。

第 09 篇将时间与空间相关、随机驱动的频谱、涨落与响应、重排和密度统计联系起来；生物案例包括果蝇背板与斑马鱼尾部。04、05、08 在相应位置交叉引用。

经典汇合 Vertex Model、非汇合界面模型、接触图、连续流变和运输分别标明自由度与适用范围。原图在论文图来源文件中登记校验值；两段蛹翅视频来自 eLife 正式论文的原始媒体。

## 复现图示

需要 Python 3、NumPy 和 Matplotlib：

```bash
python3 -m venv /tmp/vertex-model-figures
/tmp/vertex-model-figures/bin/pip install -r notes/vertex_model/requirements.txt
/tmp/vertex-model-figures/bin/python notes/vertex_model/generate_figures.py
/tmp/vertex-model-figures/bin/python notes/vertex_model/generate_active_matter_figures.py
/tmp/vertex-model-figures/bin/python notes/vertex_model/generate_fluctuation_figures.py
```

脚本导出 SVG 和相应 300 dpi PNG，数学字体为 Computer Modern，文字优先 Arial。前两个脚本不使用随机采样，包括正多边形几何、真实仿射六边形路径、线性模态、预应力弹簧、激活率分布、Maxwell 响应、自由主动粒子与参考长度可塑性，以及接触角力平衡、铰链约束矩阵、重排活性与方向偏置、孔隙中的溶质守恒等计算。第二个脚本生成主动堆积一篇的减速失稳、六重取向、网络拓扑与向列绕数图。第三个脚本生成相关函数、频谱响应、相关驱动与圆窗数目方差；OU 样本使用固定随机种子 20261020，其余曲线由公式或数值积分计算。各图注明方程、参数及数据来源。

## 论文原图与展示规范

原图文件保持下载时的原始字节，来源、DOI、原文图号、检索日期与 SHA-256 见 [paper_figure_sources.json](paper_figure_sources.json)。图注保留作者与来源；使用作者稿的图像在正文标明版本，包括主动堆积和张力涨落文献的作者稿图。

需要检查或重新获取原图时运行：

```bash
python3 notes/vertex_model/download_paper_figures.py
```

脚本优先检查本地校验值，仅下载缺失文件；远端文件与记录不一致时不覆盖已有图。图片权利与许可按原文保留。

图片的网页宽度与类型记录在仓库 [data/vertex_figures.json](../../data/vertex_figures.json)。宽幅多面板论文图与单一几何示意分别设置展示宽度，点击均可查看原图；使用与 qbio_course 相同的长文、公式滚动、目录和图片展示方式。

## 论文视频

正文嵌入线虫皮层切割、果蝇顶点滑移、蛹翅切割和发育成像、斑马鱼刚性簇图谱、A-PSM 细胞起伏、孔隙与信号响应、胚盘铺展和接触动态，以及主动泡沫模型动画。第 08 篇加入 MDCK 上皮缺陷与细胞挤出的成像分析，第 09 篇加入果蝇背板的局部邻居交换。显微成像、成像叠加分析与模拟分别标明；每段附标记、原文视频编号、时间尺度和物理联系。

使用原期刊媒体服务的远程地址与 HTML5 原生控件，无自动播放，不在仓库中存储或转换视频。`preload="metadata"` 用于加载尺寸、时长及浏览器可用的首帧；点击后播放原媒体。兼容性或网络问题时，可使用每段下方的原文和原始视频链接。

来源、DOI、原文视频号、显示尺寸与检索日期记录在 [paper_video_sources.json](paper_video_sources.json)。专题目录提供逐段跳转入口。视频沿用原始画面与播放内容，实验和模型的判断以原论文测量、参数和图例为准。
