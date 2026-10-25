# 三维细胞模拟：从相场、网格到 Vertex

五部分、九篇长文，网站入口 `/notes/cell_simulation_3d/`。正文结合 24 张论文图、19 幅计算图和 56 篇文献，沿用 Vertex Model 专题的长文布局、目录、公式滚动和原尺寸图片链接。

| 章节 | 正文 |
| --- | --- |
| 01 · 几何、能量与模型选择 | [foundations/body.md](foundations/body.md) |
| 02 · 多相场 | [phase_field/body.md](phase_field/body.md) |
| 03 · 表面网格 | [surface_mesh/body.md](surface_mesh/body.md) |
| 04 · 三维 Vertex | [vertex_3d/body.md](vertex_3d/body.md) |
| 05 · 格点与 Cellular Potts | [lattice_potts/body.md](lattice_potts/body.md) |
| 06 · 中心粒子与亚细胞元 | [particles/body.md](particles/body.md) |
| 07 · 连续材料与流固耦合 | [continuum_fluids/body.md](continuum_fluids/body.md) |
| 08 · 三维形态发生 | [morphogenesis/body.md](morphogenesis/body.md) |
| 09 · 参数反演与数值验证 | [validation/body.md](validation/body.md) |

`content/zh/notes/morphogenesis/cell_simulation_3d/` 保存页面元信息，正文通过 `shared-note` shortcode 读取。首页与笔记总目录提供入口。专题不依赖新插件或客户端脚本。

## 正文结构

每篇以 `Part / Note` 标记所属部分与篇序，正文收拢为 5—6 个大节，仅在较长推导中保留 2—3 个必要小节，分别使用 `1.` 与 `1.1` 连续编号。论文图与短论点直接接入正文，不再各占一个标题。开头讨论模型思想的来源、变量选择与物理近似，结尾集中讨论本模型的适用范围、失效条件与扩展方向。

第 08 篇以运动学、胚胎约束、定向重排与流动、主动弯曲与增长、液压、机制比较六节组织。核受限迁移移入第 03 篇，球状团增长和守恒输运移入第 06 篇；原位置保留交叉引用。目录额外配对相场与 CPM、表面网格与 Vertex、粒子与连续材料，篇序和公开 URL 不变。

所有既有标题锚点均保留，撤去的标题改为零高度定位点，目录和跨篇深层链接不因合并而失效。每篇的局限讨论可通过 `#model-limits` 直接访问。篇末列出本篇实际引用的参考文献，并提供上一篇、专题目录、下一篇的链接；最后一篇以相关 Vertex 专题作为后续入口。

正文将模型假设、力学推导、实验量与论文讨论放在对应位置。补充内容包括表面应力与接触成熟、相场耗散和接触归一化、网格材料与接触功、三维重连动力学、随机复制的时间映射、多体压缩、有限变形与胞内排液，以及视杯反馈和类器官体积预算。低维教学模型与原论文模型分别注明。

## 生成配图

需要 Python 3、NumPy、Matplotlib：

```bash
python3 -m venv /tmp/cell-simulation-3d-figures
/tmp/cell-simulation-3d-figures/bin/pip install -r notes/cell_simulation_3d/requirements.txt
MPLCONFIGDIR=/tmp/cell-simulation-3d-mpl /tmp/cell-simulation-3d-figures/bin/python notes/cell_simulation_3d/generate_figures.py
MPLCONFIGDIR=/tmp/cell-simulation-3d-mpl /tmp/cell-simulation-3d-figures/bin/python notes/cell_simulation_3d/generate_depth_figures.py
MPLCONFIGDIR=/tmp/cell-simulation-3d-mpl /tmp/cell-simulation-3d-figures/bin/python notes/cell_simulation_3d/generate_mechanism_figures.py
```

三个脚本依次输出共 19 幅 SVG 与对应的 300 dpi PNG，并更新 `data/cell_simulation_3d_figures.json` 中的展示宽度和原始尺寸。第一个脚本保留论文图策略，第二个增加深入算例，第三个增加共同内陷比较；重新生成时按上述顺序运行。无随机抽样。

第一组包含双阱剖面、正八面体细分、棱柱极小值、格点方向代价、Morse 力、本构响应、腔体松弛和合成反演。腔体低维方程使用固定步长 RK4；脚本检查网格几何收敛、相场面积归一化、Morse 力与腔体稳态。

第二组包含三张力平衡、薄界面速度、三角弹簧网的等面积剪切、三维形状指数、CPM 体积能缩放、球形营养分布、持续随机运动 MSD、Maxwell 频率响应、多孔弹性时间标度、薄片屈曲与曲率误差。脚本检查力平衡、剪切系数、体积能单位换算等关系。图注中的数值均为教学参数，公式和假设在对应正文给出。

第三组将同一轴对称内陷模式投影为三个低维动力学模型，比较形成、共同压平探针与关闭来源后的响应；明确区分增长停止和压缩载荷释放。`mechanism_comparison/` 保存逐点 CSV、参数和 RK4 步长减半检查。这个算例没有预测内陷位置与波长，也没有执行相场与显式表面的完整三维跨表示基准。

这些图没有执行完整三维组织动力学，也不使用实验拟合数据。它们用于展示表示、推导与数值基准。正文中概括的教学能量不等同于所链接软件的全部实现。

## 文献与软件

56 篇论文、综述与方法分析，以及官方入口，在各篇相应段落给出，集中登记于 [references.json](references.json)，核对日期为 2026-10-25。记录区分文献类型，并列出实际引用的章节。期刊与卷页信息使用注册书目元数据及原始出版页面核对，元数据来源在相应条目中登记。模型来源部分补充组织重建、弥散界面与 Potts 格点模型的经典工作。机制综合补入 Guirao 与 Etournay 的组织运动学及蛹翅研究、Mongera 的体轴材料分区研究；Guirao 依期刊引用记为 2015，正式版本更新于 2016。Van Liedekerke 等的高分辨率模型区分 2019 在线发表与 2020 卷期年；其他具体书目信息以原文为准。

## 论文图与形态解读

24 张论文图穿插在各篇相应问题旁边，专题首页提供按形态阅读的入口。每幅图注明作者、年份、图号、来源和权利信息，并说明实验成像、分割重建、模拟或力学推断的区别。正文为独立撰写的中文读图说明，不照译原图注。

[paper_figure_sources.json](paper_figure_sources.json) 保存每张图的 DOI、来源链接、下载日期、尺寸、展示宽度、权利信息和 SHA-256。23 张保留发布源图像文件的字节；Okuda 等的网页缩略图清晰度不足，因此从京都大学机构库保存的正式论文 PDF 中完整提取 Fig. 1 的嵌入图像，不裁切、重采样、改色或改动面板，并同时登记源 PDF 校验值与提取方式。恢复这一幅图使用 `pypdf==6.10.0`，本地离线校验不依赖该包。

```bash
python3 notes/cell_simulation_3d/download_paper_figures.py --check-only
# 缺失时从登记来源恢复；文件或源文件校验值不符时停止，不覆盖现有图：
python3 notes/cell_simulation_3d/download_paper_figures.py
```

原图与教学计算分别标为 `paper` 和 `calculation`。页面提供按屏幕尺寸生成的加载版本，点击图片可打开完整尺寸；仓库中保留上述来源文件。没有嵌入远程视频。几何、材料、动力学与算法的适用条件分别在正文说明；软件表用于阅读与选择起点，不作性能排名。
