# Hierarchical Learning 专题重构交付记录

核查日期：2026-10-29。范围：当前仓库本地修改和验证；本任务未执行远程推送或网站发布。

## 1. 修改摘要与审计结论

原有五篇 RHM 正文已经详细覆盖随机语法、同义不变性、样本阈值、语言统计、BP/扩散、解析和共享表示。缺少的是对固定谱模态的训练时间、Fourier 支持可达性、非线性中间表示、有效维数和跨层反馈的成体系解释。已有 `neural_network_theory/hierarchical_learning` 则是预设乘法架构中的可解机制模型，不能代替整个领域的综述。

本次将入口改为以机制组织的 Hierarchical Learning 总论，增加四篇理论正文和一篇综合比较。原五篇 RHM 的 URL、主要推导、图表、参考文献和显式章节锚点保留，调整为独立编号的 RHM 子系列。RHM 01 新增完整第 8 节，讨论 Delrocq 2026 的输入 masking、局部目标、配对信息及数值证据。专题现为总论加十篇正文。

实际源文件根据 wrapper、shared-note shortcode 和网站配置定位，并非从 URL 猜测。未发现适用的 AGENTS.md。沿用现有 Hugo、Markdown、KaTeX、篇内编号、桌面/手机目录、篇末文献与连续导航。未引入项目构建依赖；临时补充 Go 仅用于已有 Hugo 模块加载。

## 2. 必选文献的逐篇覆盖

所有条目均核对了原始正文的模型、算法、关键定义/方程或定理，以及相应适用限制；并非仅根据摘要。全文 PDF 与提取文本放在本机临时研究目录，未加入仓库。

| 指定序号 | 文献与核对版本 | 实质讨论位置 | 核查重点 |
| --- | --- | --- | --- |
| 1 | Saxe, McClelland & Ganguli，2014；arXiv:1312.6120 | [理论 01 §1–3](spectral_dynamics/body.md) | 白化矩阵梯度流、两侧协作、logistic 模态、深度幂次、竞争与初始化 |
| 2 | Saxe, McClelland & Ganguli，PNAS 2019 | [理论 01 §4](spectral_dynamics/body.md#semantic-hierarchy) | 树生成相关谱、粗细分化、隐层 Gram、典型性与语义错觉 |
| 3 | Abbe et al.，NeurIPS 2021 | [理论 02 §2](gradient_accessibility/body.md#staircase) | Fourier staircase 定义、系数条件与逐层随机坐标算法 |
| 4 | Abbe, Boix-Adserà & Misiakiewicz，COLT 2022 | [理论 02 §3](gradient_accessibility/body.md#merged-staircase) | MSP、维数无关动力学、强 O(d) 标度与一般位置限制 |
| 5 | 同上，COLT 2023 | [理论 02 §4–5](gradient_accessibility/body.md#leap-complexity) | Boolean/Gaussian leap、弱信号、在线步数与总运算成本、鞍区逃逸；一般刻画为猜想 |
| 6 | Dandi et al.，LoFi arXiv:2605.13612v2，2026-09-04，预印本 | [理论 04 §1–2](coupled_feature_learning/body.md#forward-lofi) | 局部梯度展开、相关算子、滤波/提升、核递归、有效维数与噪声 |
| 7 | Dandi et al.，BackLoFi arXiv:2610.32849v1，2026-10-23，预印本 | [理论 04 §4–5](coupled_feature_learning/body.md#backlofi-local) | 精确梯度、冻结周围网络的局部定理、配对谱算法、自洽三层替代过程 |
| 8 | Delrocq et al.，arXiv:2605.18557v1，2026-05-18，预印本 | [RHM 01 §8](generative_learning/body.md#local-learning-rules) | DFA/FA 的输入 masking、CLAPP/LPL、配对监督、样本与时间证据 |
| 9 | Allen-Zhu & Li，COLT 2023 扩展摘要；完整 arXiv:2001.04413v6 | [理论 04 §3](coupled_feature_learning/body.md#backward-feature-correction) | 信息间隙、残差代偿、逐次误差修正与规定 SGD 变体 |
| 10 | Dandi et al.，NeurIPS 2025 | [理论 03 §4](nonlinear_depth/body.md#dandi-depth) | SIGHT/MIGHT、增长潜在维数、球面训练/预条件/独立数据、各阶段样本指数 |
| 11 | Saxe, Sodhani & Lewallen，ICML 2022 | [理论 01 §5](spectral_dynamics/body.md#neural-race) | 路径条件统计、共享速率、神经赛跑及门控假设 |
| 12 | Nichani, Damian & Lee，NeurIPS 2023 | [理论 03 §2](nonlinear_depth/body.md#nichani) | 一步梯度提取核投影、非线性读出、特定深浅分离 |
| 13 | Wang, Nichani & Lee，ICLR 2024 | [理论 03 §3](nonlinear_depth/body.md#wang) | 扩散 Hermite 特征、近似 Stein、瓶颈与两阶段 GD；精度余量 |
| 14 | Tabanelli et al.，ICML 2026，PMLR 306 | [理论 03 §5](nonlinear_depth/body.md#hierarchical-spectral) | 随机张量、Hermite 矩阵、谱清理、投影界与恢复条件 |
| 15 | Chen et al.，NeurIPS 2020 | [理论 03 §6.1](nonlinear_depth/body.md#chen-random-representation) | 固定随机表示、白化、二阶 Taylor 模型，非任务驱动表示发现 |
| 16 | Daniely，ICML 2026，PMLR 306 | [理论 03 §6.2](nonlinear_depth/body.md#daniely-label-hierarchy) | 可观察多标签层级、残差网络与逐层学习，非单标签任意电路恢复 |

完整作者、正式发表链接、预印本版本和核对章节记录在 [sources.json](sources.json)，各篇末尾另有可直接访问的原始文献。

## 3. 主要理论判断

1. 必须分别定义生成、目标组合、架构、表示和训练过程的层级；网络层数不规定表示形成的训练顺序。
2. 分阶段曲线至少对应三种机制：已有模态速率不同、先行表示改变梯度可达性、跨层共同演化改写方向/响应/残差。三者有不同干预预测。
3. 非线性中间表示的价值可以来自可组合的低频种子、潜在搜索维数下降，或原始高阶统计在新坐标中变成低阶统计。形成表示的成本和误差不能省去。
4. 上层促进底层学习已有 BFC 与 BackLoFi 等理论实例；新的问题应聚焦促进、延迟、竞争、噪声与稳定性的可预测条件。
5. RHM 让语法、同义类与推断参照可控；其 CNN 匹配结构、标签或正配对也提供知识，不能将其样本阈值直接套成训练时间公式。
6. LoFi 局部推导、硬谱选择算法、受控替代动力学与常规网络实验是不同证据。BackLoFi 的闭合仍含函数态，尚非少数标量的普适 SGD 理论。
7. Chen 的随机表示和 Daniely 的丰富标签需要单独定位。表达、统计上界、计算保证及完整训练轨迹也不可互相替代。

综合篇给出跨文献机制表、三类分阶段机制的干预表，并具体讨论集体变量与 RG 需要满足的条件。

## 4. 保留的开放问题与限制

- 同样的反馈促进何时造成平滑增长、明显延迟、突发或不稳定？极限次序与有限初始化如何影响判断？
- 共享残差、有限容量和冲突任务何时产生真正竞争，何时只改变拟合次序？
- 怎样沿自适应表示轨迹同时控制有限样本算子误差、谱隙、模态切换和时间？
- 哪些模型能够进一步约化到少数重合和响应变量，并对新初始化或新反馈强度作出预测？
- 相近准确率与样本效率是否来自相同内部计算，需要怎样的功能干预？

本次不独立重证所有文献定理，也未重现各论文的神经网络训练实验。正文明确区分证明、受控极限、局部近似、算法设计和数值观察。未将逐层变换视作严格 RG；有效方程只作为明示的研究设想。三篇指定 2026 预印本按上表版本讨论，后续版本或发表状态需重新核对。必选文献均已取得可用原文，无未完成正文页面。

## 5. 文件清单

路径相对于仓库根目录；前两项目录中的 `body.md` 是正文，`index.md` 是本地化元数据 wrapper。

| 类别 | 文件 |
| --- | --- |
| 总论 | `notes/hierarchical_data/body.md`；`content/zh/notes/learning_and_generalization/hierarchical_data/_index.md` |
| 新增理论/综合正文 | `notes/hierarchical_data/{spectral_dynamics,gradient_accessibility,nonlinear_depth,coupled_feature_learning,mechanism_comparison}/body.md` |
| 新增页面 wrapper | `content/zh/notes/learning_and_generalization/hierarchical_data/{spectral_dynamics,gradient_accessibility,nonlinear_depth,coupled_feature_learning,mechanism_comparison}/index.md` |
| 保留并更新的 RHM 正文 | `notes/hierarchical_data/{generative_learning,language_statistics,inference_generation,parsing_latents,shared_representations}/body.md` |
| RHM 元数据与顺序 | 上述五个 slug 对应的 `content/zh/notes/learning_and_generalization/hierarchical_data/<slug>/index.md` |
| 专题元数据 | `notes/hierarchical_data/series.json`、`sources.json`、`verification.json` |
| 已有相关章节的回链 | `notes/neural_network_theory/hierarchical_learning/body.md` |
| 仅作用于综合页宽表的样式 | `assets/css/custom.css` 中 `.hierarchical-comparison-note` 规则 |
| 验证与交付记录 | `notes/hierarchical_data/review_checks.py`；本文件 |

未改动原有 URL 路由表、图片资源、构建配置或项目依赖。工作期间仓库存在其他并行编辑，本任务没有修改或回退其内容；验证使用独立快照并只叠加本专题相关工作树文件。

## 6. 实际验证与复现命令

验证输出汇总见 [verification.json](verification.json)。独立快照与构建结果位于 `/private/tmp/hierarchical-learning-20261029/`。临时目录不属于交付所依赖的项目文件。

### 6.1 Hugo 全站构建

构建快照由当前本地 HEAD 加本专题正文、wrapper、相关回链与专用表格 CSS 组成。使用项目工作流相同的 Hugo Extended 0.136.5，最后一次执行：

```sh
PATH=/private/tmp/hierarchical-learning-20261029/toolchain/go/bin:$PATH \
/private/tmp/fancheng-site-tools/hugo \
  --source /private/tmp/hierarchical-learning-20261029/source \
  --destination /private/tmp/hierarchical-learning-20261029/site \
  --cacheDir /private/tmp/fancheng-site-tools/hugo-cache \
  --noBuildLock --minify
```

通过，exit 0；中文 284 页、英文 29 页。最初缺少 Go 的构建错误已通过临时官方工具链解决，未更改 `go.mod`、`go.sum` 或网站构建配置。本次验证 Hugo 构建，不包括发布流程中的搜索索引部署。

### 6.2 页面、锚点、文献覆盖与公式算例

```sh
python3 notes/hierarchical_data/review_checks.py \
  --site /private/tmp/hierarchical-learning-20261029/site

git diff --check
```

均通过。检查总论及十篇正文、1058 个页面链接项、图片目标、篇内/跨页锚点、全部 16 篇指定文献的覆盖位置、15 条旧入口，以及五篇原 RHM 的显式锚点和正文保留。检查还拒绝数学表达式被误解析为 HTML 标签。1058 项包含外部链接；HTTP 可访问性并非由该静态脚本验证，外部核心文献另经原文检索与下载核查。

对说明性深线性 logistic 解作有限差分残差检查，最大误差约 `7.51e-11`；四项目层级相关矩阵的各特征向量关系通过。与原文方程的人工复核包括 Neural Race 路径系数、BackLoFi 自洽式和 Wang 的精度/样本余量，并渲染相应 PDF 页检查。

### 6.3 KaTeX 与浏览器

使用主题已经包含的 KaTeX 0.16.9，并以 `throwOnError: true` 解析全部专题源码中的数学片段，实际执行：

```sh
/Users/ldepression/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node \
  /private/tmp/hierarchical-learning-20261029/check_math.cjs
```

1224 个片段通过，零解析错误。此检查脚本是临时 QA 文件，不是新增的构建依赖。

通过本地静态服务器在 Codex 浏览器检查全部 11 页：桌面有效宽度 1600、手机宽度 390；逐页核对公式、页面宽度及章节数，并展开全部十篇正文的手机目录。验证了目录的键盘锚点跳转。宽表在专用容器内滚动，页面本身无横向溢出；视口覆盖最后恢复。最终核对源码与浏览器渲染的公式数量一致。

浏览器检查发现并修复新增定义与原语言统计页中 `<` 被当作 HTML 的错误，以及宽表的可访问数学内容逃逸问题。最终版本无已知公式解析、旧 URL、内部锚点或页面溢出问题。
