# 层次结构数据专题

正文保存在本目录的 `body.md` 与五个子目录中。网站入口位于
`content/zh/notes/learning_and_generalization/hierarchical_data/`，通过 `url`
保留 `/notes/hierarchical_data/` 地址。

五篇共享的机制是：边界隐状态将子树内部与外部条件分隔，对外条件分布
因而可由有效状态的后验权重表达；统计分离、有限样本与算法再决定这种
表示能否被辨认并递归复用。总览的机制图由 `mechanism_checks.py` 生成。

## 复现本轮计算

新计算使用 Python、NumPy、Matplotlib 和 Pillow：

```sh
python notes/hierarchical_data/mechanism_checks.py
MPLCONFIGDIR=/tmp/hierarchical-data-mpl python notes/hierarchical_data/bp_mean_field_comparison.py
```

`mechanism_checks.py` 验证有效状态矩阵分解、跨位置传递的互信息与仿射残差，
并生成 `figures/effective_states.svg` 和 `mechanism_checks.json`。

`bp_mean_field_comparison.py` 默认取 `v=16, m=4, s=2`，24 套独立随机字典，
每个字典和深度 16 棵参考树，深度为 2、4、6、8，软证据网格间距为 0.025，
随机种子为 20261025。可用 `--grammars`、`--samples`、`--seed` 改变抽样设置。
其输出为 `figures/bp_mean_field.svg`、PNG 与 JSON；JSON 保留逐字典均值、
字典均值间标准差、均值标准误和闭合误差。

该对照使用与标量闭合相同的参考样本软证据。观测量是参考根状态的后验质量，
不是实际噪声通道中的两次后验采样重合率。根只需向上 BP 即可精确计算；
实现另外用微型树全枚举验证根与向下传播后的叶后验。图中阴影为字典均值间
的一个标准差，其中包含每字典参考树的有限采样误差。

正文中自然语言统计与 Transformer 解析的实验结果来自所引原论文；本轮没有
重新训练这些网络。`mechanism_revision_20261025.json` 记录改动范围与验证结果。
