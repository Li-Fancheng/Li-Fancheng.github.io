# 线虫专题

以秀丽隐杆线虫（Caenorhabditis elegans）为模型的研究笔记与文献综述。

胚胎形态发生系列保留原入口 `/notes/celegans/embryonic_morphogenesis/`，分为五篇：

1. [细胞内移与力传递](embryonic_morphogenesis/gastrulation/body.md)
2. [表皮重排与腹侧闭合](embryonic_morphogenesis/epithelial_morphogenesis/body.md)
3. [早期伸长与各向异性力学](embryonic_morphogenesis/early_elongation/body.md)
4. [肌肉驱动与形态记忆](embryonic_morphogenesis/late_elongation/body.md)
5. [三维体形与模型辨识](embryonic_morphogenesis/geometry_and_inference/body.md)

正文通过 `shared-note` 与长篇布局展示，总览使用按 weight 排序的系列布局。各篇有独立目录、局部图号、文献和上一篇／下一篇链接。论文与计算图集中保存在 `embryonic_morphogenesis/figures/`，避免重复素材。

35 幅论文原图保留面板和标注，来源、权利信息、校验值与逐篇图号见 `embryonic_morphogenesis/media_sources.json`；33 条文献元数据见 `literature_sources.json`。显示宽度与原始尺寸见 `data/celegans_figures.json`，网页保留原图链接。

8 幅计算图使用文中明确给出的几何关系或无量纲模型。生成脚本为 `figures/generate_figures.py` 与 `figures/generate_model_figures.py`，需要 NumPy 和 Matplotlib。论文原图不经过这些脚本修改。

各阶段开头嵌入 7 段原始显微影片，第二篇分别展示背侧插入、腹侧包裹和腹侧裂隙封闭；原有激光切割影片保留。影片直接使用研究机构或期刊的远程 MP4，保留播放控制、采样间隔和原始来源，记录见 media_sources.json 中的 embedded_stage_movies。竖向影片限制显示宽度，避免在长篇正文中过度放大。
