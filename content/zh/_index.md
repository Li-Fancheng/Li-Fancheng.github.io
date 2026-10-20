---
title: ""
type: landing
translationKey: home
summary: "李凡成，西湖大学交叉科学中心博士生。记录统计物理、生物系统与机器学习的研究兴趣、课程笔记和阅读思考。"
design:
  spacing: "4rem"
sections:
  - block: resume-biography-3
    id: about
    content:
      username: admin
      button:
        text: 下载简历
        url: /uploads/resume.pdf
    design:
      css_class: 'dark personal-intro'
      background:
        color: '#111827'
        image:
          filename: mybackground.jpg
          filters:
            brightness: 0.45
          size: cover
          position: center
          parallax: false
  - block: markdown
    id: research
    content:
      title: '📚 我的研究'
      text: |-
        我关心如何借助物理与数学，为真实系统中的复杂现象建立可分析、可检验的模型。目前的探索围绕三个相互联系的方向展开：

        - **统计物理与复杂系统**：随机动力学、非平衡过程，以及微观相互作用如何形成集体行为。
        - **生物系统的动力学**：单细胞生长与群体增长、基因调控网络，以及生物组织中的自组织现象。
        - **AI for Science 与 Science for AI**：利用深度学习研究 RNA 结构，也尝试用物理思想理解和改进神经网络。

        [了解研究与学习经历](/projects/)

        欢迎交流与合作 😃
    design:
      columns: '1'
  - block: collection
    id: papers
    content:
      title: 精选论文
      filters:
        folders:
          - publication
        featured_only: true
    design:
      view: article-grid
      columns: 2
  - block: collection
    id: recent-publications
    content:
      title: 近期论文
      filters:
        folders:
          - publication
        exclude_featured: false
    design:
      view: citation
  - block: collection
    id: talks
    content:
      title: 近期与即将举行的报告
      filters:
        folders:
          - event
    design:
      view: article-grid
      columns: 1
  - block: note-highlights
    id: notes
    content:
      title: 学习与阅读笔记
      text: 整理课程学习与文献阅读中的模型、推导和问题，记录对统计物理、生物系统与机器学习的理解。
      items:
        - page: /notes/neural_network_theory
          label: 神经网络理论 · 四卷 11 篇
          title: 神经网络理论与可解模型
          text: 从谱动力学与高维泛化出发，推导 NTK、深度线性网络、特征学习、Grokking 与集体动力学，附完整推导和数值图。
        - page: /notes/scaling_laws
          label: Scaling Laws · 六篇专题
          title: Scaling Laws 机制与学习曲线
          text: 从统计估计、几何与谱分辨率，到离散技能、特征学习、ER 和 BA 随机图及计算配置，逐步展开不同标度机制。
        - page: /notes/hierarchical_data
          label: 层次结构数据 · 五篇专题
          title: 层次结构数据与学习理论
          text: 五篇详细笔记，从随机语法、语言统计与扩散推断，到可变树解析和跨语言共享表示，附直观例子、完整推导与计算图。
        - page: /notes/qbio_course
          label: 课程笔记 · 2026 秋季
          title: 定量生物学导论
          text: 单细胞生长与群体增长、基因调控网络，以及从随机游走到集体聚集的趋化模型。
        - page: /notes/vertex_model
          label: 组织力学 · 理论专题
          title: Vertex Model：从细胞几何到组织力学
          text: 7 篇笔记，从能量与几何约束出发，推导组织刚性、细胞重排与流变，并讨论形态发生中的适用边界。
        - page: /notes/celegans
          label: 线虫专题 · 发育与力学
          title: 线虫胚胎形态发生
          text: 从细胞内移、表皮重排与腹侧闭合，到表皮和肌肉共同驱动的身体伸长。
        - page: /notes/ising_course
          label: 统计物理 · 数值计算
          title: Ising 模型的数值模拟
          text: 从 Ising 模型出发，学习统计物理中的数值模拟方法。
        - page: /notes/geometry_of_DGM
          label: 机器学习 · 理论笔记
          title: 生成模型的几何结构
          text: 整理生成模型中的几何观点，并联系扩散模型的相关学习。
      archive:
        page: /notes
        text: 浏览全部笔记 →
    design:
      columns: 2
  - block: markdown
    id: thoughts
    content:
      title: 一些想法
      text: |-
        这里记录研究与阅读中尚未成形的想法、值得追问的问题，以及笔记之外的零散思考。

        **2026.10.18 · [关于细胞命运与生物理论研究的困惑](/thoughts/biological-theory/)**

        宽泛的框架需要具体的约束，才能对细胞行为作出可检验的判断。

        [浏览一些想法](/thoughts/) →
    design:
      columns: '1'
---
