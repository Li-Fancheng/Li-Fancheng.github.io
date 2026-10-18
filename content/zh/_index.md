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
        - page: /notes/qbio_course
          label: 课程笔记 · 2026 秋季
          title: 定量生物学导论
          text: 单细胞生长与群体增长、基因调控网络，以及从随机游走到集体聚集的趋化模型。
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
    id: contact
    content:
      title: 联系我
      text: |-
        欢迎交流统计物理、生物建模与 AI for Science 相关的问题，也欢迎指出笔记中的疏漏。

        [电子邮箱](mailto:lifancheng@westlake.edu.cn) · [GitHub](https://github.com/Li-Fancheng) · [个人经历](/experience/)
    design:
      columns: '1'
---
