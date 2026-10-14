---
title: ""
date: 2026-10-14
type: landing
cascade:
  reading_time: false
summary: "Fancheng Li — PhD student at Westlake University, exploring statistical physics, biological systems, and AI for science."
design:
  spacing: "4rem"
sections:
  - block: resume-biography-3
    id: about
    content:
      username: admin
      button:
        text: Download CV
        url: /uploads/resume.pdf
    design:
      css_class: dark
      background:
        color: black
        image:
          filename: mybackground.jpg
          filters:
            brightness: 0.7
          size: cover
          position: center
          parallax: false
  - block: markdown
    id: research
    content:
      title: Research Interests
      text: |-
        I am interested in building theoretical models of real systems, drawing on statistical physics, biology, and machine learning. My current interests are exploratory and span three connected directions:

        - **Statistical physics and complex systems.** Stochastic dynamics, nonequilibrium processes, and how microscopic interactions shape collective behavior.
        - **Biological dynamics.** Cell-size regulation, cell-cycle correlations, and the physical principles underlying biological organization.
        - **AI for Science & Science for AI.** Using deep learning to study RNA structures and using ideas from physics to understand and improve neural networks.

        My earlier work at Wuhan University included computational condensed matter physics and RNA structure prediction. I hope to bring methods from different disciplines together to answer concrete scientific questions.
    design:
      columns: '1'
  - block: collection
    id: notes
    content:
      title: Research Notes
      text: "Reading notes and worked derivations on physics, biology, and stochastic dynamics."
      count: 6
      filters:
        folders:
          - notes
      order: desc
    design:
      view: article-grid
      columns: 2
  - block: markdown
    id: contact
    content:
      title: Get in Touch
      text: |-
        I welcome conversations about statistical physics, biological modeling, and AI for science.

        [Email me](mailto:lifancheng@westlake.edu.cn) · [GitHub](https://github.com/Li-Fancheng) · [Browse all notes](/en/notes/) · [中文主页](/)
    design:
      columns: '1'
---
