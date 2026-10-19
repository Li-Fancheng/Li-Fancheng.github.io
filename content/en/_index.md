---
title: ""
type: landing
translationKey: home
summary: "Fancheng Li, PhD student at Westlake University. Research interests and study notes on statistical physics, biological systems, and machine learning."
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
      title: Research Interests
      text: |-
        I am interested in building theoretical models of real systems, drawing on statistical physics, biology, and machine learning. My current interests are exploratory and span three connected directions:

        - **Statistical physics and complex systems.** Stochastic dynamics, nonequilibrium processes, and how microscopic interactions shape collective behavior.
        - **Biological dynamics.** Cell-size regulation, cell-division statistics, population growth, and the physical principles underlying biological organization.
        - **AI for Science & Science for AI.** Using deep learning to study RNA structures and using ideas from physics to understand and improve neural networks.

        My earlier work at Wuhan University included computational condensed matter physics and RNA structure prediction. I hope to bring methods from different disciplines together to answer concrete scientific questions.

        [More about my research and learning](/en/projects/)
    design:
      columns: '1'
  - block: collection
    id: papers
    content:
      title: Featured Publications
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
      title: Recent Publications
      filters:
        folders:
          - publication
        exclude_featured: false
    design:
      view: citation
  - block: collection
    id: talks
    content:
      title: Recent & Upcoming Talks
      filters:
        folders:
          - event
    design:
      view: article-grid
      columns: 1
  - block: note-highlights
    id: notes
    content:
      title: Study & Reading Notes
      text: Course notes, reading notes, and worked derivations on physics, biological systems, and machine learning. Most detailed notes are in Chinese.
      items:
        - page: /notes/qbio_course
          label: Course notes · Fall 2026
          title: Introduction to Quantitative Biology
          text: Notes on single-cell growth and division, population growth, gene regulation, and chemotaxis. The detailed notes are in Chinese.
        - page: /notes/qbio_course/single_cell_to_population_growth
          language: zh
          label: Biological dynamics · 中文
          title: From Single-Cell Variation to Population Growth
          text: Notes in Chinese on cell-size regulation, division-time fluctuations, and correlations across generations, with worked derivations.
        - page: /notes
          language: zh
          label: Reading collection · 中文
          title: More Notes in Chinese
          text: Explore the larger collection on statistical physics, biological systems, and generative models.
      archive:
        page: /notes
        text: Browse all notes →
    design:
      columns: 2
  - block: markdown
    id: thoughts
    content:
      title: Thoughts
      text: |-
        A place for ideas from research and reading, open questions, and reflections beyond the notes.

        **2026.10.18 · [Reflections on Cell Fate and Biological Theory](/thoughts/biological-theory/)**

        How concrete constraints turn broad modeling frameworks into testable explanations. Written in Chinese.

        [Explore thoughts](/en/thoughts/) →
    design:
      columns: '1'
---
