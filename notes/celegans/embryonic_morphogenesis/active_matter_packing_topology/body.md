线虫的细胞需要在特定位置改变邻居，其他连接又必须继续传力；肌肉反复拉动表皮时，一部分形变需要释放，另一部分则要通过骨架与连接重塑保留下来。形态发生所要求的，是这些过程在空间和时间上的配合。

软物质理论为理解这种配合提供了几条不同路线。受限堆积研究边界与界面作用怎样选择细胞接触；主动细胞模型研究局部活动怎样改变轮廓和邻接；连续体与骨架模型则讨论方向性承载、周期负载和形态保持。这些模型的价值，要回到它们实际解释了什么形态变化，以及哪些实验支持了相应机制。

下面先讨论线虫本身的研究，再比较其他体系中可借鉴的模型。后半部分以已发表的细胞插入、力传播与伸长工作为主，通用材料关系只用于解释相关文献。一个过程的细胞在动、连接在更新，并不意味着整片组织已经失去承载能力。

<span id="celegans-studies"></span>

## 1. 线虫研究：早期排列与细胞内移 {#worm-packing-internalization}

先看线虫已有的结果：受限排列与界面作用、内移中的产力和传力、手性流及定向插入、伸长中的材料重塑。细胞质微流变单独说明，因为它测量的尺度与组织重排不同。

<span id="8-线虫中的堆积-主动驱动与材料更新"></span>

### 1.1 四细胞排列怎样同时受卵壳和界面张力控制 {#81-四细胞排列怎样同时受卵壳和界面张力控制}

早期排列是最适合检验颗粒模型的线虫问题之一。细胞数量少，分裂时间、大小和方向可以逐一指定，卵壳形状也能测量。Fickentscher 等用受限软球模型重现早期细胞轨迹，说明几何排斥和边界已经能解释相当一部分有序运动。[Fickentscher 等，2013](https://doi.org/10.1016/j.bpj.2013.09.005) 接下来的问题是：这些条件能否解释胚胎对几何扰动的适应，以及不同细胞接触程度的差异？

Yamamoto 与 Kimura 改变线虫卵壳的长宽比，并比较只有排斥的模型与加入吸引的模型。在去壳胚胎中，AB 的两个子细胞和 EMS 结合较紧，EMS 与 P2 则较松。作者用稳定中心距与两细胞半径之和的比值 $\alpha_{ij}=d_{ij}^{*}/(R_i+R_j)$ 描述这一差异，并把不同接触对的 $\alpha_{ij}$ 放入非均一吸引模型。这里的“非对称”指不同细胞对具有不同相互作用，不是 $i$ 对 $j$ 与 $j$ 对 $i$ 的作用力不相等。加入这些差异后，模型能在更宽的卵壳形状范围内保留 diamond 排列，也能出现只有排斥时遗漏的 T 形排列。[Yamamoto 与 Kimura，2017](https://doi.org/10.1242/dev.154609)

![非均一细胞吸引与卵壳长宽比共同影响线虫四细胞排列](../figures/yamamoto2017_fig5.jpg)

*图 1｜壳形变化与四细胞排列稳定性。引自 [Yamamoto 与 Kimura（2017），Fig. 5](https://pmc.ncbi.nlm.nih.gov/articles/PMC5769634/#DEV154609F5)，CC BY 3.0。A 为模型中的距离与相互作用力；B 改变卵壳长宽比后的排列；C 为不同排列出现的比例。模型对 EMS–P2 与其他细胞对采用不同的稳定距离参数。*

这个结果给“受限堆积”增加了一层含义。卵壳可以把很多排列排除在外，细胞间作用则决定剩下的排列中哪些更稳定。谱系差异因此能通过力学参数进入空间结构。一个模型若只在正常卵壳里得到正确形状，还不能说明它抓住了维持排列的机制；改变壳形以后是否仍然正确，是更有区分力的检验。

不过，$\alpha_{ij}$ 是有效几何参数。中心距离变短，可能来自黏附增强，也可能来自自由表面皮层张力升高，使细胞更愿意扩大相互接触的区域。要区分这两种机制，需要从中心粒子转向界面描述。一个可用于比较的表面能是

<div class="celegans-equation">
$$
E=\sum_i\gamma_{i0}A_{i0}
+\sum_{i\lt j}\gamma_{ij}A_{ij}
+\sum_i\frac{K_V}{2V_i^0}(V_i-V_i^0)^2.
$$
</div>

$A_{i0}$ 是细胞朝向外部介质的面积，$A_{ij}$ 是两细胞的接触面积，$\gamma_{i0}$ 和 $\gamma_{ij}$ 是相应的有效界面张力。最后一项限制体积变化，$K_V$ 具有压力单位；卵壳作为外部几何约束另行加入。对称、无外部压迫的两细胞接触中，若每一侧自由表面的张力为 $\gamma_0$，并把单侧表面与接触平面的夹角定义为 $\theta$，局部力平衡给出 $\gamma_{12}=2\gamma_0\cos\theta$。同一个接触角可以对应不同的绝对张力，因此角度主要提供张力比。

这一层区分已有直接的线虫实验支撑。Yamamoto 等结合 AFM、活体成像和张力反演，研究了两细胞与四细胞阶段。他们把 ABa、ABp、EMS 的紧密结合与自由表面张力增加联系起来，把 P2 初期相对外凸与较高的接触界面张力联系起来；还发现 cadherin 对接触张力既有直接贡献，也会通过调节 actomyosin 间接起作用。[Yamamoto 等，2026](https://doi.org/10.1016/j.devcel.2026.04.013)

![线虫四细胞胚胎的界面张力时序与卵壳内形状预测](../figures/yamamoto2026_fig5.jpg)

*图 2｜从界面张力回到卵壳内的形状。引自 [Yamamoto 等（2026），Fig. 5](https://doi.org/10.1016/j.devcel.2026.04.013)，CC BY 4.0。A 为校准后的张力时序，B、C 比较预测与实测的接触角和界面长度，E 对照显微影像、三维模型与张力分布。图中的张力具有独立实验校准，区别于只给出相对大小的形状反演。*

由此可以看出几类模型之间的关系。球模型用距离势概括细胞接触的结果；表面模型把自由表面和接触界面分开；相场模型再允许轮廓和接触区连续改变。增加自由度的价值，在于能够回答“哪一个界面的性质改变了”。若研究目标只是排列是否稳定，球模型可能已经够用；若要解释 HMR-1 或 myosin 扰动为何改变某一对细胞的接触，就需要把这些界面分别保留下来。

正确的接触结构还参与后续发育。Goldstein 的细胞分离与重组实验表明，P2 与 EMS 的接触参与内胚层诱导；改变接触方位又能改变 EMS 的分裂轴。因此，diamond 排列的生物学意义包括在合适时窗内建立特定细胞接触，而不只是形成一种稳定的外观。[Goldstein，1992](https://doi.org/10.1038/357255a0)；[Goldstein，1995](https://doi.org/10.1083/jcb.129.4.1071)

几何论文与诱导实验回答的是相连但不同的问题：前者给出哪些细胞相接、界面怎样保持，后者说明特定接触能够影响什么发育结果。现有这些结果不能直接等同于“卵壳稍有变形就会改变命运”。更准确的联系是，受限排列为信号交换提供了空间条件，而命运的最终结果还取决于接触时序与信号系统。

### 1.2 原肠胚内移怎样结合堆积能垒与力传递 {#82-原肠胚内移怎样结合堆积能垒与力传递}

Ea、Ep 内移可以沿两条模型路线理解。一条研究细胞数增加、尺寸减小和壳内拥挤怎样改变内外排列的机械代价；另一条研究顶端 actomyosin 怎样通过连接牵动细胞边界。前者解释什么样的移动在几何上有利，后者解释细胞怎样实际获得这种移动。Miao 等的受限堆积预印本展示了仅含各向同性相互作用时出现内移的可能性；Roh-Johnson 等和 Slabodnick 等的实验则表明，皮层已经发生收缩与细胞边界真正开始收缩之间，存在力学耦合的问题。[Miao 等，2021，预印本](https://arxiv.org/abs/2105.05795)；[Roh-Johnson 等，2012](https://doi.org/10.1126/science.1217869)；[Slabodnick 等，2023](https://doi.org/10.1371/journal.pgen.1010319)

Thiels 等的三维定量成像与力学建模预印本，把 Ea、Ep 的皮层流、顶端接触、邻近分裂和胚胎尺度运动放到一起讨论，为这种多细胞视角提供了新的线索。[Thiels 等，2026，预印本](https://doi.org/10.64898/2026.03.30.715391) 另一个有用的对照来自原始生殖细胞：它们可以通过 E-cadherin 介导的连接，随内胚层细胞一起进入内部。[Chihara 与 Nance，2012](https://doi.org/10.1242/dev.079863) 一个细胞发生持续位移，并不要求它本身具有主动粒子模型里的推进项；活动也可以发生在它的邻居身上，再通过接触传递。

在线虫中比较这些解释，需要同时追踪皮层流速、连接边界速度、顶端面积和邻居位移。若皮层仍在流、边界却脱耦，首先指向传力问题；若局部收缩近似保留，改变分裂时序后内移路径和邻居运动显著改变，则更能说明几何与集体松弛的贡献。

## 2. 线虫研究：手性流与定向重排 {#worm-chirality-directed-rearrangement}

### 2.1 手性皮层流怎样把局部活动变成细胞转动 {#83-手性皮层流怎样把局部活动变成细胞转动}

早期线虫还有一种主动过程，很难用“每个细胞获得一个随机推进力”来表达：分裂细胞的两半可以出现方向相反的皮层旋转流。Pimpale 等比较不同谱系，发现这类 counter-rotating flow 与 AB 谱系的分裂轴偏转有关，并通过调节 Rho 通路和 myosin 活动检验二者的关系。[Pimpale 等，2020](https://doi.org/10.7554/eLife.54930)

<figure class="celegans-video celegans-stage-video" id="worm-chiral-cortex" style="--stage-video-width: 34rem;">
<video controls loop preload="metadata" playsinline width="270" height="178" aria-label="线虫AB细胞分裂时两半皮层的反向旋转流">
<source src="https://cdn.ncbi.nlm.nih.gov/pmc/blobs/5b36/7394549/35a054ebf3fe/elife-54930-video2-pmcvs_normal.mp4" type="video/mp4">
浏览器无法播放时，可通过下方原始影片链接观看。
</video>
<figcaption><strong>AB 细胞分裂时的皮层流</strong>｜Pimpale 等（2020），Video 2，CC BY 4.0。NMY-2::GFP 标记 myosin；留意分裂沟两侧的流动方向。这里追踪的是细胞表面的物质流，不是两个子细胞中心绕圈运动。<a href="https://elifesciences.org/articles/54930#video2" target="_blank" rel="noopener">原始影片与说明</a>。</figcaption>
</figure>

这类活动首先产生的是表面切向流和力矩。细胞中心完全不动时，皮层仍然可以持续运动。一个局部薄层近似可以把沿细胞长轴 $x$ 变化的环向流 $v_\theta$ 写成

<div class="celegans-equation">
$$
\eta_s\partial_x^2v_\theta-\Gamma v_\theta
+\partial_x\tau_a(x)=0,
\qquad
\ell_\Gamma=\sqrt{\eta_s/\Gamma}.
$$
</div>

$\eta_s$ 是表面黏度，$\Gamma$ 是皮层对环境的摩擦，$\tau_a$ 表示具有方向性的主动切向应力分量；这里采用使该分量以 $+\partial_x\tau_a$ 进入力平衡的符号约定。主动应力的空间变化驱动流动，黏性把相邻区域耦合起来，摩擦使流动在长度 $\ell_\Gamma$ 上受到屏蔽。真实细胞是曲面，完整理论还要处理曲率、收缩流和活性力矩；这个一维式子只保留它们之间的基本竞争。

![线虫不同细胞谱系的手性皮层流与分裂轴偏转](../figures/pimpale2020_fig3.jpg)

*图 3｜主动流与受限几何怎样共同改变排列。引自 [Pimpale 等（2020），Fig. 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC7394549/#fig3)，CC BY 4.0。A 比较几何限制与非均匀摩擦传动的解释；B、C 为皮层流及谱系相关差异；其余面板比较偏转方向和不同扰动条件。*

皮层怎样带动整个细胞偏转，可以从总力矩平衡看得更清楚。为简化，只考虑细胞绕固定中心的一条轴 $\mathbf e$ 转动。设相对细胞本体的皮层流为 $\mathbf u_c(\mathbf s)$，本体角速度为 $\Omega$，环境暂时静止，摩擦作用于皮层相对环境的总速度。忽略其他外力矩时，

<div class="celegans-equation">
$$
\Omega=
-\frac{\displaystyle
\int_A\Gamma(\mathbf s)
[\mathbf r\times\mathbf u_c(\mathbf s)]\cdot\mathbf e\,dA}
{\displaystyle
\int_A\Gamma(\mathbf s)|\mathbf e\times\mathbf r|^2\,dA}.
$$
</div>

分子是皮层主动流通过摩擦传出的力矩，分母是本体转动的摩擦阻力。在理想对称的反向旋转流中，两半的力矩可能抵消；邻细胞接触和外部包被提供不同摩擦时，抵消就可能被打破。于是，相同的皮层流在不同接触几何中可以导致不同偏转。这也是为什么把摩擦一律当作“越小越容易运动”会失效：摩擦既耗散能量，也提供把表面流转成整体运动的传动条件。

![调节线虫手性皮层流后分裂轴偏转速率的变化](../figures/pimpale2020_fig4.jpg)

*图 4｜把主动流的强弱与偏转速率分开测量。引自 [Pimpale 等（2020），Fig. 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC7394549/#fig4)，CC BY 4.0。A–C 比较 Rho 调控条件下的手性流和朝向分裂沟的收缩流，D 比较皮层流、轴偏转及纺锤体伸长。不同运动分量的测量帮助区分“收缩整体变强”与“旋转分量改变”。*

黏附还会在皮层流中被重新分配。Khor等在两细胞胚胎中观察到HMR-1的定向流与手性分布，并将接触图案同AB收缩环的非对称闭合联系起来。[Khor 等，2025](https://doi.org/10.1016/j.cub.2025.09.052) 这使界面成为一个动态结构：比较细胞作用时，除了距离与接触面积，还需要考虑黏附分子怎样随流动重新分布。它补充了只使用固定两体距离势时未保留的一层信息。

这些结果使线虫中的 active matter 问题比“主动速度多大”更具体：需要问清活动发生在细胞中心还是皮层，产生的是净力还是力矩，哪一部分环境提供摩擦，以及黏附图案是否会被流动改写。手性也不能与后文讨论的向列缺陷混为一谈；一个描述左右不对称的动力学，一个描述方向场绕一圈后的结构。

### 2.2 背侧插入怎样区分定向迁移与组织流动性 {#84-背侧插入怎样区分定向迁移与组织流动性}

背侧插入具有很强的方向和身份信息。两侧细胞形成朝内侧的突起，并逐步取得交替排列。Rac、CDC-42 和 TIAM-1 的研究分别约束突起产生和方向控制。[Walck-Shannon 等，2015](https://doi.org/10.1242/dev.127597)；[Walck-Shannon 等，2016](https://doi.org/10.1371/journal.pgen.1006415)；[Zhu 等，2024](https://doi.org/10.1242/jcs.261509) 对这种过程，随机转向的主动粒子可以作为运动基线，但模型还需要知道细胞应朝哪一侧前进。

![线虫背侧插入中细胞尖端方向和形状可以分别改变](../figures/walck2016_fig2.png)

*图 5｜背侧插入的突起、形状与定向。引自 [Walck-Shannon 等（2016），Fig. 2](https://doi.org/10.1371/journal.pgen.1006415.g002)，CC BY 4.0。A 比较细胞轮廓，B 比较尖端方向，C 比较长宽比；这些量分别约束形变与极性，不能由同一个“活动强度”代替。*

这里还有一个容易在顶点模型里丢失的过程：细胞可以先在基底侧伸出突起，再发生顶端接触更新。若用单层多边形把整个细胞压成一个轮廓，模型中的 T1 时间可能对应后一个事件，驱动却来自前一个事件。因此，应分别记录基底前缘开始前进的时间、顶端旧接触缩短的时间，以及新接触能够持续承载的时间。它们之间的延迟能约束顶底耦合、细胞形变与连接成熟。

由此可以区分两种常被合在一起的解释。一种是组织的重排阻力降低，同样的定向驱动更容易完成插入；另一种是组织阻力近似保留，细胞的方向控制或牵引能力改变。两者都能影响插入速度。前者还应改变对小幅外力或局部释放的响应，后者则可以首先改变突起方向、运动持续性或皮层与边界的速度差。

接触消失也必须按事件分类。线虫表皮会通过 EFF-1 介导的融合形成合胞体，因此在较长的发育时间窗里，细胞节点本身并不守恒。[Mohler 等，2002](<https://doi.org/10.1016/S1534-5807(02)00129-6>)

| 观察到的变化 | 空间或拓扑含义 | 对模型的处理 |
| --- | --- | --- |
| 顶端旧接触消失，基底接触仍在 | 某一层的邻接改变，三维细胞未必分离 | 分层接触图或三维界面 |
| 四个细胞完成旧邻居与新邻居交换 | 可与 T1 过程比较 | 记录方向、持续时间与逆转 |
| 两个细胞之间的膜消失，细胞质连通 | 细胞融合，两个节点合为一个 | 更新节点及原有界面 |
| 一个母细胞产生两个子细胞 | 分裂引入新节点和新接触 | 按谱系建立对应 |

若把后三种事件全部统计成邻居保持率下降，就很容易得到“组织越来越流动”的结论，却不知道细胞究竟通过哪种方式改变结构。在线虫中，谱系、膜标记与融合信息可以帮助把这些事件分开。

### 2.3 腹侧包裹怎样从自由边界推进变成连接承载 {#85-腹侧包裹怎样从自由边界推进变成连接承载}

腹侧包裹中的活动发生在一片尚有开口的表皮前缘。Leading cells 先在中线相遇，后方 pocket 随后封闭；前缘下方的神经母细胞还提供一个会变形的接触环境。经典研究区分了 actin 依赖的推进和后续封闭，cadherin 相关实验进一步说明，初次相遇与稳定封口之间还存在连接建立的过程。[Williams-Masson 等，1997](https://doi.org/10.1242/dev.124.15.2889)；[Raich 等，1999](<https://doi.org/10.1016/S0960-9822(00)80015-9>)

![线虫腹侧包裹的前缘突起和F actin分布](../figures/bernadskaya2012_fig2.png)

*图 6｜从自由前缘到封闭表面的动态过程。引自 [Bernadskaya 等（2012），Fig. 2](https://doi.org/10.1371/journal.pgen.1002863.g002)，CC BY 4.0。A 为 actin 报告系统显示的前缘行为和局部放大，B 比较相关事件的时间间隔。前缘位置、突起分布与事件时间可以分别约束推进和封口模型。*

拓扑的改变也与 T1 不同：先前暴露于开口的两段自由边界，转成相互连接的界面，边界条件随之改变。将整个过程称为组织流化，会漏掉封口之后建立承载路径这一步。腹侧包裹因而需要同时讨论迁移、曲率、底层组织运动和连接成熟，单独增加主动粒子的速度很难代替这些机制。

包裹完成后，表皮从有自由边缘的细胞层变为连续覆盖内部组织的表面。新连接随后要承受伸长的牵拉，因而“覆盖到中线”和“形成能共同承载的表皮”是两个相继的生物学结果。线虫的迁移与连接实验为这一区分提供了直接依据；其他胚胎的闭合模型，则进一步说明材料更新怎样与连续覆盖同时发生。

## 3. 线虫研究：伸长与材料更新 {#worm-elongation-material-remodelling}

### 3.1 伸长中的 actin 取向怎样进入主动固体模型 {#86-伸长中的-actin-取向怎样进入主动固体模型}

到了早期伸长，主要问题开始从“细胞怎样更换邻居”转向“不同区域怎样分担形变”。侧表皮的收缩、背腹表皮的 actin 排列，以及连接和外部结构的传力共同作用。Vuong-Brender 等通过方向相关的激光切割分析了应力与材料各向异性，Fang 等则把主动收缩和纤维取向变化放入连续体模型。[Vuong-Brender 等，2017](https://doi.org/10.7554/eLife.23866)；[Fang 等，2021](https://doi.org/10.1126/sciadv.abg3264)

向列张量在这里仍然有用，但需要明确它代表哪一种方向。由细胞长轴计算的 $\mathsf Q_{\rm cell}$，由 actin 纤维计算的 $\mathsf Q_{\rm actin}$，以及由瞬时速度计算的取向，可以给出不同结果。细胞整体沿身体轴变长，并不要求其中每根承载纤维也沿同一方向排列。比较它们是否耦合，可以用局部方向差的 $\langle\cos2(\theta_{\rm cell}-\theta_{\rm actin})\rangle$；接近 1 表示平行，接近 −1 表示垂直，接近 0 表示该平均下没有明显偏好。

### 3.2 肌肉周期拉动怎样转成持续伸长 {#87-肌肉周期拉动怎样转成持续伸长}

后期伸长更接近一个被周期驱动的复合材料问题。肌肉缩短时向表皮传递负载，表皮与连接既要承受短时形变，又要在多次循环之后改变较长期的形状。若每个细胞都只获得持续随机推力，模型会漏掉驱动的节律、肌肉位置以及表皮的加载与卸载路径。

这种区分与线虫实验直接相关。Lardennois 等研究的 actin 黏塑性锁定涉及组织如何保留肌肉引起的形变，相关通路受损时可出现伸长后回缩。Fang 等把 actin 切断与重新束化的速率竞争写进连续体模型；Yang 等进一步研究了周期机械输入、表皮连接更新和伸长之间的关系。[Lardennois 等，2019](https://doi.org/10.1038/s41586-019-1509-4)；[Fang 等，2021](https://doi.org/10.1126/sciadv.abg3264)；[Yang 等，2025](https://doi.org/10.1016/j.devcel.2025.06.012)

![actin切断与重新束化速率共同决定线虫晚期伸长的模型结果](../figures/fang2021_fig4.jpg)

*图 7｜材料更新也可能破坏承载。引自 [Fang 等（2021），Fig. 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC8245039/#F4)，CC BY-NC 4.0。A 为有无 actin 塑性的模型比较，B、C 改变切断与重新束化速率。过快拆解而来不及重建，与适度更新产生的持续伸长，具有不同后果；图中曲线为模型输出。*

## 4. 线虫研究：力反演与细胞质力学 {#worm-force-inference-microrheology}

### 4.1 怎样从线虫轨迹得到有效相互作用和力网络 {#88-怎样从线虫轨迹得到有效相互作用和力网络}

<span id="有效力反演的动力学选择"></span>

前面的模型往往从势能出发预测轨迹，也可以反过来由轨迹约束相互作用。Koyama 等利用细胞核追踪，在有效中心力与过阻尼动力学的假设下反演细胞间作用，并将方法用于线虫胚胎。他们获得了含排斥区和吸引区的有效距离势，比较不同阶段及归一化后的曲线，再用这些作用进行形状模拟。[Koyama 等，2023](https://doi.org/10.1371/journal.pcbi.1011306)

![由线虫细胞核轨迹反演有效细胞间力和距离势](../figures/koyama2023_fig4.png)

*图 8｜轨迹如何约束粗粒化模型。引自 [Koyama 等（2023），Fig. 4](https://doi.org/10.1371/journal.pcbi.1011306.g004)，CC BY 4.0。A、B 为谱系轨迹及估计尺寸，C 为推断的有效力网络，D 检查给定方法中不同初值的结果，E–G 为距离—力、距离—势及归一化比较。力采用任意单位，不能直接与 AFM 校准的张力相加。*

### 4.2 线虫细胞质的微流变为何不能直接代表组织物态 {#97-线虫细胞质的微流变为何不能直接代表组织物态}

线虫中已经可以测量发育过程中的细胞质力学。Koizumi 等的预印本使用微米级磁性液滴做主动微流变，同时用 40 nm 荧光微珠追踪局部运动。作者观察到早期发育中细胞质蠕变柔量逐渐降低，抑制胞质分裂后这一趋势仍存在，因而不能完全归因于细胞越分越小。[Koizumi 等，2026，预印本](https://doi.org/10.64898/2026.05.19.726147)

![线虫早期胚胎细胞质的主动微流变蠕变柔量与卸载恢复](../figures/koizumi2026_fig3.jpg)

*图 9｜线虫中对材料响应的直接测量。引自 [Koizumi 等（2026），Fig. 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC13228458/#F3)，预印本，CC BY-NC-ND 4.0。A 为随加载时间变化的柔量，B 比较 10 s 时的柔量，C 为拟合的时间指数，D 为卸载恢复比例。原图完整保留。*

这里的磁滴在细胞内部运动，主要探测局部细胞质对受控力的响应；Mongera 等的微滴位于细胞群之间，探测周围组织的重排与承载。相似的磁性探针，因所在位置不同，测到的是不同层次的材料。

这一差别可以用一组并不矛盾的变化说明：细胞内部变得更难形变，细胞之间的连接却更快交换，组织仍可能更快释放剪切应力；反过来，细胞质仍较黏性，持久连接和方向性骨架也可以让多细胞网络在发育时间尺度内承载。细胞质、皮层、连接与整个细胞群需要分别定义响应。柔量下降说明给定条件下更难被推动，并不等于已经观察到细胞质或组织的固化相变。

<span id="transferable-theories"></span>

## 5. 理论借鉴：受限堆积怎样选择细胞接触 {#packing-confinement-rigidity}

### 5.1 卵壳、界面作用与接触结构 {#1-从胚胎堆积到软物质问题}

<span id="11-拥挤-刚性与运动停滞怎样区分"></span>
<span id="12-卵壳不仅限制体积-还选择接触排列"></span>
<span id="13-从位置到形状-再到取向"></span>

早期排列的关键输出，是哪些细胞接触，以及这种接触在分裂和边界约束下能否保持。Giammona 与 Campàs 比较不同包裹条件中的四细胞排列，将细胞大小、分裂和界面作用放入同一几何问题。狭长的外壳会排除一部分自由空间中的排列，却不能独自决定所有接触；细胞对之间不同的界面性质仍然重要。[Giammona 与 Campàs，2021](https://doi.org/10.1371/journal.pcbi.1007994)

![早期胚胎的细胞排列界面作用和包裹约束](../figures/giammona2021_fig1.png)

*图 10｜包裹边界与细胞相互作用共同限制排列。引自 Giammona 与 Campàs（2021），[Fig. 1](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007994#pcbi-1007994-g001)，CC BY 4.0。A 比较不同物种，B–F 定义模型中的细胞作用和边界，G–I 涉及空间分区与分裂。图中的模型服务于接触排列问题。*

细胞球模型用中心、大小及相互作用描述这种选择。它能比较分裂方向、壳形和吸引强度改变后出现的邻接；三维表面模型则进一步保留界面面积与形状，可以与接触角和张力测量对照。两者的区别是保留了哪些几何量，不能只按“简单”和“复杂”排序。

这里也要区分整体稳定与关键接触稳定。两个排列的长宽比相近，P2–EMS 接触却可能不同；反过来，整体轮廓有所变化，关键接触也可能继续保持。把接触结构与前述命运诱导实验联系起来，比只比较轮廓更能说明堆积力学的发育意义。

### 5.2 刚性与重排阻碍来自什么 {#2-被动软颗粒的刚性从哪里来}

<span id="21-接触力来自重叠能量"></span>
<span id="22-接触数怎样约束运动"></span>
<span id="23-为什么同样的网络也可能有不同刚度"></span>

对给定排列，软颗粒模型从接触的变形代价求回复力。若 $U$ 是相互作用能，力为 $-\nabla U$，稳定点附近的二阶导数矩阵描述哪些集体形变较难、哪些较容易。接触很多不一定使所有形变都困难：一些连接允许转动，某些集体运动仍可比较柔软。[O’Hern 等，2003](https://doi.org/10.1103/PhysRevE.68.011306)

![接触网络约束Hessian与刚性判据的示意](../figures/20_packing_rigidity.svg)

*图 11｜接触、约束与刚性并非同一个量。四节点杆网络从正方形、增加一条斜杆到增加两条斜杆。第一条斜杆限制原有剪切软模，第二条引入冗余约束；连接数、柔软运动与自应力需要分别记录。此图是约束计数示意，不是胚胎测量。*

对线虫，周围细胞能否调整，是这套理论最有用的联系。一个内移细胞面对固定邻居，与面对可以形变和移动的邻居，会有不同的阻碍。分裂还会更新接触与细胞大小，使此前的稳定构型发生变化。因此，局部收缩相同，也可能得到不同的内移时序。

这些结论解释的是给定材料附近的稳定性，尚未规定细胞怎样主动跨过重排路径。进入表皮重排以后，细胞形状、接触和活动方向必须继续作为动态变量保留。

## 6. 理论借鉴：主动细胞模型保留哪些空间自由度 {#active-particles-confluent-tissues}

### 6.1 主动粒子描述驱动与集体运动 {#3-主动驱动怎样改变笼困与重排}

<span id="31-主动软粒子的最小动力学"></span>
<span id="32-自由扩散与笼中运动的解析比较"></span>
<span id="33-从单一笼子到集体软模态"></span>
<span id="34-主动温度-mips-与解堵分别解释什么"></span>

高密度主动粒子研究表明，即使接触作用不变，持续的自驱动也会改变笼困、集体运动和重排。一个常用关系是 $\dot{\mathbf r}_i=-\mu\nabla_iU+v_0\mathbf p_i$：位置随相互作用和主动方向共同变化。活动的大小与持续时间是两个不同因素；相同的平均速度，可以对应短暂、分散的推动，也可以对应持续作用于同一路径的负载。[Henkes 等，2011](https://doi.org/10.1103/PhysRevE.84.040301)

在线虫中，这有助于区分主动驱动与被邻居带动，但粒子中心轨迹没有保留前缘突起、细胞侧面和顶端连接。它适合回答“活动怎样推动拥挤网络”，不足以单独解释背侧细胞怎样伸出楔形尖端、穿过邻居并完成插入。自由扩散和谐势笼困的解析关系描述时间尺度；形态怎样改变，还取决于轮廓和接触。

### 6.2 顶点与 Voronoi 模型保留形状和邻接 {#4-填满空间的细胞怎样改变刚性}

<span id="41-面积与周长之间的几何不相容"></span>
<span id="42-主动-voronoi-模型怎样把运动加进来"></span>
<span id="43-顶点和-voronoi-的刚性结论为何不能直接互换"></span>
<span id="44-三维堆积怎样增加新的形状和拓扑变量"></span>

在满铺组织中，细胞之间几乎没有可用空隙，运动需要通过形状改变与接触更新完成。顶点模型把细胞表示为共享边界的多边形，面积与周长通常进入能量

<div class="celegans-equation">
$$
E=\sum_i\left[\frac{K_A}{2}(A_i-A_{0i})^2
+\frac{K_P}{2}(P_i-P_{0i})^2\right].
$$
</div>

这个写法规定形状改变的代价；顶点的运动、边收缩与T1规则，才决定接触怎样更新。SPV模型则以细胞中心生成Voronoi分区，再用分区的面积、周长和主动运动计算演化。它保留了满铺与邻接，但轮廓受Voronoi构造限制，不能像独立皮层模型那样直接描述两侧膜的滑动。[Bi 等，2015](https://doi.org/10.1038/nphys3471)；[Bi 等，2016](https://doi.org/10.1103/PhysRevX.6.021011)

![主动Voronoi细胞模型的几何与运动状态](../figures/bi2016_fig1.jpg)

*图 12｜细胞轮廓与轨迹在同一空间模型中演化。引自 Bi 等（2016），[Fig. 1](https://doi.org/10.1103/PhysRevX.6.021011)，CC BY 3.0。模型中的位置、形状和邻接彼此联系，和仅给单个方向角写动力学方程有不同的描述能力。*

![主动Voronoi模型中的驱动目标形状及持续性依赖](../figures/bi2016_fig3.jpg)

*图 13｜活动与目标形状共同改变动力学。引自 Bi 等（2016），[Fig. 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC5619672/#F3)，CC BY 3.0。各参数切片及三维示意采用该文的运动学判据，不能直接作为线虫胚胎的材料相图。*

目标形状 $p_0=P_0/\sqrt{A_0}$ 是模型参数，实际形状指数 $P/\sqrt A$ 是组织已经变成怎样的测量。外部拉伸也能提高实际形状指数，因而不能把观察值直接代入一个通用阈值。后续研究还指出，顶点与Voronoi模型即使使用相似能量，静态刚性的结论也不完全相同。[Sussman 与 Merkel，2018](https://doi.org/10.1039/c7sm02127e)

三维模型另需保留细胞体积、界面面积及不同高度的接触。Scutoid等几何研究说明，顶面与底面的邻居可以不同；但几何上容许这种结构，并没有同时解释驱动来自哪里。[Merkel 与 Manning，2018](https://doi.org/10.1088/1367-2630/aaaa13)；[Gómez-Gálvez 等，2018](https://doi.org/10.1038/s41467-018-05376-1)

### 6.3 定向T1怎样产生组织形变 {#5-邻接网络的拓扑怎样参与组织重排}

<span id="51-t1-改变的是哪几条接触"></span>
<span id="52-欧拉关系为什么使平均边数接近六"></span>
<span id="53-重排率还需要能垒与驱动"></span>
<span id="54-三维邻接图需要记录接触出现在哪一层"></span>

一次平面T1使一对细胞失去旧接触，另一对细胞取得新接触。它可以释放或积累剪切形变；多次事件能否共同改变组织形状，还取决于方向、空间次序与是否逆转。

![四个细胞的一次T1邻接交换以及边数变化](../figures/22_t1_topology.svg)

*图 14｜T1改变接触图。A–B旧接触消失，C–D新接触形成，外侧接触保留。图中的连线是邻接关系，不是细胞膜，也未指定实际重排由哪种分子活动驱动。*

二维铺砌中的平均边数受到欧拉关系约束，三叉顶点的无边界周期铺砌平均为六；这类统计本身没有给出边张力与重排时间。对正在融合、分裂或具有自由边缘的线虫表皮，还必须分别处理这些拓扑变化。

背侧插入尤其需要保留事件发生在哪一层。基底突起已进入邻居之间时，顶端旧接触仍可能存在。一个平面T1把两者合为同一步，因此可概括最终邻接变化，却不能解释二者之间的时差。独立皮层与多细胞模型进一步保留了轮廓和接触，能讨论方向角方程中没有包含的传力过程。

## 7. 理论借鉴：组织怎样加载、松弛与继续承载 {#solid-fluid-transition}

### 7.1 材料性质需要对应具体的受力过程 {#91-固态和流态究竟在区分什么}

<span id="93-刚性丧失与活动驱动的解笼怎样区分"></span>

固流讨论首先针对组织对扰动的响应。小幅形变引起多大回复、应力在多久内释放，以及较大负载是否引起不可逆重排，分别对应刚度、松弛与屈服。它们可以共同改变，也可以在同一组织中彼此分开。

对线虫，背侧插入与后期伸长给出了两种很不同的使用情形。前者需要邻接沿特定方向更新，后者需要表皮在肌肉拉动时承载、在周期之间更新部分结构。用一个标量“流动性”同时描述它们，会丢掉方向、位置和负载条件。

| 文献中的判据 | 对应的操作与观测 | 与形态发生的联系 |
| --- | --- | --- |
| 小形变刚度 | 在当前状态附近施加小变形，测量额外回复 | 哪些连接或区域仍能传递负载 |
| 应力松弛 | 快速变形后保持形状，观察额外应力下降 | 该时间窗内能否释放已有的机械记忆 |
| 蠕变与屈服 | 保持外力，比较有限变形与持续重排 | 新邻接是否需要超过一定负载才出现 |
| 结构更新 | 追踪接触的形成、消失、方向和逆转 | 哪些局部重排真正贡献了组织形变 |

相同的细胞位移可以出现在不同加载条件下。无外力时的主动爬行、被邻居拉动时的位移，以及外力撤除后的回弹，各自反映不同的机制。因而，主动组织的流变测量通常需要从同阶段的背景活动中分离额外响应，而不能用自然发育速度替代受力实验。

在顶点与主动Voronoi文献中，这一区分同样影响“解堵”的含义。改变材料参数可以改变某个构型的回复与重排代价；保持材料关系不变、增加持续驱动，也能使细胞跨过已有阻碍。两条路线在运动状态上可以相似，物理上却分别涉及能量地形和驱动。[Bi 等，2015](https://doi.org/10.1038/nphys3471)；[Bi 等，2016](https://doi.org/10.1103/PhysRevX.6.021011)

### 7.2 短时承载与长时响应 {#92-用应力松弛-蠕变和频率响应判断材料}

松弛模量 $G(t)=\delta\sigma(t)/\gamma_0$ 表示在小应变阶跃 $\gamma_0$ 后，额外应力怎样变化。Maxwell模型给出 $G(t)=G_0e^{-t/\tau_R}$；加入一条长期弹性路径后，可写成 $G(t)=G_\infty+G_re^{-t/\tau_R}$。前者最终释放该方向的应力，后者在模型中保留回复平台。

![固态长松弛流体与普通流体的材料响应比较](../figures/27_solid_fluid_rheology.svg)

*图 15｜同一观察窗口中的不同材料响应。标准线性材料关系的综合示意：A为松弛，B为频率响应，C为恒应力蠕变。取 $G_0=1$；两种Maxwell流体的 $\tau_R=1,100$，固体取 $G_\infty=0.3,G_r=0.7,\tau_R=1$。此图解释判据，不是线虫材料参数的拟合。*

松弛很慢的流体可以在一次较快的形态变化中维持负载；一片长期仍能承载的材料，也可以通过少数局部连接更新改变形状。因此，材料名称需要同时保留方向、时间窗和加载方式。频率响应中的储能与损耗交点，首先提供时间尺度，不能单凭这个交点判定发育中发生了相变。

Kelvin–Voigt关系 $\sigma=E\varepsilon+\eta\dot\varepsilon$ 中的 $\eta/E$ 是恒应力下应变建立的迟滞时间，而Maxwell关系中的 $\eta/G$ 是应力松弛时间。两者虽有相同量纲，描述的加载实验不同。分子驻留或FRAP恢复时间则描述组分交换，也不能直接替代任一机械时间。

蠕变实验补上另一种信息：外力保持时，材料是逐渐达到有限形变，还是通过结构更新持续变形。加载—卸载的回线则记录耗散，但慢回弹也会留下短时的剩余形变，不能把每一条回线都当作永久重塑。文献中的“残余”因此要连同卸载等待时间一起理解。

激光切割又有不同的控制方式。它释放已有负载，最先读到的回缩涉及切割前应力和局部阻力；后续开口涉及材料及周围约束。把切割、自然形变和分子交换资料放在一起，能够约束各自没有直接测到的量，但同一个初速度不能同时独立给出绝对张力、刚度和黏性。这也是下列组织研究结合多类观测的原因。

### 7.3 闭合组织可以更新形状而继续传力 {#94-屈服之后流动与固流转变有什么关系}

Tah等研究果蝇背侧闭合时，面对的正是“明显改变形状却不破坏连续承载”的组织。Amnioserosa细胞随开口缩小而变窄，形状指数升高；若把参考周长固定，顶点模型会走向不同的刚性状态。允许参考周长随发育更新，再保留细胞间差异，模型可以在闭合过程中维持非零剪切模量。[Tah 等，2025](https://doi.org/10.1073/pnas.2322732121)

![闭合中的细胞形状与参考状态更新](../figures/tah2025_fig1.jpg)

*图 16｜形状变化与模型状态同时记录。引自 Tah 等（2025），[Fig. 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC11725931/#fig01)，CC BY-NC-ND 4.0，完整原图未修改。A、B比较实验和模型的空间形变；其他面板记录实际形状与模型参考周长。*

![闭合过程中的刚性和参考周长更新](../figures/tah2025_fig3.jpg)

*图 17｜参考状态影响闭合中的承载。引自 Tah 等（2025），[Fig. 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11725931/#fig03)，CC BY-NC-ND 4.0，完整原图未修改。比较固定与更新参考周长时的模型刚性及形态轨迹。*

这项工作的价值是将当前形状与材料偏好的形状分开，并把它们放回真实闭合过程。论文还用切割资料检验模型给出的张力变化。它为线虫表皮提供的是材料更新与承载可以共存的实例；果蝇amnioserosa的细胞行为和线虫表皮包裹并不相同，不能直接转用其参数。

这一解释还有一个可与实验比较的结果：张力在闭合中不是简单地持续升高或降低。原文结合模型、连接直度与切割初始回缩，比较了非单调的张力变化，因而不只是在终点重现一个较窄的开口。

![闭合顶点模型与激光切割及连接几何的比较](../figures/tah2025_fig2.jpg)

*图 18｜材料模型接受形状以外的检验。引自 Tah 等（2025），[Fig. 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC11725931/#fig02)，CC BY-NC-ND 4.0，完整原图。A–D为形状与方向统计，E、G将模型中的张力变化同切割初始回缩及连接直度比较。这里检验的是变化趋势，绝对张力仍涉及耗散标定。*

Fielding等的连续本构模型从另一侧分开刚度与塑性松弛：周长和取向参与应力，T1相关的内部变量控制松弛。与受剪切顶点模型的比较说明，有限变形、初始小形变刚度与动态屈服可以有不同依赖。[Fielding 等，2023](https://doi.org/10.1103/PhysRevE.108.L042602) 这里重要的是分别记录材料怎样承载、又怎样更新，不必给每一种响应重复建立一个固流边界。

这种本构层次的区分，有助于理解同一组织为什么先储存形变、随后重排。拉伸或剪切开始时，细胞形状与取向先改变，已有网络可以变得更难继续变形；随着T1增加，部分应力又得到释放。Fielding等比较的正是启动后的应力—应变曲线和持续加载下的流动曲线，两者不必由一个固定松弛时间概括。它是对顶点模型流变的总结，不是线虫某个发育阶段的直接测量。

### 7.4 直接施力怎样把材料状态定位到胚胎区域 {#95-其他胚胎怎样建立固流转变的证据}

Mongera等将斑马鱼体轴中的细胞运动与磁性微滴加载结合。不同区域对微滴形变的回复、屈服和残余形变不同，容易重排的尾端区域与较前方提供支撑的区域因此有了直接力学对应。[Mongera 等，2018](https://doi.org/10.1038/s41586-018-0479-2)

![斑马鱼不同体轴区域的受控加载与卸载响应](../figures/mongera2018_fig2.png)

*图 19｜材料性质落实到位置与加载历史。引自 Mongera 等（2018），[Fig. 2](https://www.nature.com/articles/s41586-018-0479-2/figures/2)，完整原图未修改，原出版方权利保留。a–c包含加载、卸载及残余形变；其他面板比较位置与黏附干预。*

Petridou等则将接触网络中的刚性簇与胚盘流变测量对照，说明几何连通和机械贯通需要分别讨论。[Petridou 等，2021](https://doi.org/10.1016/j.cell.2021.02.017) 两类研究共同推进的一步，是把组织某处的细胞运动、接触结构与受力响应放在同一证据链里。它们没有给出一条可直接套在线虫全身的区域划分。

接触网络研究还区分了“细胞仍相连”和“刚性贯穿组织”。相互连接的细胞群可能由可以彼此转动的较硬小块组成；网络没有断开，整体仍可以容易形变。Petridou等将这类结构分析与胚胎材料响应及黏附扰动对照，使刚性簇成为可以联系到流变的对象，而非仅将邻居数量作统计。

这些斑马鱼结果对线虫的意义，在于展示证据怎样从细胞接触走到材料行为。线虫细胞数、谱系相关接触、卵壳与方向性表皮都不同，不能按身体位置直接给它套上同一张分区图。线虫已有的工作更适合先沿具体运动来读：该处需要什么结构变化，什么组织提供驱动，哪些路径继续承担负载。

### 7.5 线虫文献中的材料分工 {#96-在线虫各阶段应该怎样提出固流问题}

<span id="98-怎样区分有限时间交叉与真正的转变"></span>

对背侧插入，直接证据首先定位了突起、极性与邻接的改变；对早期伸长，分区切割与纤维资料定位了应力和刚度各向异性；对肌肉参与后的伸长，回缩表型、骨架重塑和连接交换资料定位了形态保持中的不同环节。[Walck-Shannon 等，2015](https://doi.org/10.1242/dev.127597)；[Vuong-Brender 等，2017](https://doi.org/10.7554/eLife.23866)；[Lardennois 等，2019](https://doi.org/10.1038/s41586-019-1509-4)

腹侧包裹把这种分工表现为一次边界条件改变。前缘推进时，细胞需要同下方表面和邻居发生相对运动；到达中线后，新的界面又要开始承受伸长的牵拉。HMR-1相关实验中，接近、接触与稳定连接可以分别受损，因而相同的终点开口未必来自相同的材料问题。[Raich 等，1999](https://doi.org/10.1016/S0960-9822(00)80015-9)

在后期伸长中，材料状态的后果则体现为回缩与继续伸长的区别。每轮最大形变反映当时的输入和承载，周期之间保留下来的长度还依赖结构更新。现有骨架与连接研究用不同资料建立这两层联系，不能只用“是否动得快”把正常伸长和受损后的回缩排在同一条轴上。以下空间模型综述将分别说明这些文献实际保留了哪些位置、结构与时序。

这些资料支持的是具体过程中的局部材料分工。某些接触需要让出重排路径，其他接触仍然承担负载；某些应力需要松弛，另一些方向性结构需要保留。有限胚胎又在持续发育，因而从短影片中的运动状态推到一个严格的长期相态，通常超出了现有测量。

## 8. 空间模型综述：插入、力传播与身体伸长 {#面向线虫的候选计算与检验}

### 8.1 线虫研究 内移中的局部驱动与邻居响应 {#堆积能垒与-clutch-的耦合算例}

Ea、Ep内移要让细胞失去表面暴露，同时使邻居覆盖原来的位置。皮层已有收缩而边界尚未跟随的实验，说明活动能否传到边界是一个环节；周围细胞是否让位，则是另一个环节。两者共同影响细胞怎样从外表面进入内部。[Roh-Johnson 等，2012](https://doi.org/10.1126/science.1217869)；[Slabodnick 等，2023](https://doi.org/10.1371/journal.pgen.1010319)

细胞球研究保留了整个细胞群的位置、分裂和边界，能描述局部动作作用于怎样的堆积环境；但它将轮廓压缩为中心与大小，不能直接给出顶面怎样缩小。可变形表面、相场和整胚模型继续加入界面、体积与邻居形变，处理的是另一层空间自由度。Thiels等将整胚定量成像与早期原肠胚模型结合的工作目前仍为预印本，应与已发表的局部耦合实验分别看待。[Fickentscher 等，2013](https://doi.org/10.1016/j.bpj.2013.09.005)；[Thiels 等，2026，预印本](https://doi.org/10.64898/2026.03.30.715391)

从这些文献得到的物理图像是，内移的难易取决于主动结构、接触传力和周围组织的共同响应。只画一条内移坐标的势能曲线，可以解释稳定性，却没有同时解释细胞形状与邻居位移。[第一篇](/notes/celegans/embryonic_morphogenesis/gastrulation/)按模型分辨率详细比较了这些研究。

### 8.2 其他体系的借鉴 突起与接触更新怎样形成组织插入 {#方向偏置的主动粒子模型}

线虫背侧插入的结果，是具有左右身份的细胞从两列逐步形成交错排列。Walck-Shannon等在同一过程中观察细胞轮廓、基底外侧突起和方向性，并通过遗传干预区分突起形成与定向。[Walck-Shannon 等，2015](https://doi.org/10.1242/dev.127597)

![背侧插入的连续影像基底突起与方向分布](../figures/intercalation2015_fig1.jpg)

*图 20｜线虫背侧插入的空间过程。引自 Walck-Shannon 等（2015），[Fig. 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC4631769/#DEV127597F1)，完整原图。DIC序列、嵌合F-actin标记和重切片分别显示整体排列、单侧突起及其高度；它们共同界定需要解释的形态变化。*

与此最有用的模型参照，并不是单个细胞方向角的稳定性，而是保留轮廓、相邻细胞与动态接触的插入模型。Brodland的有限元研究较早把突起产力放到可变形细胞群和外部约束中，讨论组织伸长与细胞拉长怎样随负载改变。[Brodland，2006](https://doi.org/10.1387/ijdb.052040gb) 后续独立皮层与主动连接模型，进一步将局部活动接到空间重排。

#### 8.2.1 独立皮层与黏附使局部活动转成共同运动 {#published-crawling-cortex}

Weng等研究爪蟾脊索时，没有让相邻细胞永远共用一条不可分离的边。每个细胞用高分辨率多边形表示，两侧皮层可以分别运动，再由动态黏附相连。模型因此同时保留突起位置、皮层收缩、轮廓变形以及细胞之间的相对运动。[Weng 等，2022](https://doi.org/10.1016/j.celrep.2022.110666)

![爬行收缩与动态黏附模型的细胞和组织时序模拟](../figures/weng2022_fig5.jpg)

*图 21｜从四细胞插入到多细胞组织伸长。引自 Weng 等（2022），[Fig. 5](https://pmc.ncbi.nlm.nih.gov/articles/PMC9119128/#F5)，PMC作者稿原图，CC BY-NC-ND 4.0。A–E给出独立轮廓、黏附与活动，F、G为四细胞和27细胞的形态演化，H–L比较活动配合时的步进。*

模型的每个细胞用90个顶点保留局部轮廓，黏附连接相邻皮层，爬行与脉冲收缩再施加在相应位置。这样，突起能否带动邻居取决于接触传力，而不是在规定方向以后直接赋予一次成功的T1。

这套空间表示产生了比“方向正确就能插入”更具体的结果：爬行与收缩同时作用时，局部插入更有效；黏附过低时，四细胞仍能换邻居，较大细胞群却出现空隙，整体伸长减弱。局部重排与完整组织变形在这里被真正分开了。

![黏附降低后局部插入与组织伸长的不同结果](../figures/weng2022_fig6.jpg)

*图 22｜邻居交换完成，组织仍可能伸长不足。引自 Weng 等（2022），[Fig. 6](https://pmc.ncbi.nlm.nih.gov/articles/PMC9119128/#F6)，CC BY-NC-ND 4.0，完整原图。A、B比较低黏附下的四细胞和27细胞；空隙及组织轮廓直接显示局部活动未能有效连接成整体形变。*

论文还用Arvcf缺失后的步进变化检验爬行与收缩的协作。这为线虫提供的参照，是突起需要通过接触带动主体与邻居；爪蟾的活动布局不能直接等同于线虫的基底突起，也没有解释线虫顶底两层接触的时差。

#### 8.2.2 组织负载与主动连接共同选择T1 {#published-active-t1}

Sknepnek等的主动顶点模型沿另一条路线研究方向选择：外部应力造成力学各向异性，张力与myosin驻留之间的反馈再使特定连接收缩。细胞形状、连接张力和myosin在同一组织中演化，T1按模型的接触更新规则完成。[Sknepnek 等，2023](https://doi.org/10.7554/eLife.79862)

![主动T1发生时的细胞形态张力与myosin时序](../figures/sknepnek2023_fig2.jpg)

*图 23｜局部重排与周围承载组织同时保留。引自 Sknepnek 等（2023），[Fig. 2](https://elifesciences.org/articles/79862#fig2)，CC BY 4.0。模型包含主动、缓冲与被动细胞，连续帧记录T1、连接张力及myosin；外部拉力建立初始各向异性。*

不同活动与拉力组合产生的组织构型并不相同，收敛延伸只出现在相应参数范围。外部负载过弱或过强、活动不足或过强，都不能被归入同一条“活动越多，组织越容易正确重排”的趋势。

![活动与外部拉力共同影响多细胞形态的模型结果](../figures/sknepnek2023_fig6.jpg)

*图 24｜比较的是组织构型，不只是一个松弛时间。引自 Sknepnek 等（2023），[Fig. 6](https://elifesciences.org/articles/79862#fig6)，CC BY 4.0。各格为不同外部拉力与活动组合下的模拟组织；收敛延伸位于相应中间范围。*

Tetley等在果蝇中的研究则将细胞身份与条带边界纳入连接活动，说明定向重排还取决于组织的空间图案。[Tetley 等，2016](https://doi.org/10.7554/eLife.12094) 这些工作解释了接触重排怎样与承载共存；它们没有证明线虫背侧插入由相同的顶端收缩机制驱动。对线虫，基底活动先行而顶端接触后续更新，仍是区分模型层次的重要实验事实。

<figure class="celegans-video" id="active-review-published-intercalation">
<video controls loop preload="metadata" playsinline width="500" height="400" aria-label="已发表果蝇顶点模型的条带重排模拟">
<source src="https://static-movie-usa.glencoesoftware.com/mp4/10.7554/119/2afd92d2b201d28e34225b00378a697faa1afcc2/elife-12094-media5.mp4" type="video/mp4">
浏览器暂不支持内嵌播放，请打开下方原始视频。
</video>
<figcaption><strong>模型动画｜细胞身份与界面作用共同组织重排</strong><br>
Tetley 等（2016）Video 5，展示 Fig. 7G 的 simulation 4。颜色区分模型中的细胞身份；观察边界怎样移动，以及细胞交换邻居后是否仍维持条带次序。动画为果蝇研究中的计算结果。<br>
<span class="celegans-video-credit"><a href="https://elifesciences.org/articles/12094#video5">原文与视频说明</a> · <a href="https://doi.org/10.7554/eLife.12094.030">原始视频</a> · CC BY 4.0。</span></figcaption>
</figure>

### 8.3 其他体系的借鉴 新接触的几何怎样改变力传播 {#移动边界与新连接成熟}

Weng等的后续研究把问题推进到四细胞顶点解析之后：新接触不仅需要长出来，还需要在相邻拉力作用下调整方向，才能使局部插入与更远处的组织运动连接起来。作者结合爪蟾实验、静态力平衡和动态顶点模型，比较不同接触几何中的力传播。[Weng 等，2025](https://doi.org/10.1016/j.cub.2024.10.070)

在其模型中，近似六角排列与砖墙式排列对同一局部刺激给出不同的传播路径；动态模型中，使一处新接触开始转动，可以使原先停滞区域中的重排向外展开。Arvcf缺失后的接触取向与组织变形资料，为这条联系提供了实验对照。论文的[力分布图与时序模拟](https://pmc.ncbi.nlm.nih.gov/articles/PMC11706704/#F4)展示了这里的核心空间结果。

这项工作的启发是，T1的次数不是组织伸长的充分描述。事件留下的接触角度和排列，会改变后续的负载路径。它与线虫的问题相通之处，是局部重排以后仍须维持有效传力；具体的接触旋转机制及其在背侧插入中的作用，仍应与线虫实验区分。

### 8.4 线虫研究 骨架更新、连接变化与整体伸长 {#主动固体与固定邻接下的形状变化}

后期伸长中的关键联系，是肌肉产生周期性负载，表皮既要承受每轮形变，又要逐步保留身体长度。有关线虫的文献已经分别从网络重塑、连接动力学与三维身体结构讨论这件事，但各自保留的空间尺度不同。

#### 8.4.1 回缩表型与骨架参考状态更新 {#published-worm-actin-plasticity}

<span id="周期负载与黏塑性参考状态"></span>

Lardennois等将PAK-1与SPC-1相关网络异常、actin结构改变和全身回缩联系起来。这个表型说明，肌肉能够继续施加负载，表皮却未必还能维持此前形成的身体长度。作者的黏塑性模型与actin切断、重新束化的分子分析，共同解释形变如何得到稳定。[Lardennois 等，2019](https://doi.org/10.1038/s41586-019-1509-4)

![骨架相关机制受损后线虫伸长与回缩的实验轨迹](../figures/lock2019_fig2.jpg)

*图 25｜形变发生与形态保持可以分别受损。引自 Lardennois 等（2019），[公开作者稿 Fig. 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8086821/#F2)，完整原图，原权利保留。D、E比较正常、联合缺失及肌肉受损条件下的体长；作者稿图号与期刊版编排不同。*

Fang等进一步将seam的主动收缩、背腹表皮的方向性纤维，以及肌肉相关的切断与重新束化放入同一圆柱力学模型。重新束化既改变能承载的网络比例，也更新周向参考长度，再通过体积和力平衡转成轴向伸长。[Fang 等，2021](https://doi.org/10.1126/sciadv.abg3264)

![线虫伸长模型的侧表皮背腹纤维与骨架更新分工](../figures/fang2021_fig1.jpg)

*图 26｜局部材料过程进入整体伸长模型。引自 Fang 等（2021），[Fig. 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC8245039/#F1)，CC BY 4.0。模型区分两个伸长阶段、seam产力、背腹纤维与后期重塑；薄壁圆柱与区域均匀化是所用近似。*

它解释了为什么更新速率与承载能力需要一起看：拆解增加可以促进重塑，却也减少完整纤维；重新束化不足时，组织可能无法保留形态。论文用伸长、回缩及直径资料比较模型结果。不过，这类分区圆柱模型把各区的行为平均化，没有逐个预测每根actin在哪个位置断裂，也未给出完整表皮上的局部重排图。

#### 8.4.2 周期负载中的连接变化 {#published-worm-junction-turnover}

<span id="周期材料响应与实验比较"></span>

Yang等用光片成像追踪肌肉活动与表皮连接，观察到连接在折叠、较松与拉直、较紧的状态间反复变化；FRAP与遗传条件进一步将这种机械输入同E-cadherin周转联系起来。[Yang 等，2025](https://doi.org/10.1016/j.devcel.2025.06.012)

![线虫肌肉输入连接形变与持续伸长的同时观测](../figures/yang2025_fig2.jpg)

*图 27｜快速的连接形变与缓慢的长度增加并存。引自 Yang 等（2025），[Fig. 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC12549353/#fig2)，CC BY-NC 4.0。光片影像与局部长度轨迹记录周期变化，其他面板比较方向及肌肉受损条件。*

其模型把表皮细胞作为三维主动材料，联立三个方向的形变、连接浓度与线张力。连接伸长会稀释分子，运输和交换又改变浓度及力学作用；因此不是把一条FRAP曲线直接当成整个组织的松弛时间。原文用方向、位置与肌肉异常条件共同约束这个解释，连接运输受影响也不总导致伸长完全停止。

这项研究说明连接动态参与肌肉输入的效果，但并未把actin重塑和全部表皮接触纳入一个逐细胞空间模型。将它与Lardennois、Fang的结果并读，可以看到承载网络和连接动力学各自提供了哪一部分解释，而不必把它们改写成同一种“参考长度更新”。

#### 8.4.3 三维组织方位决定弯曲与扭转 {#published-worm-three-dimensional-shape}

Dai与Ben Amar把表皮、外周纤维、内部组织和四个肌肉区域放入有截面方向的三维身体。肌肉位于不同侧面，活动相位不同时，合力矩不会互相抵消；纤维倾角又使收缩与扭转耦合。[Dai 与 Ben Amar，2024，正式版](https://doi.org/10.7554/eLife.90505.3)

![肌肉方位和纤维方向对胚胎三维弯扭构型的影响](../figures/dai2023_fig7.jpg)

*图 28｜空间分区带来全长曲线之外的结果。引自 Dai 与 Ben Amar（2024），[Fig. 7](https://elifesciences.org/articles/90505#fig7)，CC BY 4.0。A–C、E为不同肌肉和actin条件下的模型构型，D为论文引用的胚胎影像，F、G比较曲率与扭转。*

这类模型可以讨论力作用在哪一侧、纤维朝哪个方向，以及局部作用怎样累积成全身弯扭。论文还通过能量分配规则将收缩与长期伸长相连；该规则包含对未显式处理过程的近似，不能仅凭形态相似就认为局部重塑路径已被唯一确定。

因此，现有线虫理论已经把“有形变”推进到“由哪些材料和几何产生形变”，但完整连接局部网络更新、表皮接触与每轮三维运动，仍不是任一单篇模型已经全部完成的工作。这里的不同研究是互补的，适合按其实际空间分辨率比较。

## 9. 文献综合：定向重排怎样与承载共存 {#10-怎样比较这些模型与线虫实验}

### 9.1 插入效率取决于活动怎样传到邻居 {#101-从现象到参数再到可区分的结果}

Weng等的独立皮层模型中，黏附降低使较大的细胞群出现空隙，局部交换与整体伸长因此分离；后续接触几何研究又显示，四细胞顶点完成解析以后留下的排列会影响力向远处传播。两类工作共同说明，活动的大小、细胞是否换邻居和组织是否有效伸长，不能相互替代。

线虫的基底外侧突起与顶端接触时序，使这项联系更具体：前缘、主体与邻接不必同步移动。其他体系的模型给出了传力和几何如何影响插入的实例，线虫实验则限定了这些过程实际发生在什么位置与细胞层中。

### 9.2 参考状态更新需要与承载结构一起理解 {#102-怎样从影片中分开整体运动与重排}

Tah等的闭合研究中，参考周长更新与持续刚性并存；线虫actin重塑文献中，重塑既关系到新的偏好形状，也关系到还有多少结构能够承载。它们从不同体系说明，形态变化可以依靠局部材料更新，并不要求所有连接同时失去回复能力。

这种联系也解释了为什么一条整体长度曲线不够。Lardennois的回缩表型区分形变发生与保留，Yang的局部连接资料区分快速周期变化与长期伸长，Dai与Ben Amar的三维构型则保留了活动方位。文献之间能相互补充的，正是这些不同的空间与时间观测。

### 9.3 各类已发表模型已经连接到哪一步 {#103-适合线虫的模型应保留哪些异质性}

| 文献与研究体系 | 保留的空间结构 | 论文建立的联系 | 尚未覆盖的层次 |
| --- | --- | --- | --- |
| Weng 2022，爪蟾 | 独立细胞轮廓、动态黏附、局部爬行与收缩 | 细胞插入、组织空隙与收敛延伸 | 线虫顶底极性及左右身份序列 |
| Sknepnek 2023，一般主动组织模型 | 活动区与周围细胞、张力及myosin | 外部负载与主动反馈选择重排 | 线虫基底突起的具体传力结构 |
| Weng 2025，爪蟾 | 接触取向、排列与组织力分布 | 顶点解析、力传播和组织形变传播 | 线虫表皮中的对应机制 |
| Fang 2021，线虫 | seam与背腹表皮的分区圆柱 | 方向性材料、重塑与伸长或回缩 | 逐位置的actin更新与表皮接触演化 |
| Yang 2025，线虫 | 三维细胞与方向不同的连接 | 肌肉输入、分子交换与连接伸长 | 完整多细胞表皮的邻接与负载网络 |
| Dai与Ben Amar 2024，线虫 | 分层截面、四组肌肉与三维中心线 | 活动方位、纤维方向与全身弯扭 | 局部重塑的分子过程及能量分配 |

这些工作支持一个共同的认识：形态发生需要某些接触沿特定方向重排，同时其他路径继续传力；也需要部分形变能够释放，另一些变化被结构更新保存。受限几何、极性、主动应力和材料重塑分别参与其中，具体分工随发育阶段改变。

## 10. 扩展阅读：主动向列与力网络拓扑 {#topology-nematics-force-networks}

<span id="11-参考文献"></span>

### 10.1 向列取向适合解释什么 {#6-取向场与主动向列液晶}

<span id="61-从细胞长轴或纤维方向得到-q-张量"></span>
<span id="62-活性怎样通过应力进入流动"></span>
<span id="63-正负半整数缺陷怎样定义"></span>
<span id="64-闭合曲面的总绕数与细胞边数是两种约束"></span>
<span id="65-缺陷附近的流动并不只由正负号决定"></span>

向列理论描述没有头尾之分的方向场，可用于比较纤维或细胞长轴的局部排列；缺陷则记录方向绕行时不能平滑消去的结构。培养细胞中的缺陷与局部流动、积聚或挤出已有实验联系，但这类证据应与线虫的谱系、曲率和定向骨架分别讨论。[Saw 等，2017](https://doi.org/10.1038/nature21718)；[Kawaguchi 等，2017](https://doi.org/10.1038/nature22321)

![取向图案中培养细胞的缺陷与集体运动](../figures/defects2025_fig1.png)

*图 29｜用明确的取向图案检验主动材料理论。引自 Zhao 等（2025），[Fig. 1](https://www.nature.com/articles/s41467-025-57783-w#Fig1)，CC BY-NC-ND 4.0，完整原图。体系为培养神经前体细胞，不是线虫胚胎。*

<figure class="celegans-video celegans-stage-video" id="active-defect-motion" style="--stage-video-width: 34rem;">
<video controls loop preload="metadata" playsinline width="560" height="560" aria-label="取向缺陷图案中培养细胞的集体运动" aria-describedby="active-defect-motion-caption">
<source src="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-57783-w/MediaObjects/41467_2025_57783_MOESM4_ESM.mp4" type="video/mp4">
浏览器暂不支持内嵌播放，请打开下方原始视频。
</video>
<figcaption id="active-defect-motion-caption"><strong>动态影像｜取向图案中的集体运动</strong><br>
培养神经前体细胞在诱导取向图案中的动态，可与上方论文图一起观察细胞运动与核心积聚。此处展示的是其他细胞体系；影片中的空间尺度与时间沿用原始标注。<br>
<span class="celegans-video-credit">来源：Zhao 等（2025），Supplementary Video 2，Nature Communications，CC BY-NC-ND 4.0。<a href="https://www.nature.com/articles/s41467-025-57783-w">原始页面与说明</a>。</span></figcaption>
</figure>

在线虫伸长中，周向actin排列可以直接进入方向性承载的描述，无须先把每处取向不均匀都解释成缺陷。手性皮层流又具有左右不对称的动力学含义，与无头尾方向场的绕数不是同一种量。取向图像、流动与具体受力结构相互对照时，这套理论才提供额外解释。

### 10.2 接触网络与力网络的拓扑 {#7-力网络与持续同调}

<span id="71-接触网络与力网络有不同的边权"></span>
<span id="72-图上的环路并不等于组织里的空洞"></span>
<span id="73-持续同调为什么比单一阈值多保留一些信息"></span>

接触图记录谁与谁相邻，力网络还需要为接触赋予有物理意义的权重。持续同调追踪不同权重阈值下连通分量与环路怎样出现，可以保留单一阈值丢掉的结构信息；它最初在颗粒力网络中得到系统发展。[Kramár 等，2013](https://doi.org/10.1103/PhysRevE.87.042207)

这种统计不自动确定传力机制。一个图上的三边环路不一定是组织中的空洞，荧光强度也不直接等于连接拉力。对这里讨论的线虫文献，细胞身份、接触方向、局部材料与实际负载已经解释了许多空间差别；拓扑工具适合作为进一步描述这些结构的方法，不宜取代前面的形态机制。

## 参考文献 {#references}

1. Corey S. O’Hern, Leonardo E. Silbert, Andrea J. Liu, Sidney R. Nagel (2003). [Jamming at zero temperature and zero applied stress: The epitome of disorder](https://doi.org/10.1103/PhysRevE.68.011306). *Physical Review E*.

2. James Giammona, Otger Campàs (2021). [Physical constraints on early blastomere packings](https://doi.org/10.1371/journal.pcbi.1007994). *PLOS Computational Biology*.

3. Silke Henkes, Yaouen Fily, M. Cristina Marchetti (2011). [Active jamming: Self-propelled soft particles at high density](https://doi.org/10.1103/PhysRevE.84.040301). *Physical Review E*.

4. J. Tailleur, M. E. Cates (2008). [Statistical Mechanics of Interacting Run-and-Tumble Bacteria](https://doi.org/10.1103/PhysRevLett.100.218103). *Physical Review Letters*.

5. Yaouen Fily, M. Cristina Marchetti (2012). [Athermal Phase Separation of Self-Propelled Particles with No Alignment](https://doi.org/10.1103/PhysRevLett.108.235702). *Physical Review Letters*.

6. Dapeng Bi, J. H. Lopez, J. M. Schwarz, M. Lisa Manning (2015). [A density-independent rigidity transition in biological tissues](https://doi.org/10.1038/nphys3471). *Nature Physics*.

7. Dapeng Bi, Xingbo Yang, M. Cristina Marchetti, M. Lisa Manning (2016). [Motility-Driven Glass and Jamming Transitions in Biological Tissues](https://doi.org/10.1103/PhysRevX.6.021011). *Physical Review X*.

8. Daniel M. Sussman, Matthias Merkel (2018). [No unjamming transition in a Voronoi model of biological tissue](https://doi.org/10.1039/c7sm02127e). *Soft Matter*.

9. Matthias Merkel, M Lisa Manning (2018). [A geometrically controlled rigidity transition in a model for confluent 3D tissues](https://doi.org/10.1088/1367-2630/aaaa13). *New Journal of Physics*.

10. Pedro Gómez-Gálvez, Pablo Vicente-Munuera, Antonio Tagua, Cristina Forja, Ana M. Castro, Marta Letrán, Andrea Valencia-Expósito, Clara Grima, Marina Bermúdez-Gallardo, Óscar Serrano-Pérez-Higueras, Florencia Cavodeassi, Sol Sotillos, María D. Martín-Bermudo, Alberto Márquez, Javier Buceta, Luis M. Escudero (2018). [Scutoids are a geometrical solution to three-dimensional packing of epithelia](https://doi.org/10.1038/s41467-018-05376-1). *Nature Communications*.

11. Felix C. Keber, Etienne Loiseau, Tim Sanchez, Stephen J. DeCamp, Luca Giomi, Mark J. Bowick, M. Cristina Marchetti, Zvonimir Dogic, Andreas R. Bausch (2014). [Topology and dynamics of active nematic vesicles](https://doi.org/10.1126/science.1254784). *Science*.

12. Thuan Beng Saw, Amin Doostmohammadi, Vincent Nier, Leyla Kocgozlu, Sumesh Thampi, Yusuke Toyama, Philippe Marcq, Chwee Teck Lim, Julia M. Yeomans, Benoit Ladoux (2017). [Topological defects in epithelia govern cell death and extrusion](https://doi.org/10.1038/nature21718). *Nature*.

13. Kyogo Kawaguchi, Ryoichiro Kageyama, Masaki Sano (2017). [Topological defects control collective dynamics in neural progenitor cell cultures](https://doi.org/10.1038/nature22321). *Nature*.

14. Zihui Zhao, He Li, Yisong Yao, Yongfeng Zhao, Francesca Serra, Kyogo Kawaguchi, Hepeng Zhang, Masaki Sano (2025). [Integer topological defects offer a methodology to quantify and classify active cell monolayers](https://doi.org/10.1038/s41467-025-57783-w). *Nature Communications*.

15. M. Kramar, A. Goullet, L. Kondic, K. Mischaikow (2013). [Persistence of force networks in compressed granular media](https://doi.org/10.1103/PhysRevE.87.042207). *Physical Review E*.

16. Herbert Edelsbrunner, David Letscher, Afra Zomorodian (2002). [Topological Persistence and Simplification](https://doi.org/10.1007/s00454-002-2885-2). *Discrete & Computational Geometry*.

17. Michael L. Falk, James S. Langer (1998). [Dynamics of viscoplastic deformation in amorphous solids](https://doi.org/10.1103/PhysRevE.57.7192). *Physical Review E*.

18. Fickentscher R, Struntz P, Weiss M (2013). [Mechanical Cues in the Early Embryogenesis of Caenorhabditis elegans](https://doi.org/10.1016/j.bpj.2013.09.005). *Biophys J*.

19. Kazunori Yamamoto, Akatsuki Kimura (2017). [An asymmetric attraction model for the diversity and robustness of cell arrangement in nematodes](https://doi.org/10.1242/dev.154609). *Development*.

20. Kazunori Yamamoto, Sacha Ichbiah, Matthieu Perez, Joana Borrego-Pinto, Fabrice Delbary, Nathan W. Goehring, Paul Northrop, Hervé Turlier, Guillaume Charras (2026). [Spatiotemporal mapping of the contractile and adhesive forces sculpting early C. elegans embryos](https://doi.org/10.1016/j.devcel.2026.04.013). *Developmental Cell*.

21. Roh-Johnson M, Shemer G, Higgins CD, McClellan JH, Werts AD, Tulu US, Gao L, Betzig E, Kiehart DP, Goldstein B (2012). [Triggering a cell shape change by exploiting preexisting actomyosin contractions](https://doi.org/10.1126/science.1217869). *Science*.

22. Slabodnick MM, Tintori SC, Prakash M, Zhang P, Higgins CD, Chen AH, Cupp TD, Wong T, Bowie E, Jug F, Goldstein B (2023). [Zyxin contributes to coupling between cell junctions and contractile actomyosin networks during apical constriction](https://doi.org/10.1371/journal.pgen.1010319). *PLoS Genet*.

23. Wim Thiels, Michiel Vanslambrouck, Casper Van Bavel, Keyu Xiao, Jef Vangheel, Bart Smeets, Rob Jelier (2026). [Integrated quantitative imaging and biomechanical modeling of early gastrulation in C. elegans](https://doi.org/10.64898/2026.03.30.715391). *bioRxiv preprint*.

24. Chihara D, Nance J (2012). [An E-cadherin-mediated hitchhiking mechanism for C. elegans germ cell internalization during gastrulation](https://doi.org/10.1242/dev.079863). *Development*.

25. Lokesh G. Pimpale, Teije C. Middelkoop, Alexander Mietke, Stephan W. Grill (2020). [Cell lineage-dependent chiral actomyosin flows drive cellular rearrangements in early Caenorhabditis elegans development](https://doi.org/10.7554/eLife.54930). *eLife*.

26. Mi Jing Khor, Yuxuan Rain Xiong, Gaganpreet Sangha, Kenji Sugioka (2025). [Cadherin flow translates cortical chirality into left-right asymmetric cell division dynamics in Caenorhabditis elegans embryos](https://doi.org/10.1016/j.cub.2025.09.052). *Current Biology*.

27. Walck-Shannon E, Reiner D, Hardin J (2015). [Polarized Rac-dependent protrusions drive epithelial intercalation in the embryonic epidermis of C. elegans](https://doi.org/10.1242/dev.127597). *Development*.

28. Walck-Shannon E, Lucas B, Chin-Sang I, Reiner D, Kumfer K, Cochran H, Bothfeld W, Hardin J (2016). [CDC-42 Orients Cell Migration during Epithelial Intercalation in the Caenorhabditis elegans Epidermis](https://doi.org/10.1371/journal.pgen.1006415). *PLoS Genet*.

29. Zhu Y, Tesone Z, Tan M, Hardin J (2024). [TIAM-1 regulates polarized protrusions during dorsal intercalation in the Caenorhabditis elegans embryo through both its GEF and N-terminal domains](https://doi.org/10.1242/jcs.261509). *J Cell Sci*.

30. Vuong-Brender TT, Ben Amar M, Pontabry J, Labouesse M (2017). [The interplay of stiffness and force anisotropies drives embryo elongation](https://doi.org/10.7554/eLife.23866). *eLife*.

31. Fang C, Wei X, Shao X, Lin Y (2021). [Force-mediated cellular anisotropy and plasticity dictate the elongation dynamics of embryos](https://doi.org/10.1126/sciadv.abg3264). *Sci Adv*.

32. Lardennois A, Pásti G, Ferraro T, Llense F, Mahou P, Pontabry J, Rodriguez D, Kim S, Ono S, Beaurepaire E, Gally C, Labouesse M (2019). [An actin-based viscoplastic lock ensures progressive body-axis elongation](https://doi.org/10.1038/s41586-019-1509-4). *Nature*.

33. Yang X, Ferraro T, Molnar K, Pontabry J, Malanda SR, Maghelli N, Royer L, Grill SW, Myers G, Grigolon S, Labouesse M (2025). [Repeated extrinsic and anisotropic mechanical inputs promote C. elegans polarized adherens junction elongation](https://doi.org/10.1016/j.devcel.2025.06.012). *Dev Cell*.

34. Hiroshi Koyama, Hisashi Okumura, Atsushi M. Ito, Kazuyuki Nakamura, Tetsuhisa Otani, Kagayaki Kato, Toshihiko Fujimori (2023). [Effective mechanical potential of cell–cell interaction explains three-dimensional morphologies during early embryogenesis](https://doi.org/10.1371/journal.pcbi.1011306). *PLOS Computational Biology*.

35. Miao J, Guan G, Tang C (2021; v2 2022). [Spontaneous mechanical and energetic state transitions during Caenorhabditis elegans gastrulation](https://doi.org/10.48550/arXiv.2105.05795). *arXiv preprint*.

36. William A. Mohler, Gidi Shemer, Jacob J. del Campo, Clari Valansi, Eugene Opoku-Serebuoh, Victoria Scranton, Nirit Assaf, John G. White, Benjamin Podbilewicz (2002). [The type I membrane protein EFF-1 is essential for developmental cell fusion](<https://doi.org/10.1016/S1534-5807(02)00129-6>). *Developmental Cell*.

37. Williams-Masson EM, Malik AN, Hardin J (1997). [An actin-mediated two-step mechanism is required for ventral enclosure of the C. elegans hypodermis](<https://doi.org/10.1242/dev.124.15.2889>). *Development*.

38. Raich WB, Agbunag C, Hardin J (1999). [Rapid epithelial-sheet sealing in the Caenorhabditis elegans embryo requires cadherin-dependent filopodial priming](<https://doi.org/10.1016/S0960-9822(00)80015-9>). *Curr Biol*.

39. Bernadskaya YY, Wallace A, Nguyen J, Mohler WA, Soto MC (2012). [UNC-40/DCC, SAX-3/Robo, and VAB-1/Eph polarize F-actin during embryonic morphogenesis by regulating the WAVE/SCAR actin nucleation complex](<https://doi.org/10.1371/journal.pgen.1002863>). *PLoS Genet*.

40. Alessandro Mongera, Payam Rowghanian, Hannah J. Gustafson, Elijah Shelton, David A. Kealhofer, Emmet K. Carn, Friedhelm Serwane, Adam A. Lucio, James Giammona, Otger Campàs (2018). [A fluid-to-solid jamming transition underlies vertebrate body axis elongation](https://doi.org/10.1038/s41586-018-0479-2). *Nature*.

41. Nicoletta I. Petridou, Bernat Corominas-Murtra, Carl-Philipp Heisenberg, Edouard Hannezo (2021). [Rigidity percolation uncovers a structural basis for embryonic tissue phase transitions](https://doi.org/10.1016/j.cell.2021.02.017). *Cell*.

42. Suzanne M. Fielding, James O. Cochran, Junxiang Huang, Dapeng Bi, M. Cristina Marchetti (2023). [Constitutive model for the rheology of biological tissue](https://doi.org/10.1103/PhysRevE.108.L042602). *Physical Review E*.

43. Saaya Koizumi, Ayama Tokuyasu, Akinori W. M. Miyamoto, Takayuki Torisawa, Hirokazu Tanimoto, Akatsuki Kimura (2026). [Developmentally programmed changes in cytoplasmic mechanics revealed by active microrheology in C. elegans embryos](https://doi.org/10.64898/2026.05.19.726147). *bioRxiv preprint*.

44. Tah I, Haertter D, Crawford JM, Kiehart DP, Schmidt CF, Liu AJ (2025). [A minimal vertex model explains how the amnioserosa avoids fluidization during Drosophila dorsal closure](https://doi.org/10.1073/pnas.2322732121). *PNAS*.


45. Goldstein B (1992). [Induction of gut in Caenorhabditis elegans embryos](https://doi.org/10.1038/357255a0). *Nature*.

46. Goldstein B (1995). [Cell contacts orient some cell division axes in the Caenorhabditis elegans embryo](https://doi.org/10.1083/jcb.129.4.1071). *Journal of Cell Biology*.

47. Brodland GW (2006). [Do lamellipodia have the mechanical capacity to drive convergent extension?](https://doi.org/10.1387/ijdb.052040gb). *International Journal of Developmental Biology*.

48. Weng S, Huebner RJ, Wallingford JB (2022). [Convergent extension requires adhesion-dependent biomechanical integration of cell crawling and junction contraction](https://doi.org/10.1016/j.celrep.2022.110666). *Cell Reports*.

49. Sknepnek R, Djafer-Cherif I, Chuai M, Weijer C, Henkes S (2023). [Generating active T1 transitions through mechanochemical feedback](https://doi.org/10.7554/eLife.79862). *eLife*.

50. Tetley RJ, Blanchard GB, Fletcher AG, Adams RJ, Sanson B (2016). [Unipolar distributions of junctional Myosin II identify cell stripe boundaries that drive cell intercalation throughout Drosophila axis extension](https://doi.org/10.7554/eLife.12094). *eLife*.

51. Weng S, Hayashi M, Inoue Y, Wallingford JB (2025). [Planar polarized force propagation integrates cell behavior with tissue shaping during convergent extension](https://doi.org/10.1016/j.cub.2024.10.070). *Current Biology*.

52. Dai A, Ben Amar M (2024). [Cyclic muscle contractions reinforce the actomyosin motors and mediate the full elongation of C. elegans embryo](https://doi.org/10.7554/eLife.90505.3). *eLife*.

← [05 三维体形与模型辨识](/notes/celegans/embryonic_morphogenesis/geometry_and_inference/) · [系列总览](/notes/celegans/embryonic_morphogenesis/) · [07 亚细胞骨架与流动](/notes/celegans/embryonic_morphogenesis/subcellular_cortex_and_flows/) →
