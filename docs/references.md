# 参考文献（相位图校验、无衍射光束与旁瓣抑制、AI 辅助光场调控）

> 所有条目均通过 Crossref REST API 按 DOI 逐条核验（检索日期 2026-09-16），
> 标题、作者、期刊、卷期页码取自 Crossref 原始记录。
> 🔹 = 与本案例的方法或结论直接对应；无 DOI 的经典文献已单独标注。


## A. 轴锥镜与无衍射（贝塞尔）光束

- **[1]** McLeod (1954). The Axicon: A New Type of Optical Element. *Journal of the Optical Society of America*, **44**, 592. doi: [10.1364/josa.44.000592](https://doi.org/10.1364/josa.44.000592) 🔹  
  轴锥镜（axicon）原始文献，定义锥面波前
- **[2]** Durnin et al. (1987). Diffraction-free beams. *Physical Review Letters*, **58**, 1499-1501. doi: [10.1103/physrevlett.58.1499](https://doi.org/10.1103/physrevlett.58.1499) 🔹  
  首次实验实现无衍射光束
- **[3]** Durnin (1987). Exact solutions for nondiffracting beams I The scalar theory. *Journal of the Optical Society of America A*, **4**, 651. doi: [10.1364/josaa.4.000651](https://doi.org/10.1364/josaa.4.000651) 🔹  
  无衍射光束的标量理论解（J₀ 型 Bessel 光束）
- **[4]** Scott (1992). Efficient generation of nearly diffraction-free beams using an axicon. *Optical Engineering*, **31**, 2640. doi: [10.1117/12.60017](https://doi.org/10.1117/12.60017) 🔹  
  轴锥镜产生准无衍射光束的效率与实用形式
- **[5]** McGloin & Dholakia (2005). Bessel beams: Diffraction in a new light. *Contemporary Physics*, **46**, 15-28. doi: [10.1080/0010751042000275259](https://doi.org/10.1080/0010751042000275259) 🔹  
  Bessel 光束综述：无衍射区长度、自重建、旁瓣结构
- **[6]** Gerchberg, R. W., & Saxton, W. O. (1972). A practical algorithm for the determination of phase from image and diffraction plane pictures. *Optik*, **35**, 237–246. （无 DOI，DOI 时代之前的经典文献） —— 迭代相位恢复（GS 算法）原始文献，相位型 DOE 设计的起点 🔹

## B. 贝塞尔光束的整形、切趾与旁瓣抑制

- **[7]** Gori et al. (1987). Bessel-Gauss beams. *Optics Communications*, **64**, 491-495. doi: [10.1016/0030-4018(87)90276-8](https://doi.org/10.1016/0030-4018(87)90276-8) 🔹  
  Bessel-Gauss 光束：高斯切趾对环带包络的影响
- **[8]** Jiang (1996). Super-Gaussian-Bessel beam. *Optics Communications*, **125**, 207-210. doi: [10.1016/0030-4018(95)00740-7](https://doi.org/10.1016/0030-4018(95)00740-7) 🔹  
  超高斯切趾的 Bessel 光束，直接用于压低声旁瓣
- **[9]** Mori (2015). Side lobe suppression of a Bessel beam for high aspect ratio laser processing. *Precision Engineering*, **39**, 79-85. doi: [10.1016/j.precisioneng.2014.07.008](https://doi.org/10.1016/j.precisioneng.2014.07.008) 🔹  
  面向激光加工的 Bessel 光束旁瓣抑制（与本项目目标最接近）
- **[10]** Bélanger & Rioux (1978). Ring pattern of a lens–axicon doublet illuminated by a Gaussian beam. *Applied Optics*, **17**, 1080. doi: [10.1364/ao.17.001080](https://doi.org/10.1364/ao.17.001080) 🔹  
  透镜-轴锥复合系统的环状焦斑分析（几何光学极限）
- **[11]** Sheppard & Wilson (1978). Gaussian-beam theory of lenses with annular aperture. *IEE Journal on Microwaves, Optics and Acoustics*, **2**, 105. doi: [10.1049/ij-moa.1978.0023](https://doi.org/10.1049/ij-moa.1978.0023)  
  环形口径透镜的高斯光束理论：旁瓣/环带与 NA 的定量关系
- **[12]** Rogers & Zheludev (2013). Optical super-oscillations: sub-wavelength light focusing and super-resolution imaging. *Journal of Optics*, **15**, 094008. doi: [10.1088/2040-8978/15/9/094008](https://doi.org/10.1088/2040-8978/15/9/094008)  
  超振荡聚焦：突破衍射极限与旁瓣代价的理论框架

## C. 螺旋相位与涡旋光束

- **[13]** Allen et al. (1992). Orbital angular momentum of light and the transformation of Laguerre-Gaussian laser modes. *Physical Review A*, **45**, 8185-8189. doi: [10.1103/physreva.45.8185](https://doi.org/10.1103/physreva.45.8185) 🔹  
  轨道角动量与 Laguerre-Gauss 模的经典文献
- **[14]** Bazhenov et al. (1992). Screw Dislocations in Light Wavefronts. *Journal of Modern Optics*, **39**, 985-990. doi: [10.1080/09500349214551011](https://doi.org/10.1080/09500349214551011) 🔹  
  螺旋位错（涡旋）波前的早期工作
- **[15]** Yao & Padgett (2011). Orbital angular momentum: origins, behavior and applications. *Advances in Optics and Photonics*, **3**, 161. doi: [10.1364/aop.3.000161](https://doi.org/10.1364/aop.3.000161)  
  轨道角动量综述：产生、测量与应用
- **[16]** Ostrovsky et al. (2013). Generation of the “perfect” optical vortex using a liquid-crystal spatial light modulator. *Optics Letters*, **38**, 534. doi: [10.1364/ol.38.000534](https://doi.org/10.1364/ol.38.000534) 🔹  
  完美涡旋光束（半径与拓扑荷解耦），本项目的环半径标定可对照

## D. SLM / CGH 产生无衍射光束

- **[17]** Turunen et al. (1988). Holographic generation of diffraction-free beams. *Applied Optics*, **27**, 3959. doi: [10.1364/ao.27.003959](https://doi.org/10.1364/ao.27.003959) 🔹  
  全息产生无衍射光束（计算全息 + 轴锥相位）
- **[18]** Vasara et al. (1989). Realization of general nondiffracting beams with computer-generated holograms. *Journal of the Optical Society of America A*, **6**, 1748. doi: [10.1364/josaa.6.001748](https://doi.org/10.1364/josaa.6.001748) 🔹  
  用 CGH 实现一般无衍射光束（相位图案的解析构造）
- **[19]** Chattrapiban et al. (2003). Generation of nondiffracting Bessel beams by use of a spatial light modulator. *Optics Letters*, **28**, 2183. doi: [10.1364/ol.28.002183](https://doi.org/10.1364/ol.28.002183) 🔹  
  用空间光调制器产生无衍射 Bessel 光束（与本项目「给定一张相位图」的实验场景对应）
- **[20]** Leach et al. (2006). Generation of achromatic Bessel beams using a compensated spatial light modulator. *Optics Express*, **14**, 5581. doi: [10.1364/oe.14.005581](https://doi.org/10.1364/oe.14.005581)  
  消色差 Bessel 光束与 SLM 补偿（相位图案标定的实验细节）
- **[21]** Xia (2009). Three-dimensional light modulation using phase-only spatial light modulator. *Optical Engineering*, **48**, 020502. doi: [10.1117/1.3076211](https://doi.org/10.1117/1.3076211)  
  纯相位 SLM 的三维光场调制（器件层面的编码与标定）

## E. 衍射传播数值方法

- **[22]** Matsushima & Shimobaba (2009). Band-limited angular spectrum method for numerical simulation of free-space propagation in far and near fields. *Optics Express*, **17**, 19662. doi: [10.1364/oe.17.019662](https://doi.org/10.1364/oe.17.019662) 🔹  
  带限角谱法：本项目传播算子的直接依据
- **[23]** Voelz (2011). Computational Fourier Optics: A MATLAB Tutorial. doi: [10.1117/3.858456](https://doi.org/10.1117/3.858456) 🔹  
  计算傅里叶光学教材：ASM/单 FFT/卷积三种方法及其适用条件
- **[24]** Goodman, J. W. (2017). *Introduction to Fourier Optics* (4th ed.). W. H. Freeman. （教材，无 DOI） —— 标量衍射与傅里叶光学基础 🔹

## F. 相位型 DOE 的优化与逆设计

- **[25]** Wyrowski & Bryngdahl (1988). Iterative Fourier-transform algorithm applied to computer holography. *Journal of the Optical Society of America A*, **5**, 1058. doi: [10.1364/josaa.5.001058](https://doi.org/10.1364/josaa.5.001058) 🔹  
  迭代傅里叶变换算法（IFTA）用于计算全息，GS 的实用化形式
- **[26]** Chakravarthula et al. (2019). Wirtinger holography for near-eye displays. *ACM Transactions on Graphics*, **38**, 1-13. doi: [10.1145/3355089.3356539](https://doi.org/10.1145/3355089.3356539) 🔹  
  Wirtinger 全息：以解析梯度优化纯相位全息图（与伴随法同源）
- **[27]** Molesky et al. (2018). Inverse design in nanophotonics. *Nature Photonics*, **12**, 659-670. doi: [10.1038/s41566-018-0246-9](https://doi.org/10.1038/s41566-018-0246-9) 🔹  
  纳米光子学逆设计综述：伴随法、目标函数与拓扑优化
- **[28]** Hughes et al. (2018). Adjoint Method and Inverse Design for Nonlinear Nanophotonic Devices. *ACS Photonics*, **5**, 4781-4787. doi: [10.1021/acsphotonics.8b01522](https://doi.org/10.1021/acsphotonics.8b01522) 🔹  
  伴随法逆设计的完整推导（∂J/∂φ 的解析梯度形式）
- **[29]** Swanson (1989). Binary Optics Technology: The Theory and Design of Multi-Level Diffractive Optical Elements. doi: [10.21236/ada213404](https://doi.org/10.21236/ada213404)  
  多台阶二元光学 DOE 的设计与效率（相位量化的工程背景）
- **[30]** Dammann & Görtler (1971). High-efficiency in-line multiple imaging by means of multiple phase holograms. *Optics Communications*, **3**, 312-315. doi: [10.1016/0030-4018(71)90095-2](https://doi.org/10.1016/0030-4018(71)90095-2)  
  Dammann 光栅：相位型分束与振幅编码思想

## G. 应用：光镊 / 超分辨 / 激光加工

- **[31]** Garcés-Chávez et al. (2002). Simultaneous micromanipulation in multiple planes using a self-reconstructing light beam. *Nature*, **419**, 145-147. doi: [10.1038/nature01007](https://doi.org/10.1038/nature01007) 🔹  
  自重建光束的多平面光操控（Bessel 光束长焦深优势）
- **[32]** Klar et al. (2000). Fluorescence microscopy with diffraction resolution barrier broken by stimulated emission. *Proceedings of the National Academy of Sciences*, **97**, 8206-8210. doi: [10.1073/pnas.97.15.8206](https://doi.org/10.1073/pnas.97.15.8206) 🔹  
  STED 超分辨中的环形/零中心光斑需求（本器件暗核结构的典型用途）
- **[33]** Duocastella & Arnold (2012). Bessel and annular beams for materials processing. *Laser &amp; Photonics Reviews*, **6**, 607-621. doi: [10.1002/lpor.201100031](https://doi.org/10.1002/lpor.201100031) 🔹  
  Bessel/环形光束用于材料加工：长焦深与旁瓣的工艺影响
- **[34]** Fahrbach et al. (2010). Microscopy with self-reconstructing beams. *Nature Photonics*, **4**, 780-785. doi: [10.1038/nphoton.2010.204](https://doi.org/10.1038/nphoton.2010.204)  
  自重建光束显微：无衍射区长度与成像深度

## H. AI / 可微优化 + 光场调控

- **[35]** Lin et al. (2018). All-optical machine learning using diffractive deep neural networks. *Science*, **361**, 1004-1008. doi: [10.1126/science.aat8084](https://doi.org/10.1126/science.aat8084) 🔹  
  衍射深度神经网络（D2NN）：把光传播直接纳入可训练模型
- **[36]** Sitzmann et al. (2018). End-to-end optimization of optics and image processing for achromatic extended depth of field and super-resolution imaging. *ACM Transactions on Graphics*, **37**, 1-13. doi: [10.1145/3197517.3201333](https://doi.org/10.1145/3197517.3201333) 🔹  
  端到端可微光学：光学元件与后端处理联合优化
- **[37]** Horisaki et al. (2018). Deep-learning-generated holography. *Applied Optics*, **57**, 3859. doi: [10.1364/ao.57.003859](https://doi.org/10.1364/ao.57.003859) 🔹  
  深度学习生成全息图（用网络替代迭代优化）
- **[38]** So et al. (2020). Deep learning enabled inverse design in nanophotonics. *Nanophotonics*, **9**, 1041-1057. doi: [10.1515/nanoph-2019-0474](https://doi.org/10.1515/nanoph-2019-0474) 🔹  
  深度学习驱动的纳米光子逆设计综述
- **[39]** Liu et al. (2018). Training Deep Neural Networks for the Inverse Design of Nanophotonic Structures. *ACS Photonics*, **5**, 1365-1369. doi: [10.1021/acsphotonics.7b01377](https://doi.org/10.1021/acsphotonics.7b01377)  
  训练神经网络做光子结构逆设计（数据集 + 正向代理模型）
- **[40]** Raissi et al. (2019). Physics-informed neural networks: A deep learning framework for solving forward and inverse problems involving nonlinear partial differential equations. *Journal of Computational Physics*, **378**, 686-707. doi: [10.1016/j.jcp.2018.10.045](https://doi.org/10.1016/j.jcp.2018.10.045)  
  物理信息神经网络（PINN）：把物理方程作为约束嵌入训练

## I. 结构化光场综述（背景与总览）

- **[41]** Forbes et al. (2021). Structured light. *Nature Photonics*, **15**, 253-262. doi: [10.1038/s41566-021-00780-4](https://doi.org/10.1038/s41566-021-00780-4) 🔹  
  结构化光场综述（模态、自由度与器件）
- **[42]** Rosales-Guzmán & Forbes (2024). Structured Light with Spatial Light Modulators. doi: [10.1117/3.100024](https://doi.org/10.1117/3.100024) 🔹  
  SLM 结构化光场专著：相位图案生成、标定与实验实现
- **[43]** Rubinsztein-Dunlop et al. (2016). Roadmap on structured light. *Journal of Optics*, **19**, 013001. doi: [10.1088/2040-8978/19/1/013001](https://doi.org/10.1088/2040-8978/19/1/013001) 🔹  
  结构化光场路线图：SLM/DOE/超表面的能力边界

---

## 检索与核验说明

```
核验方式：Crossref REST API（https://api.crossref.org/works/{doi}）
检索日期：2026-09-16
字段来源：title / author / container-title / issued / volume / page / DOI
```

若需批量导入文献管理软件，可用 DOI 直接查询 BibTeX：

```bash
curl -LH "Accept: application/x-bibtex" https://doi.org/10.1364/oe.17.019662
```
