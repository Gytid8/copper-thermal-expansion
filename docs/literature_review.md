# 单晶铜热膨胀文献比较与模拟路线

调研日期：2026年10月17日。

建议采用Python调用LAMMPS，在三维周期性边界下建立FCC单晶铜，用Mishin 2001 EAM1势进行各温度点的零外压NPT模拟。重点评价相对晶格膨胀和300–800 K温区平均线膨胀系数。已有工作支持该方法，课程项目的价值是可靠复现和误差分析。

## 文献范围与证据强度

共整理21篇文献，包含直接热膨胀研究、实验及参考数据、势函数研究、统计和软件方法。并非21篇都直接研究300–800 K单晶铜。2021年PACE、2023年晶粒尺寸研究、Narasimhan的铜热学原文及统计指南核对了相关全文段落；其他条目主要核对出版方摘要、作者公开版本或NIST原始势条目。未获取全文的工作不据此填写原子数、步长和时长。

检索使用copper thermal expansion molecular dynamics、copper thermal properties quasiharmonic、copper EAM potential comparison等关键词，并沿原文参考文献补充经典工作。本地paper-search脚本的数据库请求受代理限制（ProxyError / WinError 10013），且已安装版本不支持--json参数。因此以下为网页逐篇核对的定向文献集，不是数据库系统综述；没有可靠的引用次数统计，不据此给作者或论文排名。

## 直接研究铜热膨胀及有限温度性质

| 编号 | 文献与链接 | 方法和结果 | 对本课题的作用 |
|---|---|---|---|
| 1 | Alper & Politzer, 1999. Molecular dynamics simulation of the temperature-dependent behavior of solid copper. https://doi.org/10.1016/S0166-1280(99)00147-5 | 多体势、NPT、300–1100 K；计算体膨胀、压缩率与热容，总体趋势和量级较好 | 支持逐温度平衡后计算物性；注意体膨胀与线膨胀区别 |
| 2 | Kazanc et al., 2006. Temperature and pressure dependence of the some elastic and lattice dynamical properties of copper: a molecular dynamics study. https://doi.org/10.1016/j.physb.2005.12.259 | EAM框架中的修正Morse形式；同时研究膨胀、弹性和声子并对照实验 | 膨胀精度需专门检验，不由单项静态性质决定 |
| 3 | Narasimhan & de Gironcoli, 2001预印本. Ab initio calculation of the thermal properties of Cu: Performance of the LDA and GGA. https://arxiv.org/abs/cond-mat/0109020 | DFT、DFPT、准谐自由能；报告300 K膨胀系数偏差约LDA -14%、GGA +14% | 第一性原理也有近似误差；不能把单温度偏差当作全温区平均误差 |
| 4 | Rudin et al., 2002. First-principles-based thermodynamic description of solid copper using the tight-binding approach. https://doi.org/10.1103/PhysRevB.65.235114 ; https://arxiv.org/abs/cond-mat/0201141 | 第一性原理拟合紧束缚模型与准谐近似，计算膨胀及其他热学性质 | 展示自由能路线，与直接MD采样路线形成对照 |
| 5 | Lysogorskiy et al., 2021. Performant implementation of the atomic cluster expansion (PACE) and application to copper and silicon. https://doi.org/10.1038/s41524-021-00559-9 | NPT对比Mishin EAM、ACE、SNAP、GTINV；文中报告600 K以下热膨胀与实验较好，高温存在小偏离 | 最直接的多模型比较；复杂模型的综合优势不等于目标温区膨胀必然更准 |
| 6 | Olsson et al., 2023. Grain Size-Dependent Thermal Expansion of Nanocrystalline Metals. https://doi.org/10.3390/ma16145032 ; https://lup.lub.lu.se/search/files/159582754/Olsson_2023c.pdf | LAMMPS、Mishin势、单晶/非晶/纳米晶与DFT；单晶/非晶约173–798 K，每25 K一点，每点平衡1 ns、采样1 ns；三次多项式处理体积 | 最接近的模拟流程；最高约20%的增强是晶粒尺寸效应，不是预测误差 |

## 实验与参考数据

| 编号 | 文献与链接 | 信息 | 使用边界 |
|---|---|---|---|
| 7 | Pojur & Yates, 1973. Thermal expansion at elevated temperatures. I. Apparatus for use in the temperature range 300–800 K: the thermal expansion of copper. https://doi.org/10.1088/0022-3735/6/1/030 | 温区与本课题直接匹配的实验工作 | 已核对书目信息；数值表应在取得原文后提取 |
| 8 | Wang & Reeber, 1996. Thermal Expansion of Copper. High Temperature and Materials Science 35, 181–186. https://www.researchgate.net/profile/Robert-Reeber/publication/235924895_Thermal_Expansion_of_Copper/links/59e2ec2a458515393d59f7be/Thermal-Expansion-of-Copper.pdf | 基于数据和半经验热力学模型的宽温区推荐结果 | 区分推荐值与直接实验测量；2021年PACE论文也采用此参考 |
| 9 | Fornasini et al., 2004. Extended x-ray-absorption fine-structure measurements of copper: Local dynamics, anharmonicity, and thermal expansion. https://doi.org/10.1103/PhysRevB.70.174301 | 4–500 K EXAFS；局域第一近邻膨胀与晶体学膨胀存在区别 | 不应直接用径向分布函数峰位或平均键长替代晶格常数 |

## 势函数横向比较

| 编号 | 文献与链接 | 主要目标 | 本项目判断 |
|---|---|---|---|
| 10 | Foiles, Baskes & Daw, 1986. Embedded-atom-method functions for the fcc metals Cu, Ag, Au, Ni, Pd, Pt, and their alloys. https://doi.org/10.1103/PhysRevB.33.7983 | 基于凝聚能、晶格常数、弹性与空位等性质的经典EAM | 可作早期经验势对照；静态拟合不保证膨胀准确 |
| 11 | Mishin et al., 2001. Structural stability and lattice defects in copper: Ab initio, tight-binding, and embedded-atom calculations. https://doi.org/10.1103/PhysRevB.63.224106 | 铜结构稳定性与缺陷；EAM1含第一性原理拟合信息 | 主势；在NIST下载Cu01.eam.alloy并记录版本与哈希 |
| 12 | Mendelev et al., 2008. Analysis of semi-empirical interatomic potentials appropriate for simulation of crystalline and liquid Al and Cu. https://doi.org/10.1080/14786430802206482 | 改善液体结构及凝固应用；相同静态性质的势可能产生不同动力学 | 若增加势函数对照，采用明确的Cu1.eam.fs版本；不是同名2012版本 |
| 13 | Sheng et al., 2011. Highly optimized embedded-atom-method potentials for fourteen fcc metals. https://doi.org/10.1103/PhysRevB.83.134118 | 第一性原理势能面及实验校准，验证热膨胀、熔化、液体动力学 | 与热学任务相关的候选势；需核查铜文件具体来源与格式 |
| 14 | Asadi et al., 2015. Two-phase solid-liquid coexistence of Ni, Cu, and Al by molecular dynamics simulations using the modified embedded-atom method. https://doi.org/10.1016/j.actamat.2014.12.010 | 2NN-MEAM与固液共存 | 高温/熔化取向；当前温区不必优先增加复杂度 |
| 15 | Etesami & Asadi, 2018. Molecular dynamics for near melting temperatures simulations of metals using modified embedded-atom method. https://doi.org/10.1016/j.jpcs.2017.09.001 | 针对近熔点性能改进MEAM，包含铜 | 证明势的适用性与训练目标有关 |
| 16 | Song et al., 2024. General-purpose machine-learned potential for 16 elemental metals and their alloys. https://doi.org/10.1038/s41467-024-54554-x | 包含Cu的UNEP-v1，验证多种结构和合金行为 | 作为前沿背景；不能直接推断比Mishin势更准确预测本温区膨胀 |

## 系综统计与软件依据

| 编号 | 文献与链接 | 本项目使用方式 |
|---|---|---|
| 17 | Martyna, Tobias & Klein, 1994. Constant pressure molecular dynamics algorithms. https://doi.org/10.1063/1.467468 | 恒压采样理论；具体实现以LAMMPS文档为准 |
| 18 | Shinoda, Shiga & Mikami, 2004. Rapid estimation of elastic constants by molecular dynamics simulation under constant stress. https://doi.org/10.1103/PhysRevB.69.134103 | LAMMPS温压控制算法相关依据；原文测试对象包含Ni，不是铜热膨胀基准 |
| 19 | Flyvbjerg & Petersen, 1989. Error estimates on averages of correlated data. https://doi.org/10.1063/1.457480 | 处理轨迹时间相关性，使用分块统计 |
| 20 | Grossfield et al., 2019. Best Practices for Quantification of Uncertainty and Sampling Quality in Molecular Simulations. https://doi.org/10.33011/livecoms.1.1.5067 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC6286151/ | 平衡、重复、自相关与采样误差方法指南 |
| 21 | Thompson et al., 2022. LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales. https://doi.org/10.1016/j.cpc.2021.108171 | 软件引用 |

## 从比较得到的方法选择

早期研究建立经验势MD求热学性质的流程，第一性原理与紧束缚工作补充自由能、声子与近似误差分析，后续研究拓展到多势比较和晶粒尺寸。该定向文献集不能用于推断全领域发表量随时间的趋势。

最值得优先阅读的顺序：Olsson 2023的方法、Lysogorskiy 2021图4及其方法、Mishin 2001及NIST说明、Wang与Reeber 1996参考数据、Grossfield 2019统计指南。

主要比较维度是目标物性、材料状态、势的拟合目标、有限温度处理及不确定性。没有统一协议重跑各模型，不能给出所有势的误差排行榜。也没有可靠的引用次数数据，不给出“最高引用论文/作者”排名。本次文献主要围绕热膨胀、铜、多体势、分子动力学和有限温度自由能展开。

## 建议模拟路线

以下是课程项目建议参数，不是某一篇论文的原样复制，也尚未实际运行验证。

1. Python调用LAMMPS共享库，Jupyter整理结果。先串行跑通，再考虑并行。Python包与共享库需匹配，Notebook内核需属于同一环境。
2. 三维FCC单晶、[100]/[010]/[001]与盒轴平行、10×10×10常规晶胞、4000原子、三向周期性边界。该模型是体相近似，不是自由纳米颗粒。
3. 使用metal单位制，Cu质量63.546，Mishin EAM1，pair_style eam/alloy，pair_coeff * * Cu01.eam.alloy Cu。初始晶格常数约3.615 Å仅作建模初值，能量和盒子共同弛豫后再作有限温度计算。
4. 300、400、500、600、700、800 K各自独立运行，3个独立初速度种子。预试验稳定后，可扩展到50 K间隔。
5. NPT，iso 0 0 1.0；1 fs时间步（metal下0.001 ps）；温度阻尼0.1 ps，压力阻尼1.0 ps。这些是起始设置，需检查振荡和收敛。npt本身积分，不对相同原子叠加nve。
6. 预试验先平衡100 ps、采样200 ps。正式计算先按平衡200 ps、采样500 ps安排；300与800 K延长到平衡1 ns、采样1 ns作时长检查。若短长结果不一致或仍漂移，延长全部相关温度点。先前50 ps加100 ps不能直接称为已验证充分。
7. 采样阶段每0.1 ps输出时间、温度、压力、势能、Lx/Ly/Lz、体积。保存未经长时间平均的标量轨迹，结构每5–10 ps一帧即可。记录种子、版本、参数与势文件校验值。
8. 晶格常数由每帧Lx/10计算后时间平均。相对膨胀以300 K为参考。主结果为300–800 K工程平均线膨胀系数；局部温度系数只作补充，不强求六点三次拟合。
9. 对轨迹做分块平均，增加块长观察均值标准误是否稳定；3次独立重复分别计算系数，再报告均值及重复间标准差，并明确误差棒定义。不要将每个时间步当作独立样本。
10. 300与800 K检查10³和12³晶胞（4000与6912原子）、1与0.5 fs（物理时长一致）、短长采样时长；结果差异若超过统计不确定性，则继续扩大尺寸/延长时长。可把平均系数相对变化约1–2%作为课程项目预设目标，不能称为通用标准。

## 系数定义与实验对照

工程温区平均：alpha_bar = [a(800)-a(300)]/[a(300)×500]，单位K^-1。

瞬时热力学线系数：alpha(T) = (1/a(T)) da/dT。部分文献用固定参考长度a_ref作分母，比较前需统一。立方体系瞬时体系数beta(T)=3alpha(T)；有限温区体积相对变化不可机械地精确除以3。

先对照epsilon(T)=[a(T)-a(300)]/a(300)，再对照同定义平均系数。Wang与Reeber是宽温区评价参考，Pojur与Yates温区匹配但需要原始表格。以图读数时保留数字化误差。不得用室温单点16–17×10^-6 K^-1的常见量级替代300–800 K实验平均值。

## Python执行结构

Notebook组织为背景、文献、参数、环境检查、建模、试运行、正式运行、统计、实验对照、讨论。Python负责按温度与种子调用同一LAMMPS模板，每次输出独立文件夹，保存原始日志与汇总表。首阶段6温度×3种子=18条主轨迹，可靠性检查另计；扩展11温度×3种子=33条主轨迹。运行时间先由短试跑测量，不预先承诺几分钟完成。

```python
from lammps import lammps

temperatures = [300, 400, 500, 600, 700, 800]
seeds = [13579, 24681, 97531]

for temperature in temperatures:
    for seed in seeds:
        # 此处创建独立输出目录并生成该任务的输入文件
        lmp = lammps()
        try:
            lmp.file(input_file_for(temperature, seed))
        finally:
            lmp.close()
```

该段仅展示组织方式，input_file_for是待实现的项目函数，不是可直接运行的完整程序。

## 官方实现文档

- Python组件与共享库：https://docs.lammps.org/Python_install.html
- Python调用方式：https://docs.lammps.org/Howto_python.html
- 温压控制：https://docs.lammps.org/fix_nh.html
- 势函数调用：https://docs.lammps.org/pair_eam.html
- 势文件来源：https://www.ctcms.nist.gov/potentials/entry/2001--Mishin-Y-Mehl-M-J-Papaconstantopoulos-D-A-et-al--Cu-1/

最终研究定位：在明确的模型与统计条件下评价Mishin EAM1对300–800 K体相单晶铜晶格膨胀和平均线膨胀系数的预测。前期调研时尚未生成结果；截至2026-10-18，24条轨迹与全部统计核查已完成，平均系数17.8883×10⁻⁶ K⁻¹，对NIST同定义参考偏差−3.8463%。完整结果见Notebook与results目录。


## 本轮补充与实际采用方法（2026-10-17）

以上是前期方法调研；下面记录本轮实际实现和数据来源，发生冲突时以此处与 `run_study.py` 为准。

|编号|补充来源|证据与用途|
|---|---|---|
|22|Hahn, 1970, Thermal Expansion of Copper from 20 to 800 K—Standard Reference Material 736，[DOI](https://doi.org/10.1063/1.1658614)|铜热膨胀标准参考物质的原始研究；最终数值采用后续修订证书|
|23|[NIST SRM 736历史证书，1990-10-07修订](https://tsapps.nist.gov/srmext/certificates/archives/736.pdf)|直接核对扫描件第1页，录入300、400、500、600、700、800 K表值；统一换算到300 K后比较|
|24|Edwards, Tildesley & Binsted, 1997, Cumulant expansion analysis of thermal disorder in fcc copper，[DOI](https://doi.org/10.1080/00268979709482724)|补充检索得到的相关热无序研究；仅核对检索元数据，未据此提取具体参数或数值|
|25|Sadigh & Grimvall, 1996，[DOI](https://doi.org/10.1103/PhysRevB.54.15742)|液态铜相关工作，属于范围边界；不能作为本项目固态300–800 K结果的直接验证|

### 横向比较后为什么这样做

|路线|能回答什么|代价或限制|本项目安排|
|---|---|---|---|
|EAM + 零外压NPT MD|有限温度下晶胞体积、晶格常数和非谐膨胀|经验势系统误差；经典原子统计|主路线；Mishin EAM1，便于课程内复现|
|DFT/DFPT + 准谐自由能|从自由能最小值获得热膨胀|交换关联近似、计算成本、准谐近似|文献比较，不声称本项目做了DFT|
|紧束缚 + 准谐方法|连接电子结构与热力学|仍有拟合和准谐误差|用于解释不同方法的误差来源|
|ACE、SNAP及其他机器学习势|更灵活描述多种局域环境|训练范围、部署和成本；复杂不保证本指标更准|作为方法背景；没有实际运行这些势|
|单晶与纳米晶比较|晶界对热膨胀的贡献|需要不同晶粒结构及更多采样|不纳入当前理想单晶模型|
|膨胀仪与X射线/EXAFS实验|宏观长度、晶格或局域键的热响应|各测量量定义不同|以宏观长度变化作同定义参考；不把局域键长当晶格常数|

2023年工作采用更长的逐温度采样，并研究不同微结构。本项目只研究理想单晶：6个温度点×3次独立初速度，共18条正式轨迹，每条200 ps平衡+500 ps采样。在300和800 K分别进行尺寸（4000→6912原子）、步长（1→0.5 fs）和时长（1 ns平衡+1 ns采样）检查，共6条验证轨迹。上述时长是本项目选定参数，并非宣称逐字复现该论文。验证相对基准r1的平均系数变化以±2%作为预设实用阈值，未通过时不称作已收敛。

统计使用1–50 ps分块敏感性、前后半段比较和独立重复。最终误差棒为3次重复的t统计区间，只有2个自由度，不能代表势函数、经典近似或实验参考的不确定度。端点工程平均为主结果；不对6点数据过度高阶拟合后宣称精确瞬时系数。

### 检索透明度

本轮paper-search经OpenAlex与Crossref返回6个去重条目；arXiv及Semantic Scholar出现429，DBLP解析失败。新增条目中，暖稠密铜、铜溴化物自旋体系及介电综述不满足固态铜热膨胀主题，未用于论证。它们不是纳入研究，更不能用来凑“研究数量”。这是一份定向文献调研，不是覆盖所有数据库的系统综述；没有做引用量排序。文献24、25仅作补充线索，证据强度低于已阅读的相关原文。

Narasimhan与de Gironcoli应按**2002年正式发表**引用：[Physical Review B 65, 064302](https://doi.org/10.1103/PhysRevB.65.064302)，2001是预印本年份。前期列表保留预印本入口方便获取原文。

原先考虑Wang–Reeber或Pojur–Yates参考；未取得可核对的原始数值表，因此没有从摘要编造数值。最终采用可公开核对的NIST修订证书。完整录入表、单位与参考温度换算见 `data/reference/SOURCE.md`。
