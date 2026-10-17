# 单晶铜热膨胀行为及线膨胀系数的分子动力学研究

材料学院计算材料学课程项目。问题很简单：**把理想单晶铜从300 K加热到800 K，它的晶格和体积怎样变化？Mishin EAM1势算出的膨胀系数与实验差多少？**

## 当前进展

正式及验证轨迹完成 **0/24**。 **模拟尚未完成，当前Notebook是方法稿，不能作为最终论文提交。** 现有两温度试跑不替代正式数据。

阶段1环境、阶段2建模、阶段3试跑已完成。阶段4正式模拟与阶段5验证按下列清单推进；阶段6统计和阶段7 Notebook论文必须等数据完成并核查后才能定稿。PPT、视频及答辩属于后续阶段，本轮没有编造这些成果。

## 队友先看哪里

1. 点开[Notebook论文](铜热膨胀研究.ipynb)：GitHub可以直接阅读文字和已有输出，不必安装软件。当前状态以本文和Notebook中的标记为准。
2. 阅读[文献对比](docs/literature_review.md)：包括25项参考来源及范围说明，解释为什么用EAM、NPT、独立重复和参数检查。
3. 阅读[实验数据说明](data/reference/SOURCE.md)：实验原表基准是293 K，必须换算到300 K再比较。
4. 看[项目推进方案](docs/project_plan.md)了解一周安排；看[阶段进展](阶段进展.md)了解实际进度。

Notebook就是把“论文文字、读取数据的代码和图表”放在一个文件里。队友在VS Code打开`.ipynb`文件即可阅读；需要执行单元格时，右上角选择本项目Python环境。运行Notebook只展示结果，不会自动重新启动长模拟。

## 文件目录树

```text
copper-thermal-expansion/
├── README.md                      项目入口和运行说明
├── 铜热膨胀研究.ipynb               Notebook论文/方法稿，见进度标记
├── 阶段进展.md                    按任务记录实际进度
├── docs/
│   ├── literature_review.md        全部参考文献、方法横向比较
│   └── project_plan.md             一周推进方案
├── data/reference/
│   ├── nist_srm736.csv             人工核对的实验表值
│   └── SOURCE.md                   来源、单位和换算公式
├── potentials/
│   ├── Cu01.eam.alloy              Mishin 2001 EAM1铜势
│   └── SOURCE.md                   下载来源和校验值
├── outputs/
│   ├── model/                     初始/弛豫结构、模型图与检查
│   ├── trials/                    300、800 K试跑，非正式研究数据
│   ├── production/                18条正式轨迹的独立文件夹
│   ├── validation/                6条端点验证的独立文件夹
│   └── study_manifest.json         24项任务的参数清单
├── results/                       全部计算完成后才生成
│   ├── run_summary.csv            单条轨迹统计
│   ├── temperature_summary.csv    按温度汇总
│   ├── alpha_replicates.csv        各次重复的平均系数
│   ├── validation_comparison.csv  参数检查结果
│   ├── metrics.json               论文使用的主要数值
│   └── figures/                   正式图表
├── config.json                    建模及试跑参数
├── md_project.py                  环境检查、建模和试跑
├── run_study.py                   正式及验证模拟
├── analyze_study.py               读取完整数据、统计和绘图
├── write_paper.py                 根据真实数据生成论文
├── execute_notebook.py            执行论文中的展示单元格
├── write_readme.py                按真实文件更新进度
├── finish_project.py              等待模拟后分析、执行论文、上传
├── requirements.txt               主要Python依赖版本
└── requirements-lock.txt          本机Windows完整依赖记录
```

目录树包含计划生成项；尚未存在的`results/`不代表结果已完成。每条模拟文件夹内有`input.in`（LAMMPS输入）、`thermo.dat`（原始时间序列）、`log.lammps`（运行日志）、`metadata.json`（参数、种子、哈希与完成状态）。结构轨迹和二进制重启文件留在本机，避免仓库过大；公开的原始标量与脚本足以重新分析并复算。

## 实际计算条件

|项目|设置与含义|
|---|---|
|模型|FCC铜，10×10×10常规晶胞，4000原子；三向周期边界，无表面/晶界|
|势|Mishin等2001 EAM1，NIST Cu01.eam.alloy，eam/alloy|
|温度|300、400、500、600、700、800 K，每点3次独立初速度|
|系综|NPT；目标压力0 bar；各向同性盒长调整|
|步长与阻尼|1 fs；温度阻尼0.1 ps，压力阻尼1 ps|
|正式时长|每条200 ps平衡+500 ps采样；只统计平衡之后|
|输出|标量每0.1 ps；结构每100 ps|
|尺寸检查|300/800 K各一条6912原子轨迹|
|步长检查|300/800 K各一条0.5 fs轨迹，保持物理时长|
|时长检查|300/800 K各1 ns平衡+1 ns采样|
|统计|1–50 ps分块敏感性；前后半段比较；3个重复的t统计区间|
|软件|LAMMPS 22Jul2025 update4，双精度OpenMP；Python调用|

温区平均系数为`[a(800)-a(300)]/[a(300)×500]`，单位K⁻¹。不是室温单点瞬时系数。验证以相对基准r1的平均系数变化±2%作为预设实用目标，并结合轨迹漂移、分块统计判断；它不是普适标准。3次重复较少，误差棒只反映统计误差，不代表经验势模型准确度。

## 在自己的电脑运行（Windows / VS Code）

使用Python 3.12，先打开项目文件夹，再在VS Code终端运行：

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m ipykernel install --prefix .venv --name copper-md --display-name "Copper MD"
```

Windows的LAMMPS wheel还需要微软MS-MPI运行库：[微软官方入口](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi)。本机已配置项目内运行库；二进制不放进Git。队友可安装官方运行库，代码也支持系统搜索路径。不要复制别人电脑的`.venv`。Linux需要与wheel匹配的MPI；本项目实际验证的是Windows环境。

按顺序执行：

```powershell
.\.venv\Scripts\python.exe md_project.py environment
.\.venv\Scripts\python.exe md_project.py build
.\.venv\Scripts\python.exe run_study.py batch --workers 2
.\.venv\Scripts\python.exe analyze_study.py
.\.venv\Scripts\python.exe write_paper.py
.\.venv\Scripts\python.exe execute_notebook.py
.\.venv\Scripts\python.exe write_readme.py
```

完整模拟可能运行数小时，取决于电脑。不要在同一目录同时重复启动批量任务。已完成任务自动跳过，未完成旧记录会触发保护并停止；断电后的不完整任务需要先另目录归档，核查原因后重跑，不能直接删除证据或标记完成。`finish_project.py`可在已经启动批量计算后等待并衔接后续步骤；只有加`--publish`才会提交和推送到main。

若只想阅读/重画图，直接执行Notebook，不必重跑模拟。重新生成论文会覆盖Notebook正文，应先保留手动修改。

## 可复现性与研究边界

势文件SHA256：`5e374260448eea5c26de59432d75e702c1778cc6999d2527a789e5387ce9dc9e`。脚本记录每条任务输入哈希、实际种子、软件版本、原子数和运行时间。无完整24条数据时分析程序拒绝生成正式结果。失败输出和硬件试验不进入正式统计。

原始课程课件、他人论文PDF、虚拟环境、GPU/MPI安装程序及个人环境路径记录不上传。势文件来源单独注明；本仓库未给第三方资料重新授予许可。项目不声称发现新的铜热膨胀规律，也不声称机器学习势已被本项目实际对比运行。
