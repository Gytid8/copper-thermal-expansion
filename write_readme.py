"""Keep the public progress statement derived from files, not assumptions."""
from pathlib import Path
import json
from datetime import datetime,timezone
from run_study import manifest,folder
ROOT=Path(__file__).resolve().parent

def main():
    states=[]
    for j in manifest():
        p=folder(j)/'metadata.json'
        try:
            s=json.loads(p.read_text(encoding='utf8')).get('status') if p.exists() else 'pending'
        except (json.JSONDecodeError,PermissionError):
            s='updating'
        states.append((j['id'],s))
    done=sum(s=='completed' for _,s in states)
    status=f'正式及验证轨迹完成 **{done}/24**。'
    mpath=ROOT/'results/metrics.json'
    if mpath.exists() and done==24:
        m=json.loads(mpath.read_text(encoding='utf8'))
        status+=f"\n\n300–800 K平均线膨胀系数：**{m['alpha_300_800_Kinv']*1e6:.4f} ×10⁻⁶ K⁻¹**；95%统计区间半宽{m['alpha_ci95_halfwidth_Kinv']*1e6:.4f} ×10⁻⁶ K⁻¹。NIST同定义参考{m['nist_alpha_300_800_Kinv']*1e6:.4f} ×10⁻⁶ K⁻¹；偏差{m['relative_error_percent']:+.2f}%。三项端点验证通过：**{m['all_practical_checks_pass']}**。"
        if not m['all_practical_checks_pass']:
            status+='\n\n**参数验证尚未全部通过，不能称作已收敛，需要进一步补算。**'
    else:
        status+=' **模拟尚未完成，当前Notebook是方法稿，不能作为最终论文提交。** 现有两温度试跑不替代正式数据。'
    if mpath.exists() and done==24:
        import pandas as pd
        check=pd.read_csv(ROOT/'results/validation_comparison.csv')
        status+='\n\n|验证项目|平均系数相对基准r1的变化|\n|---|---|\n'
        labels={'size12':'尺寸：4000→6912原子','dt05':'时间步：1→0.5 fs','long':'时长：延长至1+1 ns'}
        for _,r in check.iterrows():status+=f"|{labels[r['check']]}|{r.change_percent:+.4f}%|\n"
        status+='\n![统一300 K参考后的实验对照](results/paper_figures/figure3_reference_composite_v1_preview.png)\n\n英文论文已整合**6张新版图、2张正文表和12条参考文献**，图文引用与编号已同步更新。正文采用研究论文结构，可执行代码集中在Notebook附录。[双栏PDF阅读版](results/copper_thermal_expansion_journal.pdf)适合阅读与打印；[HTML阅读版](results/copper_thermal_expansion_paper.html)下载后用浏览器打开，附录代码可折叠。小组成员姓名、单位与实际分工仍需组员填写。'
    readme='''# 单晶铜热膨胀行为及线膨胀系数的分子动力学研究

材料学院计算材料学课程项目。问题很简单：**把理想单晶铜从300 K加热到800 K，它的晶格和体积怎样变化？Mishin EAM1势算出的膨胀系数与实验差多少？**

## 当前进展

'''+status+'''

阶段1至7的环境、建模、模拟、验证、统计与Notebook写作已完成。论文按“引言—计算方法—结果与讨论—结论”组织，参考文献已逐项核对，详见引用审阅记录。PPT、视频及答辩准备属于后续工作。

## 队友先看哪里

1. 点开[Notebook论文](铜热膨胀研究.ipynb)：GitHub可以直接阅读文字和已有输出，不必安装软件。当前状态以本文和Notebook中的标记为准。
2. 阅读[文献对比](docs/literature_review.md)：包括25项参考来源及范围说明，解释为什么用EAM、NPT、独立重复和参数检查。
3. 阅读[实验数据说明](data/reference/SOURCE.md)：实验原表基准是293 K，必须换算到300 K再比较。
4. 看[项目推进方案](docs/project_plan.md)了解一周安排；看[阶段进展](阶段进展.md)了解实际进度。
5. 阅读[纯文字稿](docs/manuscript.md)、[写作参照比较](docs/paper_writing_comparison.md)和[引用与统计审阅](docs/reference_audit.md)，核对正文依据。

Notebook就是把“论文文字、读取数据的代码和图表”放在一个文件里。队友在VS Code打开`.ipynb`文件即可阅读；需要执行单元格时，右上角选择本项目Python环境。运行Notebook只展示结果，不会自动重新启动长模拟。

## 文件目录树

```text
copper-thermal-expansion/
├── README.md                      项目入口和运行说明
├── 铜热膨胀研究.ipynb               Notebook论文/方法稿，见进度标记
├── 阶段进展.md                    按任务记录实际进度
├── docs/
│   ├── literature_review.md        全部参考文献、方法横向比较
│   ├── project_plan.md             一周推进方案
│   ├── recovery_record.md          中断恢复、边界清理与最终核查
│   ├── manuscript.md               六图英文论文正文、引用与附录源稿
│   ├── figure_manifest.json        六张选定图、图注、来源和校验值
│   ├── paper_layout_review.md     排版与整合核查记录
│   ├── paper_style.css             标题、正文、图注和表格的统一样式
│   ├── manuscript_revision_english.md 英文修订及检查记录
│   ├── reference_verification_20261019.json 本轮引用核查
│   ├── paper_writing_comparison.md  同类论文结构与写法对比
│   ├── reference_audit.md          本轮引用与统计审阅说明
│   └── reference_verification.json 实时文献核查记录
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
│   ├── copper_thermal_expansion_journal.pdf 双栏论文阅读版
│   ├── copper_thermal_expansion_paper.html 六图HTML阅读版
│   ├── run_summary.csv            单条轨迹统计
│   ├── temperature_summary.csv    按温度汇总
│   ├── alpha_replicates.csv        各次重复的平均系数
│   ├── validation_comparison.csv  参数检查结果
│   ├── metrics.json               论文使用的主要数值
│   ├── figures/                   原始分析图表
│   └── paper_figures/              六张正式图、PDF/PNG和逐图数据核查
├── config.json                    建模及试跑参数
├── md_project.py                  环境检查、建模和试跑
├── run_study.py                   正式及验证模拟
├── resume_study.py                从平衡重启点恢复完整采样
├── prepare_recovery.py            本次中断归档与恢复准备记录
├── repair_completed_recovery.py   完成输出的重复边界严格核查
├── show_progress.py               只读实时进度窗口
├── analyze_study.py               读取完整数据、统计和绘图
├── paper_layout/main.tex         双栏PDF的可编辑排版源文件
├── plot_figure2_v2.py             图2：趋势、波动与体积矩阵
├── plot_figure3_v1.py             图3：实验对照与残差
├── plot_figure4_v1.py             图4：设置矩阵与参数敏感性
├── plot_figure5_v1.py             图5：3D轨迹与端点诊断
├── plot_figure6_v1.py             图6：分块统计敏感性
├── export_journal_pdf.py          导出双栏论文PDF（需要XeLaTeX）
├── export_paper_html.py           导出六图HTML阅读版
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
.\\.venv\\Scripts\\python.exe -m pip install -r requirements.txt
.\\.venv\\Scripts\\python.exe -m ipykernel install --prefix .venv --name copper-md --display-name "Copper MD"
```

Windows的LAMMPS wheel还需要微软MS-MPI运行库：[微软官方入口](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi)。本机已配置项目内运行库；二进制不放进Git。队友可安装官方运行库，代码也支持系统搜索路径。不要复制别人电脑的`.venv`。Linux需要与wheel匹配的MPI；本项目实际验证的是Windows环境。

按顺序执行：

```powershell
.\\.venv\\Scripts\\python.exe md_project.py environment
.\\.venv\\Scripts\\python.exe md_project.py build
.\\.venv\\Scripts\\python.exe run_study.py batch --workers 2
.\\.venv\\Scripts\\python.exe analyze_study.py
.\\.venv\\Scripts\\python.exe write_paper.py
.\\.venv\\Scripts\\python.exe execute_notebook.py
.\\.venv\\Scripts\\python.exe export_paper_html.py
.\\.venv\\Scripts\\python.exe write_readme.py
```

完整模拟可能运行数小时，取决于电脑。不要在同一目录同时重复启动批量任务。已完成任务自动跳过，未完成旧记录会触发保护并停止；断电后的不完整任务需要先另目录归档，核查原因后重跑，不能直接删除证据或标记完成。`finish_project.py`可在已经启动批量计算后等待并衔接后续步骤；只有加`--publish`才会提交和推送到main。

直接阅读已保存的Notebook或PDF不必安装环境，也不必重跑模拟。Notebook附录默认核查数据而不重画图；需要重画时将附录的REBUILD_FIGURES设为True。重新生成论文会覆盖Notebook正文，应先保留手动修改。双栏PDF可运行export_journal_pdf.py重建，需要XeLaTeX；本机编译路径配置见paper_layout/README.md。

## 可复现性与研究边界

2026-10-17两条步长验证曾被外部中断，现从200 ps平衡重启点重新采满500 ps，原中断记录完整归档。恢复方法及核查见[中断与恢复记录](docs/recovery_record.md)。已完成的20条轨迹不重算。

势文件SHA256：`5e374260448eea5c26de59432d75e702c1778cc6999d2527a789e5387ce9dc9e`。脚本记录每条任务输入哈希、实际种子、软件版本、原子数和运行时间。无完整24条数据时分析程序拒绝生成正式结果。失败输出和硬件试验不进入正式统计。

原始课程课件、他人论文PDF、虚拟环境、GPU/MPI安装程序及个人环境路径记录不上传。势文件来源单独注明；本仓库未给第三方资料重新授予许可。项目不声称发现新的铜热膨胀规律，也不声称机器学习势已被本项目实际对比运行。
'''
    (ROOT/'README.md').write_text(readme,encoding='utf8')
    progress='# 阶段进展\n\n'+status+'\n\n最后核查时间（UTC）：'+datetime.now(timezone.utc).isoformat()+'\n\n|任务|状态|\n|---|---|\n'+'\n'.join(f'|{n}|{s}|' for n,s in states)
    (ROOT/'阶段进展.md').write_text(progress,encoding='utf8')

if __name__=='__main__':main()
