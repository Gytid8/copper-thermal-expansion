"""Generate a readable paper; numerical conclusions require complete verified data."""
from pathlib import Path
import json
import nbformat as nb
import pandas as pd
ROOT=Path(__file__).resolve().parent

def main():
    complete=(ROOT/'results/metrics.json').exists()
    cells=[]
    def md(s):cells.append(nb.v4.new_markdown_cell(s))
    def code(s):cells.append(nb.v4.new_code_cell(s))
    md('# 单晶铜热膨胀行为及线膨胀系数的分子动力学研究\n\n小组成员：【请填写姓名及实际分工】\n\n计算材料学课程项目。本文由原始计算数据生成；研究定位为方法复现与势函数适用性评价。')
    md('## 摘要\n\n采用LAMMPS和Mishin等人2001年铜EAM1势，对三向周期性FCC铜进行近零外压NPT分子动力学模拟，研究300–800 K的晶格常数与体积。设置6个温度、每点3次独立重复，并在温区端点检查尺寸、时间步与时长。实验参考采用NIST SRM 736的1990修订表，统一以300 K长度为参考。')
    if complete:
        m=json.loads((ROOT/'results/metrics.json').read_text(encoding='utf8'))
        md(f"计算结果：温区平均线膨胀系数为 **{m['alpha_300_800_Kinv']*1e6:.4f} ×10⁻⁶ K⁻¹**，95%统计区间半宽为{m['alpha_ci95_halfwidth_Kinv']*1e6:.4f} ×10⁻⁶ K⁻¹；NIST同定义参考为{m['nist_alpha_300_800_Kinv']*1e6:.4f} ×10⁻⁶ K⁻¹，相对偏差{m['relative_error_percent']:+.2f}%。尺寸、时间步和时长三项检查{'全部满足预设实用阈值' if m['all_practical_checks_pass'] else '存在未通过项'}。统计区间不包括模型误差。")
    else:
        md('**当前为论文方法稿：正式模拟尚未全部完成，结果、摘要数值和最终结论暂不填写。现有试跑仅用于运行检查，不作为18条正式轨迹的替代。**')
    md('## 1. 研究背景与文献比较\n\n铜受热后的平均原子位置间距增加，反映势能的非谐性。谐振动模型本身不能给出普通热膨胀。NPT允许体积随温度调整，适合获得近零外压的平衡晶格参数。\n\nAlper与Politzer（1999）以及Kazanc等（2006）展示经验势有限温度研究；Narasimhan与de Gironcoli（2002）、Rudin等（2002）采用电子结构与准谐自由能路线。Lysogorskiy等（2021）直接比较铜的多类势；Olsson等（2023）用Mishin势研究单晶及晶界影响。它们支持本项目路线，但不能据此宣称我们已获得相同精度。\n\n本项目保留原选题的Mishin EAM1，不额外开展机器学习势或DFT计算。详细的25项来源、阅读证据和方法对比见 [文献调研](docs/literature_review.md)。')
    md('## 2. 模型与计算方法\n\nFCC常规晶胞沿[100]、[010]、[001]排列，10×10×10晶胞共4000原子，三向周期性边界，无自由表面和晶界。初始a=3.615 Å，先做原子和晶胞能量最小化。Cu质量63.546 u。势文件为NIST `Cu01.eam.alloy`，调用`eam/alloy`；来源与哈希见potentials/SOURCE.md。\n\n使用metal单位制：长度Å、时间ps、能量eV、压力bar。NPT设温300、400、500、600、700、800 K，外压0 bar，等比例调整盒长；时间步0.001 ps（1 fs），温度阻尼0.1 ps，压力阻尼1 ps。每点200 ps平衡、500 ps采样，间隔0.1 ps输出未作长时间平均的标量，结构每100 ps输出。\n\n3个重复的种子基数为13579、24681、97531，实际种子=基数+19×温度，每个温度独立初始化。端点各增加12³晶胞、0.5 fs步长和1 ns平衡+1 ns采样三项检查，总计24条轨迹。所有具体参数见run_study.py及outputs/study_manifest.json。')
    md('### 2.1 数据处理\n\n每帧由a=Lx/n取得晶格常数，再对平衡期之后的样本求均值。立方盒保持等比例变化，体积同时从LAMMPS直接读取。主结果定义为\n\n$$\\bar\\alpha_{300-800}=\\frac{\\bar a(800)-\\bar a(300)}{500\\bar a(300)}.$$\n\n这是工程温区平均，区别于瞬时系数$\\alpha(T)=a(T)^{-1}da/dT$。不把体膨胀系数直接当作线膨胀系数，也不把局域平均键长等同于晶格常数。\n\n每条轨迹比较前后半段，采用1、2、5、10、20、50 ps分块检查时间相关性，以20 ps分块报告单轨迹均值标准误。每温度汇总3条独立轨迹；系数的不确定度以3个重复端点系数的标准误和t分布（自由度2）得到。重复数较少，统计区间需谨慎解释。\n\n端点验证与基准r1比较，预设平均系数变化±2%为课程项目的实用目标，不能视作普适收敛标准。单次检查不足以精确量化所有系统误差。')
    md('### 2.2 实验参考与公平比较\n\n采用[NIST历史修订证书](https://tsapps.nist.gov/srmext/certificates/archives/736.pdf)第1页表值，原表参考293 K，换算为$\\epsilon_{300}(T)=[1+e_{293}(T)]/[1+e_{293}(300)]-1$。同定义300–800 K参考平均约18.604×10⁻⁶ K⁻¹。原始表与转换说明见[data/reference/SOURCE.md](data/reference/SOURCE.md)。实验为高纯铜试样，MD为理想周期晶体；两者差别属于讨论内容。')
    code("from pathlib import Path\nimport json\nimport pandas as pd\nfrom IPython.display import display, Image\nroot=Path.cwd()\nassert (root/'config.json').exists(), '请在项目根目录打开Notebook'\ndisplay(pd.read_csv('data/reference/nist_srm736.csv'))")
    md('## 3. 计算结果与可靠性检查')
    code("from run_study import manifest, folder\nrecords=[]\nfor j in manifest():\n    p=folder(j)/'metadata.json'\n    state=json.loads(p.read_text(encoding='utf8')).get('status') if p.exists() else 'pending'\n    records.append({'任务':j['id'],'状态':state,'温度/K':j['temperature']})\ndisplay(pd.DataFrame(records))")
    if complete:
        summary=pd.read_csv(ROOT/'results/temperature_summary.csv')
        checks=pd.read_csv(ROOT/'results/validation_comparison.csv')
        runs=pd.read_csv(ROOT/'results/run_summary.csv')
        low=summary.loc[summary.temperature==300].iloc[0];high=summary.loc[summary.temperature==800].iloc[0]
        md(f"### 3.1 晶格和体积的温度响应\n\n300 K的晶格常数为{low.a_A:.7f} Å，800 K为{high.a_A:.7f} Å。六个温度点的平均晶格常数随温度单调增加；相对300 K，800 K线膨胀为{high.strain_from_300*100:.4f}%，体积增加{high.relative_volume*100:.4f}%。这里的体积是瞬时体积的时间平均，不强制令平均体积等于平均晶格常数的立方。\n\n两者方向一致，符合正热膨胀。实验归一化长度曲线同样上升；本模型预测的温区平均系数低于NIST参考，偏差为{m['relative_error_percent']:+.2f}%。因此可认为模型再现了趋势和量级，同时保留了可分辨的系统偏差，不能说与实验完全一致。")
        md(f"### 3.2 平衡表现与参数敏感性\n\n24条轨迹的采样平均温度相对目标的最大绝对偏差为{runs.temperature_offset_K.abs().max():.3f} K；平均压力范围为{runs.P_mean_bar.min():.2f}至{runs.P_mean_bar.max():.2f} bar。NPT的零外压指控制目标，不要求每帧压力为零。晶格常数前后半段差异除以其分块标准误的最大值为{runs.halves_z.max():.3f}，未触发预设强漂移复核阈值4；该诊断不能单独证明采样充分。\n\n尺寸由4000增至6912原子、步长由1降至0.5 fs、时长延至1 ns平衡+1 ns采样时，平均系数相对基准r1的变化分别为{checks.iloc[0].change_percent:+.4f}%、{checks.iloc[1].change_percent:+.4f}%和{checks.iloc[2].change_percent:+.4f}%。这些变化{'均小于预设±2%目标' if m['all_practical_checks_pass'] else '未全部满足预设目标'}，范围仅限已检查的端点和参数组合。")
        code("metrics=json.loads(Path('results/metrics.json').read_text(encoding='utf8'))\ndisplay(pd.read_csv('results/temperature_summary.csv'))\ndisplay(pd.read_csv('results/alpha_replicates.csv'))\ndisplay(pd.read_csv('results/validation_comparison.csv'))\nprint(json.dumps(metrics,ensure_ascii=False,indent=2))")
        for name,caption in [('lattice_temperature','图1 晶格常数及3次重复统计区间'),('expansion_comparison','图2 统一300 K参考后的线膨胀与体积变化'),('alpha_comparison','图3 同定义温区平均系数对照'),('equilibration_300','图4 300 K温度、压力、能量、体积轨迹'),('equilibration_800','图5 800 K对应轨迹'),('block_sensitivity','图6 分块长度敏感性'),('validation','图7 尺寸、步长、时长检查')]:
            md(caption);code(f"display(Image(filename='results/figures/{name}.png'))")
        blocks=pd.read_csv(ROOT/'results/block_statistics.csv').pivot(index='job',columns='block_ps',values='a_sem_A')
        ratios=blocks[50]/blocks[20]
        md(f'## 4. 讨论\n\n20 ps增至50 ps分块时，均值标准误比值的中位数为{ratios.median():.3f}，范围为{ratios.min():.3f}–{ratios.max():.3f}，没有所有轨迹一致增长的现象。块数减少也会使误差估计自身波动，因此不能把单条曲线的局部平台看作严格证明。最终统计区间主要来自3次独立重复，自由度仅为2。\n\n端点平均系数的三项参数敏感性远小于当前约3.85%的实验偏差，表明在已检查范围内，单纯改变这些数值参数不足以消除实验差异。这支持进一步关注势函数适用性及物理模型差别，但不构成对误差来源的唯一归因。\n\n尺寸、步长和时长检查只针对300与800 K，中间温度没有完整扫描。Mishin EAM1的拟合目标、经典原子统计、理想晶体与真实试样的区别均可能产生偏差；本研究没有计算各因素的单独贡献。统计误差棒没有包含实验不确定度、势函数系统误差和经典近似误差。\n\n本项目没有自行运行DFT、机器学习势或第二种EAM势。文献横向比较是方法与既有证据的比较，不能将其写成我们统一协议下跑出的势函数精度排行榜。')
        if not m['all_practical_checks_pass']:
            md('**至少一项预设验证未通过：当前结果应视为待进一步验证，不能写成所有参数均已收敛。应依据validation_comparison.csv延长或扩大相关计算。**')
        md(f"## 5. 结论\n\n当前协议得到300–800 K平均线膨胀系数{m['alpha_300_800_Kinv']*1e6:.4f}×10⁻⁶ K⁻¹，相对NIST参考偏差{m['relative_error_percent']:+.2f}%。该结果针对Mishin EAM1、近零压、周期理想单晶及本项目温区，不推广到晶界、自由表面、熔化或其他势。可复现实验文件和全部统计表随仓库提供。")
    else:
        md('## 4. 待数据完成的讨论与结论\n\n尚无完整正式数据，不给出最终系数、误差排名或“已收敛”判断。完成后按实际数据生成晶格常数曲线、温压稳定性、分块误差、三项参数验证及实验偏差。')
    md('## 参考文献与复现\n\n完整参考文献和横向方法比较见[文献调研](docs/literature_review.md)，势文件说明见[来源](potentials/SOURCE.md)。运行步骤、目录树及未完成事项见[README](README.md)。Notebook只读取结果和展示图，不在打开时自动启动数小时的模拟。')
    if complete:
        md('### 计算中断与数据核查说明\n\n两条0.5 fs验证在中断后从200 ps平衡重启点重新连续采满500 ps，中断的采样后缀没有混入结果。重启边界自动输出的200 ps行与原前缀重复，严格检查确认只有这一处重复且值一致后移除，原始文件与修复前metadata均保留。该点不在大于200 ps的正式采样区间，清理不改变生产样本。检查同时核对最终重启点为1400000步、700 ps、4000原子。详见[恢复记录](docs/recovery_record.md)。')
    book=nb.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Copper MD (Python)','language':'python','name':'copper-md'}})
    nb.validate(book);nb.write(book,ROOT/'铜热膨胀研究.ipynb')
    print('Paper generated; complete data:',complete)

if __name__=='__main__':main()
