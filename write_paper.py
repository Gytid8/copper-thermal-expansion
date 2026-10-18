"""Build the concise manuscript from reviewed prose and reproducible code cells."""
from pathlib import Path
import re
import nbformat as nb
ROOT=Path(__file__).resolve().parent

CODE={
'setup':'''from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from IPython.display import display, Image
from run_study import manifest, folder
from analyze_study import read_job, reference
from md_project import plotting
jobs = manifest()
checked = [read_job(job)[0] for job in jobs]
assert len(checked) == 24
for job in jobs:
    meta = json.loads((folder(job)/"metadata.json").read_text(encoding="utf8"))
    text = (folder(job)/"input.in").read_text(encoding="ascii")
    assert hashlib.sha256(text.encode("ascii")).hexdigest() == meta["input_sha256"]
figdir = Path("results/paper_figures"); figdir.mkdir(exist_ok=True)
plt = plotting()
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})''',
'protocol':'''# 展示实际使用的建模与积分设置，不在论文中重跑MD。
for file in ["outputs/model/build.in", "outputs/production/prod_T300_r1/input.in"]:
    text = Path(file).read_text(encoding="ascii")
    keys = ("lattice ", "region ", "boundary ", "pair_style ", "pair_coeff ", "timestep ", "fix integrator ", "run ")
    print("\\n" + file)
    print("\\n".join(line for line in text.splitlines() if line.startswith(keys)))''',
'statistics':'''runs = pd.DataFrame(checked)
prod = runs[runs.group == "production"]
a = prod.pivot(index="temperature", columns="replicate", values="a_A")
v = prod.pivot(index="temperature", columns="replicate", values="volume_A3")
eps = a.div(a.loc[300], axis=1) - 1
alpha = eps.loc[800] / 500
tcrit = 4.3026527299
estimate = alpha.mean()
halfwidth = tcrit * alpha.std(ddof=1) / np.sqrt(3)
ref, refalpha = reference()
bias = 100 * (estimate / refalpha - 1)
checks = pd.read_csv("results/validation_comparison.csv")
assert np.isclose(estimate * 1e6, 17.8882558243, rtol=1e-8)
print(f"平均系数 = {estimate*1e6:.3f} ± {halfwidth*1e6:.3f} ×10⁻⁶ K⁻¹（95%统计区间）")
print(f"NIST = {refalpha*1e6:.4f} ×10⁻⁶ K⁻¹；相对偏差 = {bias:.2f}%")''',
'figure1':'''fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
for panel, data, scale, ylabel in [(ax[0], a, 1, "晶格常数 / Å"), (ax[1], v, 1e-3, r"体积 / $10^3$ Å$^3$")]:
    panel.errorbar(data.index, data.mean(axis=1)*scale, yerr=tcrit*data.std(axis=1)/np.sqrt(3)*scale,
                   fmt="o-", color="#a85a2a", capsize=3, ms=4)
    panel.set(xlabel="温度 / K", ylabel=ylabel)
    panel.grid(alpha=.18)
for panel, tag in zip(ax, ["(a)", "(b)"]): panel.set_title(tag, loc="left", fontsize=11)
fig.tight_layout()
fig.savefig(figdir/"figure1_lattice_volume.png", dpi=220)
display(Image(filename=str(figdir/"figure1_lattice_volume.png"))); plt.close(fig)''',
'figure2':'''mean_eps = eps.mean(axis=1)
ci_eps = tcrit * eps.std(axis=1) / np.sqrt(3)
reference_eps = ref.set_index("temperature_K").strain_from_300.reindex(eps.index)
fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
ax[0].errorbar(eps.index, mean_eps*100, yerr=ci_eps*100, fmt="o-", color="#a85a2a", ms=4, label="MD：三次重复")
ax[0].plot(eps.index, reference_eps*100, "s--", color="#226b8e", ms=4, label="NIST：1990修订")
ax[0].set(xlabel="温度 / K", ylabel="相对300 K的线膨胀 / %")
ax[0].legend(fontsize=8)
ax[1].axhline(0, color="gray", lw=.8)
ax[1].errorbar(eps.index, (mean_eps-reference_eps)*1e6, yerr=ci_eps*1e6,
               fmt="o-", color="#a85a2a", capsize=3, ms=4)
ax[1].set(xlabel="温度 / K", ylabel="MD减NIST的线膨胀差 / ppm")
for panel, tag in zip(ax, ["(a)", "(b)"]):
    panel.set_title(tag, loc="left", fontsize=11); panel.grid(alpha=.18)
fig.tight_layout()
fig.savefig(figdir/"figure2_reference_residual.png", dpi=220)
display(Image(filename=str(figdir/"figure2_reference_residual.png"))); plt.close(fig)''',
'result_table':'''table = pd.DataFrame({"项目": ["MD：三次重复", "NIST SRM 736修订表"],
    "平均系数": [f"{estimate*1e6:.3f}", f"{refalpha*1e6:.4f}"],
    "95%统计区间半宽": [f"{halfwidth*1e6:.3f}", "未作联合传播"],
    "相对NIST偏差": [f"{bias:.2f}%", "参考基准"]})
display(table.style.hide(axis="index"))''',
'figure3':'''fig, ax = plt.subplots(figsize=(7.2, 3.1))
bars = ax.bar(["尺寸：4000→6912原子", "步长：1→0.5 fs", "时长：延至1+1 ns"], checks.change_percent,
              color=["#226b8e", "#a85a2a", "#557764"], width=.48)
ax.axhline(0, color="gray", lw=.8)
ax.bar_label(bars, labels=[f"{x:+.4f}%" for x in checks.change_percent], padding=5, fontsize=9)
ax.set(ylabel="平均系数相对基准r1的变化 / %", ylim=(-.027, .09))
ax.set_title("预设目标：绝对变化不超过2%（纵轴放大）", fontsize=10)
ax.grid(axis="y", alpha=.18)
fig.tight_layout(); fig.savefig(figdir/"figure3_sensitivity.png", dpi=220)
display(Image(filename=str(figdir/"figure3_sensitivity.png"))); plt.close(fig)''',
'supplement1':'''fig, axes = plt.subplots(4, 2, figsize=(9, 8), sharex=True)
keys = ["temperature_K", "pressure_bar", "potential_eV", "volume_A3"]
labels = ["温度 / K", "压力 / GPa", "势能 / (eV/atom)", r"体积 / $10^3$ Å$^3$"]
for col, temp in enumerate([300, 800]):
    job = next(j for j in jobs if j["id"] == f"prod_T{temp}_r1")
    frame = pd.read_csv(folder(job)/"thermo.dat", sep=r"\\s+", comment="#", names=
        ["step", "time_ps", "temperature_K", "pressure_bar", "potential_eV", "lx_A", "ly_A", "lz_A", "volume_A3"])
    for row, (key, factor) in enumerate(zip(keys, [1, 1e-4, 1/4000, 1e-3])):
        panel = axes[row, col]; y = frame[key] * factor
        panel.plot(frame.time_ps, y, lw=.35, color="#c0cdd4")
        panel.plot(frame.time_ps, y.rolling(50, min_periods=1).mean(), lw=.75, color="#226b8e")
        panel.axvline(200, color="black", ls="--", lw=.8)
        panel.set_ylabel(labels[row], fontsize=9)
    axes[0, col].set_title(f"{temp} K，正式重复r1", fontsize=10)
    axes[-1, col].set_xlabel("时间 / ps")
fig.tight_layout(); fig.savefig(figdir/"figureS1_endpoint_traces.png", dpi=180)
display(Image(filename=str(figdir/"figureS1_endpoint_traces.png"))); plt.close(fig)''',
'supplement2':'''blocks = pd.read_csv("results/block_statistics.csv")
fig, axes = plt.subplots(1, 2, figsize=(9, 3.2))
for panel, temp in zip(axes, [300, 800]):
    block = blocks[blocks.job == f"prod_T{temp}_r1"]
    panel.plot(block.block_ps, block.a_sem_A*1e6, "o-", ms=4, color="#226b8e")
    panel.set(xlabel="分块长度 / ps", ylabel=r"均值标准误 / $10^{-6}$ Å", title=f"{temp} K，正式重复r1")
    panel.grid(alpha=.18)
fig.tight_layout(); fig.savefig(figdir/"figureS2_block_statistics.png", dpi=200)
display(Image(filename=str(figdir/"figureS2_block_statistics.png"))); plt.close(fig)'''
}

def main():
    source=(ROOT/'docs/manuscript.md').read_text(encoding='utf8')
    parts=re.split(r'<!-- code:(\w+) -->',source)
    cells=[]
    for i,part in enumerate(parts):
        if i%2:
            assert part in CODE,part
            compile(CODE[part],'<notebook:'+part+'>','exec')
            cells.append(nb.v4.new_code_cell(CODE[part]))
        elif part.strip():cells.append(nb.v4.new_markdown_cell(part.strip()))
    book=nb.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Copper MD (Python)','language':'python','name':'copper-md'},
        'language_info':{'name':'python'},'manuscript_revision':'2026-10-18 concise scientific rewrite'})
    nb.validate(book)
    nb.write(book,ROOT/'铜热膨胀研究.ipynb')
    nlines=sum(len(code.splitlines()) for code in CODE.values())
    print('Reviewed manuscript built; code lines:',nlines,'code cells:',len(CODE))

if __name__=='__main__':main()
