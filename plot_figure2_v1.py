"""Plot production MD measurements; never launch simulations or alter raw data."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FormatStrFormatter, FixedLocator

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results/paper_figures'
STEM = 'figure2_lattice_volume_v1'
TCRIT = 4.3026527299  # Two-sided 95% Student-t quantile, df = 2.


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary = ROOT/'results/run_summary.csv'
    before = sha256(summary)
    runs = pd.read_csv(summary)
    prod = runs.loc[runs.group.eq('production')].copy()
    assert len(prod) == 18
    assert not prod.duplicated(['temperature', 'replicate']).any()
    assert sorted(prod.temperature.unique()) == [300, 400, 500, 600, 700, 800]
    assert prod.groupby('temperature').size().eq(3).all()
    assert prod.samples.eq(5000).all()
    raw_records = []
    columns = ['step', 'time_ps', 'temperature_K', 'pressure_bar', 'potential_eV',
               'lx_A', 'ly_A', 'lz_A', 'volume_A3']
    # Independently check plotted means against the recorded production windows.
    for row in prod.itertuples(index=False):
        folder = ROOT/'outputs/production'/row.id
        meta = json.loads((folder/'metadata.json').read_text(encoding='utf8'))
        assert meta['status'] == 'completed'
        raw = folder/'thermo.dat'
        frame = pd.read_csv(raw, sep=r'\s+', comment='#', names=columns)
        sample = frame.loc[frame.time_ps > row.equil_ps]
        assert len(sample) == row.samples
        assert np.isfinite(sample.to_numpy()).all()
        assert np.isclose((sample.lx_A/row.cells).mean(), row.a_A, rtol=1e-12)
        assert np.isclose(sample.volume_A3.mean(), row.volume_A3, rtol=1e-12)
        raw_records.append({'file': str(raw.relative_to(ROOT)), 'sha256': sha256(raw)})

    a = prod.pivot(index='temperature', columns='replicate', values='a_A').sort_index()
    volume = prod.pivot(index='temperature', columns='replicate', values='volume_A3').sort_index()
    relative_volume = 100*(volume.div(volume.loc[300], axis=1)-1)
    linear_pct = float((100*(a.loc[800]/a.loc[300]-1)).mean())
    volume_pct = float(relative_volume.loc[800].mean())
    result = pd.DataFrame({
        'temperature_K': a.index,
        'lattice_mean_A': a.mean(axis=1).to_numpy(),
        'lattice_ci95_halfwidth_A': (TCRIT*a.std(axis=1,ddof=1)/np.sqrt(3)).to_numpy(),
        'relative_volume_mean_percent': relative_volume.mean(axis=1).to_numpy(),
        'relative_volume_ci95_halfwidth_percent': (TCRIT*relative_volume.std(axis=1,ddof=1)/np.sqrt(3)).to_numpy(),
    })
    assert np.isclose(result.lattice_mean_A.iloc[0],3.63193,atol=5e-6,rtol=0)
    assert np.isclose(result.lattice_mean_A.iloc[-1],3.66441,atol=5e-6,rtol=0)
    assert result.relative_volume_mean_percent.iloc[0] == 0
    assert result.relative_volume_ci95_halfwidth_percent.iloc[0] == 0
    notes = []
    if f'{volume_pct:.4f}' != '2.7074':
        notes.append(f'Repeat-normalized volume increase is {volume_pct:.4f}%, rather than 2.7074%. The calculated value is plotted.')
    font_names = {f.name for f in font_manager.fontManager.ttflist}
    family = 'Times New Roman' if 'Times New Roman' in font_names else 'STIXGeneral'
    if family != 'Times New Roman':
        notes.append('Times New Roman unavailable; STIXGeneral was substituted consistently.')
    plt.rcParams.update({
        'font.family': family, 'font.size':9,
        'axes.labelsize':10, 'axes.labelcolor':'#303840',
        'axes.edgecolor':'#303840', 'axes.linewidth':.7,
        'xtick.labelsize':9, 'ytick.labelsize':9,
        'xtick.color':'#303840', 'ytick.color':'#303840',
        'text.color':'#303840', 'mathtext.fontset':'stix',
        'axes.spines.top':False, 'axes.spines.right':False,
        'axes.unicode_minus':True, 'figure.facecolor':'white',
        'axes.facecolor':'white', 'pdf.fonttype':42, 'ps.fonttype':42,
        'svg.fonttype':'none',
    })
    fig, axes = plt.subplots(1,2,figsize=(180/25.4,78/25.4),sharex=True)
    fig.subplots_adjust(left=.095,right=.983,bottom=.205,top=.91,wspace=.30)
    colors = ['#438EAC','#62A99A']
    series = [('lattice_mean_A','lattice_ci95_halfwidth_A','Lattice parameter / Å'),
              ('relative_volume_mean_percent','relative_volume_ci95_halfwidth_percent','Relative volume change / %')]
    artists = []
    for axis,color,(mean,ci,label),tag in zip(axes,colors,series,['a','b']):
        axis.set_axisbelow(True)
        axis.grid(axis='y',color='#E8EDF0',lw=.6)
        axis.errorbar(result.temperature_K,result[mean],yerr=result[ci],
            fmt='o-',color=color,lw=1.5,ms=4.5,mec='white',mew=.55,
            ecolor=color,elinewidth=.8,capsize=2.2,capthick=.8,zorder=3)
        axis.set(xlim=(280,820),xlabel='Temperature / K',ylabel=label)
        axis.xaxis.set_major_locator(FixedLocator([300,400,500,600,700,800]))
        axis.tick_params(direction='out',length=3,width=.65,pad=4)
        axis.set_xlabel('Temperature / K',labelpad=7)
        axis.set_ylabel(label,labelpad=7)
        artists.append(axis.text(-.17,1.045,tag,transform=axis.transAxes,
                                 fontsize=12,fontweight='bold',va='bottom',ha='left'))
    axes[0].set_ylim(3.629,3.670)
    axes[0].set_yticks([3.63,3.64,3.65,3.66,3.67])
    axes[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axes[1].set_ylim(-.14,3.12)
    axes[1].set_yticks([0,1,2,3])
    axes[1].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    artists.append(axes[0].text(.065,.83,f'Linear expansion: {linear_pct:.4f}%',
        transform=axes[0].transAxes,color='#B77D55',fontsize=8.5,va='top'))
    artists.append(axes[1].text(.065,.83,f'Volume increase: {volume_pct:.4f}%',
        transform=axes[1].transAxes,color='#B77D55',fontsize=8.5,va='top'))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Check all essential text extents against the full export canvas.
    for artist in artists+[obj for axis in axes for obj in
            [axis.xaxis.label,axis.yaxis.label,*axis.get_xticklabels(),*axis.get_yticklabels()]]:
        box = artist.get_window_extent(renderer)
        assert box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.width and box.y1<=fig.bbox.height, (artist.get_text(),box)
    # Check uncertainty bars remain within the chosen data limits.
    for axis,(mean,ci,_) in zip(axes,series):
        lo,hi = axis.get_ylim()
        assert (result[mean]-result[ci]).min() > lo
        assert (result[mean]+result[ci]).max() < hi
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f'{STEM}.pdf',facecolor='white')
    fig.savefig(OUT/f'{STEM}.png',dpi=600,facecolor='white')
    plt.close(fig)
    result.to_csv(OUT/f'{STEM}_data.csv',index=False)
    caption = (
        '**Figure 2. Thermal expansion of periodic single-crystal copper from 300 to 800 K.** '
        '(a) Mean lattice parameter and (b) relative volume change with respect to 300 K. '
        'Each point summarizes three independently initialized production trajectories. '
        'Error bars represent two-sided 95% Student-t confidence intervals with two degrees of freedom '
        'and may be smaller than the markers. Relative volume changes were calculated from the '
        'directly averaged cell volumes and normalized separately within each repeat before aggregation. '
        'The zero value and zero interval at 300 K follow from this normalization. '
        'Connecting lines are guides to the eye, not fitted models.\n')
    (ROOT/'docs/figure2_caption_v1.md').write_text(caption,encoding='utf8')
    report = {
        'input_summary':str(summary.relative_to(ROOT)), 'summary_sha256':before,
        'raw_production_files':raw_records,'production_trajectories':18,'repeats_per_temperature':3,
        'confidence_level':.95,'student_t_df':2,'student_t_critical':TCRIT,
        'font':family,'font_file':font_manager.findfont(family),
        'size_mm':[180,78],'png_dpi':600,
        'lattice_300_A':float(a.loc[300].mean()),'lattice_800_A':float(a.loc[800].mean()),
        'mean_repeat_linear_increase_percent':linear_pct,
        'mean_repeat_volume_increase_percent':volume_pct,
        'ratio_of_pooled_mean_volume_increase_percent':float(100*(volume.loc[800].mean()/volume.loc[300].mean()-1)),
        'notes':notes,'caption_file':'docs/figure2_caption_v1.md',
        'checks':['18 production means agree with their raw production windows',
                  'Text extents inside export canvas','Confidence bars inside axes limits',
                  'No smoothing, interpolation or modified source data'],
    }
    assert sha256(summary)==before
    (OUT/f'{STEM}_audit.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps({'font':family,'lattice_300':report['lattice_300_A'],
        'lattice_800':report['lattice_800_A'],'linear_percent':linear_pct,
        'volume_percent':volume_pct,'notes':notes,'output':str(OUT/f'{STEM}.png')},indent=2))


if __name__ == '__main__':
    main()
