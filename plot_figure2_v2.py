"""Reproducible composite figure from 18 completed copper MD trajectories."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.mlab import GaussianKDE
from matplotlib.ticker import FormatStrFormatter, ScalarFormatter
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
OUT = ROOT/'results/paper_figures'
STEM = 'figure2_composite_v2'
TEMPS = [300,400,500,600,700,800]
COLORS = ['#477BA8','#4C91B0','#55A5B3','#62B5AC','#79BEA5','#95C7A4']
TCRIT = 4.3026527299


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary = ROOT/'results/run_summary.csv'
    original_notebook_hash = digest(ROOT/'铜热膨胀研究.ipynb')
    summary_hash = digest(summary)
    prod = pd.read_csv(summary).query("group == 'production'")
    assert len(prod)==18 and not prod.duplicated(['temperature','replicate']).any()
    assert prod.groupby('temperature').size().eq(3).all()
    assert sorted(prod.temperature.unique())==TEMPS
    assert prod.atoms.eq(4000).all() and prod.samples.eq(5000).all()
    fluctuations = {t:[] for t in TEMPS}
    provenance, centered_rows = [], []
    columns = ['step','time_ps','temperature_K','pressure_bar','potential_eV',
               'lx_A','ly_A','lz_A','volume_A3']
    for row in prod.itertuples(index=False):
        folder = ROOT/'outputs/production'/row.id
        assert json.loads((folder/'metadata.json').read_text())['status']=='completed'
        path = folder/'thermo.dat'
        frame = pd.read_csv(path,sep=r'\s+',comment='#',names=columns)
        sample = frame.loc[(frame.time_ps>row.equil_ps)&(frame.time_ps<=row.equil_ps+row.sample_ps)]
        assert len(sample)==5000 and row.equil_ps==200 and row.sample_ps==500
        assert np.isfinite(sample.to_numpy()).all()
        values = (sample.lx_A/10).to_numpy()
        assert np.isclose(values.mean(),row.a_A,rtol=1e-12)
        assert np.isclose(sample.volume_A3.mean(),row.volume_A3,rtol=1e-12)
        delta = values-values.mean()
        assert abs(delta.mean())<1e-12
        fluctuations[row.temperature].append(delta)
        centered_rows.append({'id':row.id,'temperature_K':row.temperature,'replicate':row.replicate,
                              'raw_samples':len(values),'centering_lattice_mean_A':values.mean()})
        provenance.append({'path':str(path.relative_to(ROOT)),'sha256':digest(path)})
    a = prod.pivot(index='temperature',columns='replicate',values='a_A').sort_index()
    volume = prod.pivot(index='temperature',columns='replicate',values='volume_A3').sort_index()
    relative_volume = 100*(volume.div(volume.loc[300],axis=1)-1)
    mean = a.mean(axis=1)
    ci = TCRIT*a.std(axis=1,ddof=1)/np.sqrt(3)
    linear_pct = float((100*(a.loc[800]/a.loc[300]-1)).mean())
    volume_pct = float(relative_volume.loc[800].mean())
    assert np.isclose(mean.loc[300],3.63193,rtol=0,atol=5e-6)
    assert np.isclose(mean.loc[800],3.66441,rtol=0,atol=5e-6)
    assert (relative_volume.loc[300]==0).all()
    available = {font.name for font in font_manager.fontManager.ttflist}
    font = 'Times New Roman' if 'Times New Roman' in available else 'STIXGeneral'
    plt.rcParams.update({
        'font.family':font,'font.size':8.5,'axes.labelsize':9.5,
        'xtick.labelsize':8,'ytick.labelsize':8.5,
        'text.color':'#303840','axes.labelcolor':'#303840','axes.edgecolor':'#303840',
        'xtick.color':'#303840','ytick.color':'#303840','axes.linewidth':.65,
        'axes.spines.top':False,'axes.spines.right':False,
        'figure.facecolor':'white','axes.facecolor':'white','mathtext.fontset':'stix',
        'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':True,
    })
    fig = plt.figure(figsize=(185/25.4,132/25.4))
    axa = fig.add_axes([.103,.14,.390,.770])
    axb = fig.add_axes([.650,.635,.330,.275])
    axc = fig.add_axes([.650,.235,.330,.220])
    cax = fig.add_axes([.650,.100,.330,.016])
    for axis in [axa,axb]:
        axis.set_axisbelow(True)
        axis.grid(axis='y',color='#E8EDF0',lw=.55)
        axis.tick_params(direction='out',length=3,width=.6,pad=4)

    # Panel a: independent trajectory means and an ensemble confidence interval.
    axa.errorbar(TEMPS,mean,yerr=ci,fmt='-',color='#365F78',lw=1.5,
                 elinewidth=.8,capsize=2.3,capthick=.8,zorder=3)
    axa.scatter(TEMPS,mean,c=COLORS,s=25,edgecolors='white',linewidths=.55,zorder=5)
    offsets = [-12,0,12]
    markers = ['s','D','^']
    for rep,offset,marker in zip([1,2,3],offsets,markers):
        axa.scatter(np.array(TEMPS)+offset,a[rep],s=18,marker=marker,
                    facecolors='none',edgecolors='#62727D',linewidths=.65,zorder=6)
    axa.set(xlim=(273,827),ylim=(3.629,3.6705),xlabel='Temperature / K',ylabel='Lattice parameter / Å')
    axa.set_xticks(TEMPS)
    axa.set_yticks([3.63,3.64,3.65,3.66,3.67])
    axa.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    axa.xaxis.labelpad=8
    axa.yaxis.labelpad=8
    axa.text(.065,.88,'300–800 K linear expansion',transform=axa.transAxes,
             fontsize=8.5,color='#B77D55',va='top')
    axa.text(.065,.83,f'{linear_pct:.4f}%',transform=axa.transAxes,
             fontsize=13,color='#B77D55',va='top')
    handles = [Line2D([],[],marker=m,color='none',markeredgecolor='#62727D',
                     markerfacecolor='none',markersize=4,label=f'r{i}')
               for i,m in enumerate(markers,1)]
    axa.legend(handles=handles,frameon=False,ncol=3,loc='lower right',
               fontsize=8,handletextpad=.25,columnspacing=.8,borderaxespad=.5)

    # Panel b: descriptive KDEs of correlated, within-trajectory-centered samples.
    distribution_rows, kde_settings = [], []
    for temp,color in zip(TEMPS,COLORS):
        values = np.concatenate(fluctuations[temp])
        assert len(values)==15000
        kde = GaussianKDE(values,bw_method='scott')
        support = np.linspace(values.min(),values.max(),256)
        density = kde(support)
        anchor = temp-7
        axb.fill_betweenx(support,anchor-39*density/density.max(),anchor,
                          facecolor=color,edgecolor=color,lw=.6,alpha=.58,zorder=3)
        box = axb.boxplot([values],positions=[temp+13],widths=15,whis=1.5,
            showfliers=False,patch_artist=True,manage_ticks=False,
            boxprops={'facecolor':color,'edgecolor':'#50616A','linewidth':.65,'alpha':.65},
            medianprops={'color':'#303840','linewidth':.85},
            whiskerprops={'color':'#50616A','linewidth':.65},
            capprops={'color':'#50616A','linewidth':.65})
        quantiles = np.quantile(values,[.25,.5,.75])
        for q,value in zip(['q25','median','q75'],quantiles):
            distribution_rows.append({'temperature_K':temp,'statistic':q,'delta_a_A':value})
        for q,value in [('min',values.min()),('max',values.max()),('sd_descriptive',values.std(ddof=1))]:
            distribution_rows.append({'temperature_K':temp,'statistic':q,'delta_a_A':value})
        kde_settings.append({'temperature_K':temp,'pooled_correlated_samples':len(values),
            'bandwidth_rule':'Scott','factor':float(kde.factor),
            'bandwidth_A':float(np.sqrt(kde.covariance[0,0])),
            'support_min_A':float(values.min()),'support_max_A':float(values.max()),'support_points':256})
    axb.axhline(0,color='#7E8B93',lw=.65,zorder=1)
    lim = max(abs(np.concatenate(fluctuations[t])).max() for t in TEMPS)*1.08
    axb.set(xlim=(245,850),ylim=(-lim,lim),xlabel='Temperature / K',ylabel='Lattice fluctuation, δa / Å')
    axb.set_xticks(TEMPS)
    axb.tick_params(axis='y',labelsize=8)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-3,3))
    axb.yaxis.set_major_formatter(formatter)
    axb.yaxis.get_offset_text().set_fontsize(8)
    axb.xaxis.labelpad=6
    axb.yaxis.labelpad=5

    # Panel c: one shared absolute color scale, no row-wise rescaling.
    cmap = LinearSegmentedColormap.from_list('volume_teal',['#F4F9F8','#D2EAE4','#92C7BB','#3F8F83'])
    matrix = relative_volume.T.to_numpy()
    vmax = 3.0
    image = axc.imshow(matrix,aspect='auto',cmap=cmap,vmin=0,vmax=vmax,interpolation='nearest')
    axc.set_xticks(range(6),labels=TEMPS)
    axc.set_yticks(range(3),labels=['Repeat 1','Repeat 2','Repeat 3'])
    axc.set_xlabel('Temperature / K',labelpad=7)
    axc.set_ylabel('Independent repeat',labelpad=6)
    axc.set_xticks(np.arange(-.5,6,1),minor=True)
    axc.set_yticks(np.arange(-.5,3,1),minor=True)
    axc.grid(which='minor',color='white',lw=1.2)
    axc.tick_params(which='minor',bottom=False,left=False)
    axc.tick_params(which='major',length=0,pad=5)
    for spine in axc.spines.values():
        spine.set_visible(False)
    for r in range(3):
        for c in range(6):
            axc.text(c,r,f'{matrix[r,c]:.3f}',ha='center',va='center',fontsize=8,
                      color='white' if matrix[r,c]>2.4 else '#303840')
    colorbar = fig.colorbar(image,cax=cax,orientation='horizontal',ticks=[0,1,2,3])
    colorbar.set_label('Volume change / %',fontsize=9,labelpad=4)
    colorbar.ax.tick_params(labelsize=8,length=2,width=.5,pad=2)
    colorbar.outline.set_visible(False)
    # Align panel tags in figure coordinates to avoid collisions with exponents.
    fig.text(.025,.940,'a',fontsize=12,fontweight='bold')
    fig.text(.565,.940,'b',fontsize=12,fontweight='bold')
    fig.text(.565,.485,'c',fontsize=12,fontweight='bold')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    essential_text = list(fig.texts)
    for axis in [axa,axb,axc,cax]:
        essential_text += [axis.xaxis.label,axis.yaxis.label,*axis.texts,
                           *axis.get_xticklabels(),*axis.get_yticklabels(),axis.yaxis.get_offset_text()]
    for text in essential_text:
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        assert box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.width and box.y1<=fig.bbox.height, text.get_text()
    assert (mean-ci).min()>axa.get_ylim()[0] and (mean+ci).max()<axa.get_ylim()[1]
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f'{STEM}.pdf')
    fig.savefig(OUT/f'{STEM}.png',dpi=600)
    # A compact proof approximates intended print size for visual inspection.
    fig.savefig(OUT/f'{STEM}_preview.png',dpi=150)
    plt.close(fig)
    summary_table = pd.DataFrame({'temperature_K':TEMPS,'lattice_mean_A':mean,
                                   'lattice_ci95_halfwidth_A':ci})
    for rep in [1,2,3]:
        summary_table[f'lattice_r{rep}_A']=a[rep]
        summary_table[f'volume_change_r{rep}_percent']=relative_volume[rep]
    summary_table.to_csv(OUT/f'{STEM}_data.csv',index=False)
    pd.DataFrame(distribution_rows).to_csv(OUT/f'{STEM}_distribution.csv',index=False)
    pd.DataFrame(centered_rows).to_csv(OUT/f'{STEM}_centering.csv',index=False)
    caption = '''**Figure 2. Mean thermal expansion, lattice fluctuations, and repeat-resolved volume changes of periodic single-crystal copper.** (a) Mean lattice parameter from three independently initialized production trajectories at each temperature. Error bars are two-sided 95% Student-t confidence intervals with two degrees of freedom and may be smaller than the markers. Open squares, diamonds, and triangles show repeats r1, r2, and r3, respectively; fixed horizontal offsets of −12, 0, and +12 K are used only for visibility, without changing their lattice-parameter values. Filled circles are plotted at the actual temperatures, and connecting lines are guides to the eye. (b) Half-violin densities and boxplots summarize production-sample fluctuations after subtracting each trajectory's own mean lattice parameter and pooling the three centered trajectories per temperature. Densities use Scott's bandwidth rule and are displayed only over the observed ranges. Violin widths are scaled to a common maximum and do not encode sample size. Boxes show the interquartile range and median; whiskers reach the most extreme observations within 1.5 interquartile ranges. Boxplot outlier symbols are omitted; the density display retains the full observed range. These correlated samples describe fluctuations, not independent replicate measurements or confidence intervals. (c) Directly averaged volumes normalized to each repeat's own 300 K volume, with values in percent. The 300 K zeros follow from normalization. Blue-to-green marker and distribution colors identify temperature in (a,b), whereas the continuous teal scale in (c) encodes volume change. All panels exclude the first 200 ps and use the subsequent 500 ps production segment.\n'''
    (ROOT/'docs/figure2_caption_v2.md').write_text(caption,encoding='utf8')
    assert digest(summary)==summary_hash
    assert digest(ROOT/'铜热膨胀研究.ipynb')==original_notebook_hash
    for record in provenance:
        assert digest(ROOT/record['path'])==record['sha256']
    report = {'summary_source':'results/run_summary.csv','summary_sha256':summary_hash,
        'production_sources':provenance,'production_trajectories':18,
        'production_window_ps':'200 < time <= 700','samples_per_trajectory':5000,
        'confidence_interval':'Mean ± t(0.975, df=2) * sample SD / sqrt(3)',
        't_critical':TCRIT,'repeat_x_offsets_K':offsets,
        'fluctuations':'a(t) - trajectory-specific production mean, pooled descriptively by temperature',
        'density_method':'matplotlib.mlab.GaussianKDE; Scott bandwidth; 256 points within observed min/max; equal maximum displayed half-width',
        'density_settings':kde_settings,'boxplot':'median, IQR; 1.5-IQR whiskers; outlier symbols omitted',
        'volume_normalization':'100 * (mean V_r(T) / mean V_r(300) - 1)',
        'heatmap_color_limits_percent':[0,vmax],'heatmap_annotations_decimals':3,
        'font':font,'font_file':font_manager.findfont(font),'size_mm':[185,132],
        'png_dpi':600,'preview_dpi':150,
        'endpoints':{'lattice_300_A':float(mean.loc[300]),'lattice_800_A':float(mean.loc[800]),
                     'linear_increase_percent':linear_pct,'volume_increase_percent':volume_pct},
        'notes':[] if font=='Times New Roman' else ['Times New Roman unavailable; STIXGeneral used.'],
        'verification':['18 means checked against raw data','All source SHA256 hashes unchanged',
                        'Notebook SHA256 unchanged','Text extents inside figure canvas','Confidence bars inside axes limits']}
    if f'{linear_pct:.4f}'!='0.8944' or f'{volume_pct:.4f}'!='2.7074':
        report['notes'].append('Calculated endpoint percentages differ from expected rounded checks; actual values retained.')
    (OUT/f'{STEM}_audit.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps({'font':font,'endpoints':report['endpoints'],'notes':report['notes'],
                      'figure':str(OUT/f'{STEM}.png')},indent=2))


if __name__=='__main__':
    main()
