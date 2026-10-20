"""Recompute correlated-data block statistics and plot their duration sensitivity."""
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
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/paper_figures'
STEM='figure6_block_sensitivity_v1'
BLOCKS=[1,2,5,10,20,50]
TEMPS=[300,400,500,600,700,800]
COLORS=['#477BA8','#4C91B0','#55A5B3','#62B5AC','#79BEA5','#95C7A4']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    protected_paths=[ROOT/'results/block_statistics.csv',ROOT/'results/run_summary.csv',
                     ROOT/'outputs/study_manifest.json',ROOT/'铜热膨胀研究.ipynb']
    # Preserve all existing figure artifacts, including prior versioned designs.
    protected_paths += [p for p in OUT.iterdir() if p.is_file() and not p.name.startswith(STEM)]
    protected={str(p.relative_to(ROOT)):digest(p) for p in protected_paths}
    manifest=json.loads((ROOT/'outputs/study_manifest.json').read_text(encoding='utf8'))
    by_id={job['id']:job for job in manifest}
    old=pd.read_csv(ROOT/'results/block_statistics.csv').set_index(['job','block_ps'])
    order=[f'prod_T{temp}_r{rep}' for temp in TEMPS for rep in [1,2,3]]
    order += [f'val_{kind}_T{temp}' for kind in ['size12','dt05','long'] for temp in [300,800]]
    assert len(order)==24 and set(order)==set(by_id)
    labels=[f'{temp} K · r{rep}' for temp in TEMPS for rep in [1,2,3]]
    labels += [f'{name} · {temp} K' for name in ['Size','Step','Long'] for temp in [300,800]]
    records,provenance,mapping=[],[],[]
    cols=['step','time_ps','temperature_K','pressure_bar','potential_eV','lx_A','ly_A','lz_A','volume_A3']
    max_relative_difference=0.
    for row_index,(job_id,label) in enumerate(zip(order,labels)):
        job=by_id[job_id]
        folder=ROOT/'outputs'/job['group']/job_id
        meta_path=folder/'metadata.json'
        raw_path=folder/'thermo.dat'
        meta=json.loads(meta_path.read_text(encoding='utf8'))
        assert meta['status']=='completed'
        for field in ['cells','dt','equil_ps','sample_ps','temperature']:
            assert meta[field]==job[field]
        assert meta['final_atoms']==4*job['cells']**3
        assert meta['sample_interval_ps']==.1
        is_long=job_id.startswith('val_long_')
        assert job['equil_ps']==(1000 if is_long else 200)
        assert job['sample_ps']==(1000 if is_long else 500)
        frame=pd.read_csv(raw_path,sep=r'\s+',comment='#',names=cols)
        sample=frame.loc[(frame.time_ps>job['equil_ps'])&
                         (frame.time_ps<=job['equil_ps']+job['sample_ps'])]
        assert len(sample)==round(job['sample_ps']/.1)
        assert np.isfinite(sample.to_numpy()).all()
        assert np.allclose(np.diff(sample.time_ps),.1,atol=1e-8,rtol=0)
        values=(sample.lx_A/job['cells']).to_numpy()
        mapping.append({'row':row_index,'label':label,'job':job_id,'group':job['group'],
                        'temperature_K':job['temperature'],'cells':job['cells'],
                        'equilibration_ps':job['equil_ps'],'production_ps':job['sample_ps']})
        for duration in BLOCKS:
            block_samples=round(duration/.1)
            count=len(values)//block_samples
            discarded=len(values)-count*block_samples
            block_means=values[:count*block_samples].reshape(count,block_samples).mean(axis=1)
            sem=float(block_means.std(ddof=1)/np.sqrt(count))
            recorded=old.loc[(job_id,duration)]
            difference=abs(sem-recorded.a_sem_A)/recorded.a_sem_A
            max_relative_difference=max(max_relative_difference,difference)
            assert count==recorded.blocks and np.isclose(sem,recorded.a_sem_A,rtol=1e-8,atol=1e-15), (job_id,duration,sem,recorded.a_sem_A)
            records.append({'job':job_id,'label':label,'group':job['group'],
                'temperature_K':job['temperature'],'block_ps':duration,'blocks':count,
                'samples_per_block':block_samples,'discarded_samples':discarded,'a_sem_A':sem})
        provenance.append({'job':job_id,'cells':job['cells'],'production_window_ps':
                           [job['equil_ps'],job['equil_ps']+job['sample_ps']],
                           'window_convention':'left boundary excluded; right boundary included',
                           'samples':len(values),'files_sha256':{
                               str(p.relative_to(ROOT)):digest(p) for p in [raw_path,meta_path]}})
    table=pd.DataFrame(records)
    absolute=table.pivot(index='job',columns='block_ps',values='a_sem_A').reindex(index=order,columns=BLOCKS)
    ratios=absolute.div(absolute[20],axis=0)
    assert np.isfinite(ratios.to_numpy()).all() and (ratios>0).all().all()
    assert np.array_equal(ratios[20].to_numpy(),np.ones(24))
    pair=pd.DataFrame(mapping).set_index('job').loc[order]
    pair['sem20_A']=absolute[20]
    pair['sem50_A']=absolute[50]
    pair['ratio_50_20']=ratios[50]
    median=float(ratios[50].median())
    observed_min=float(ratios[50].min())
    observed_max=float(ratios[50].max())
    expected=[1.018,.564,1.299]
    actual=[median,observed_min,observed_max]
    notes=[]
    if not np.allclose(actual,expected,rtol=0,atol=.0005):
        notes.append(f'Computed ratio summary differs from rounded checks: {actual}; actual values retained.')
    fonts={f.name for f in font_manager.fontManager.ttflist}
    font='Times New Roman' if 'Times New Roman' in fonts else 'STIXGeneral'
    if font!='Times New Roman':
        notes.append('Times New Roman unavailable; STIXGeneral used.')
    plt.rcParams.update({'font.family':font,'font.size':8.5,'axes.labelsize':9,
        'xtick.labelsize':8,'ytick.labelsize':8,'text.color':'#303840',
        'axes.labelcolor':'#303840','axes.edgecolor':'#303840',
        'xtick.color':'#303840','ytick.color':'#303840','axes.linewidth':.65,
        'figure.facecolor':'white','axes.facecolor':'white','mathtext.fontset':'stix',
        'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':True})
    fig=plt.figure(figsize=(185/25.4,160/25.4))
    axa=fig.add_axes([.170,.180,.350,.750])
    cax=fig.add_axes([.170,.080,.350,.013])
    axb=fig.add_axes([.685,.690,.285,.240])
    axc=fig.add_axes([.685,.190,.285,.32953125])
    log_ratios=np.log2(ratios.to_numpy())
    limit=float(np.abs(log_ratios).max())
    cmap=LinearSegmentedColormap.from_list('relative_sem',['#B3A3C8','#FAFBFC','#72B1A3'])
    image=axa.imshow(log_ratios,aspect='auto',cmap=cmap,vmin=-limit,vmax=limit,interpolation='nearest')
    axa.set_xticks(range(6),labels=BLOCKS)
    axa.set_yticks(range(24),labels=labels)
    axa.set_xlabel('Block duration / ps',labelpad=7)
    axa.tick_params(length=0,pad=5)
    for y in [2.5,5.5,8.5,11.5,14.5,19.5,21.5]:
        axa.axhline(y,color='white',lw=1.1)
    axa.axhline(17.5,color='#AEBAC3',lw=1.0)
    axa.add_patch(Rectangle((3.5,-.5),1,24,facecolor='none',edgecolor='#B77D55',lw=1.0,clip_on=False))
    for spine in axa.spines.values():
        spine.set_visible(False)
    candidates=[.25,.5,.75,1.,1.5,2.,4.]
    ratio_ticks=[value for value in candidates if abs(np.log2(value))<=limit]
    colorbar=fig.colorbar(image,cax=cax,orientation='horizontal',ticks=np.log2(ratio_ticks))
    colorbar.ax.set_xticklabels([f'{value:g}' for value in ratio_ticks])
    colorbar.set_label('SEM / SEM at 20 ps',fontsize=9,labelpad=4)
    colorbar.ax.tick_params(labelsize=8,length=2,width=.5,pad=2)
    colorbar.outline.set_visible(False)
    # Endpoint r1 curves use absolute SEM units and a logarithmic duration axis.
    for temp,color,marker,style in [(300,COLORS[0],'o','-'),(800,COLORS[-1],'^','--')]:
        values=absolute.loc[f'prod_T{temp}_r1']*1e6
        axb.plot(BLOCKS,values,marker=marker,ls=style,color=color,lw=1.2,
                 ms=3.8,mec='white',mew=.45,label=f'{temp} K')
        axb.scatter([20,50],values.loc[[20,50]],s=42,facecolors='none',edgecolors=color,lw=.7,zorder=5)
    axb.set_xscale('log')
    axb.xaxis.set_major_locator(FixedLocator(BLOCKS))
    axb.xaxis.set_major_formatter(FixedFormatter([str(v) for v in BLOCKS]))
    axb.xaxis.set_minor_locator(NullLocator())
    axb.set_xlim(.85,59)
    axb.set_xlabel('Block duration / ps',labelpad=6)
    axb.set_ylabel('Standard error / $10^{-6}$ Å',labelpad=6)
    axb.legend(frameon=False,fontsize=8,loc='center',bbox_to_anchor=(.76,.49),
               handlelength=1.9,handletextpad=.45)
    # Paired scatter: equal numerical limits and equal physical scale.
    lo=0
    hi=float(max(pair.sem20_A.max(),pair.sem50_A.max())*1e6*1.10)
    axc.plot([lo,hi],[lo,hi],color='#8B969E',ls=(0,(3,3)),lw=.8,zorder=1)
    for job_id,row in pair.iterrows():
        color=COLORS[TEMPS.index(int(row.temperature_K))]
        if row.group=='production':
            axc.scatter(row.sem20_A*1e6,row.sem50_A*1e6,s=21,marker='o',
                        facecolors=color,edgecolors='white',lw=.45,zorder=3)
        else:
            axc.scatter(row.sem20_A*1e6,row.sem50_A*1e6,s=27,marker='D',
                        facecolors='none',edgecolors=color,lw=.9,zorder=4)
    axc.set(xlim=(lo,hi),ylim=(lo,hi),xlabel='SEM at 20 ps / $10^{-6}$ Å',
            ylabel='SEM at 50 ps / $10^{-6}$ Å')
    axc.set_aspect('equal',adjustable='box')
    axc.text(.055,.95,f'Median ratio: {median:.3f}\nRange: {observed_min:.3f}–{observed_max:.3f}',
             transform=axc.transAxes,fontsize=8,va='top',linespacing=1.4)
    handles=[Line2D([],[],marker='o',ls='none',color='#63727B',ms=4,label='Production'),
             Line2D([],[],marker='D',ls='none',mfc='none',mec='#63727B',ms=4,label='Validation')]
    axc.legend(handles=handles,frameon=False,loc='lower right',fontsize=7.5,
                handletextpad=.25,borderaxespad=.4,labelspacing=.4)
    temperature_handles=[Line2D([],[],marker='o',ls='none',color=color,ms=3.7,label=f'{temp} K')
                         for temp,color in zip(TEMPS,COLORS)]
    fig.legend(handles=temperature_handles,loc='center',bbox_to_anchor=(.825,.601),
               frameon=False,ncol=3,fontsize=7.5,columnspacing=.7,handletextpad=.2,
               handlelength=.7,labelspacing=.6)
    for axis in [axb,axc]:
        axis.set_axisbelow(True)
        axis.grid(color='#E1E6EA',lw=.45,alpha=.65)
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
        axis.tick_params(direction='out',length=2.5,width=.6,pad=3)
    fig.text(.025,.955,'a',fontsize=12,fontweight='bold')
    fig.text(.600,.955,'b',fontsize=12,fontweight='bold')
    fig.text(.600,.555,'c',fontsize=12,fontweight='bold')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    artists=list(fig.texts)
    for axis in [axa,axb,axc,cax]:
        # Automatic locators may return additional, non-rendered ticks outside
        # the data limits. Check only tick labels actually inside each axis.
        xlo,xhi=sorted(axis.get_xlim())
        ylo,yhi=sorted(axis.get_ylim())
        xt=[label for value,label in zip(axis.get_xticks(),axis.get_xticklabels()) if xlo<=value<=xhi]
        yt=[label for value,label in zip(axis.get_yticks(),axis.get_yticklabels()) if ylo<=value<=yhi]
        artists += [axis.xaxis.label,axis.yaxis.label,*axis.texts,*xt,*yt]
    for text in artists:
        if not text.get_visible() or not text.get_text():
            continue
        box=text.get_window_extent(renderer)
        assert box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.width and box.y1<=fig.bbox.height,text.get_text()
    assert axc.get_xlim()==axc.get_ylim() and axc.get_aspect()==1.
    assert abs(axc.bbox.width/axc.bbox.height-1)<1e-8
    assert np.abs(log_ratios).max()<=limit
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f'{STEM}.pdf')
    fig.savefig(OUT/f'{STEM}.png',dpi=600)
    fig.savefig(OUT/f'{STEM}_preview.png',dpi=150)
    plt.close(fig)
    long_ratios=ratios.rename_axis(index='job',columns='block_ps').stack().rename('sem_ratio_to_20ps').reset_index()
    table.merge(long_ratios,on=['job','block_ps'],validate='one_to_one').to_csv(OUT/f'{STEM}_blocks.csv',index=False)
    pair.reset_index().to_csv(OUT/f'{STEM}_pairs.csv',index=False)
    pd.DataFrame(mapping).to_csv(OUT/f'{STEM}_labels.csv',index=False)
    caption='''**Figure 6. Sensitivity of lattice-parameter standard-error estimates to statistical block duration.** All 24 completed trajectories are analyzed separately: 18 production trajectories followed by six endpoint validation trajectories, with groups separated in (a). Each standard error is the sample standard deviation of nonoverlapping block means divided by the square root of the complete-block count. (a) Ratios to each trajectory's own 20 ps standard error, encoded by a shared symmetric log₂ ratio color scale; colorbar labels show the corresponding ratios. The copper outline marks the 20 ps column, whose ratio of one follows from normalization. Lavender and teal indicate ratios below and above one, not pass/fail classifications. (b) Absolute standard errors for production repeat r1 at 300 and 800 K; circles and triangles identify the two temperatures, open rings highlight 20 and 50 ps, and connecting segments are guides to the eye on the logarithmic duration axis. (c) Paired 20 and 50 ps estimates for all trajectories, with equal axes and a 1:1 reference line. Colors identify temperature; filled circles and open diamonds identify production and validation trajectories. The annotated range is the observed range of ratios, not a confidence interval. Standard trajectories use 500 ps after 200 ps equilibration; long-duration checks use 1,000 ps after 1,000 ps equilibration. Longer blocks reduce the number of available block means and can increase variability in the estimated error. Limited sensitivity within these tested block durations does not establish complete independence, universal convergence, or the absence of slower correlations.\n'''
    (ROOT/'docs/figure6_caption_v1.md').write_text(caption,encoding='utf8')
    for path,sha in protected.items():
        assert digest(ROOT/path)==sha
    for record in provenance:
        for path,sha in record['files_sha256'].items():
            assert digest(ROOT/path)==sha
    audit={'protected_files_sha256':protected,'sources':provenance,'block_durations_ps':BLOCKS,
        'definition':'sample SD (ddof=1) of complete nonoverlapping block means / sqrt(number of blocks)',
        'normalization':'within trajectory SEM(duration)/SEM(20ps)',
        'color_transform':'log2(SEM ratio); one symmetric scale centered on zero',
        'color_limits_log2':[-limit,limit],'colorbar_ratio_ticks':ratio_ticks,
        'discarded_samples_total_across_all_block_calculations':int(table.discarded_samples.sum()),
        'max_relative_difference_from_existing_block_file':max_relative_difference,
        'ratio_summary':{'median':median,'observed_min':observed_min,'observed_max':observed_max},
        'font':font,'font_file':font_manager.findfont(font),'size_mm':[185,160],
        'png_dpi':600,'preview_dpi':150,'notes':notes,
        'verification':['All 144 SEM values and block counts agree with previous records',
                        'All 24 raw production windows verified','All data, prior figures and Notebook unchanged',
                        'No heatmap values clipped','20ps ratios exactly one',
                        'Paired scatter equal limits and equal physical aspect','Text extents inside canvas']}
    (OUT/f'{STEM}_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
    print(json.dumps({'font':font,'ratios':audit['ratio_summary'],'color_limits_log2':[-limit,limit],
                      'discarded_samples':audit['discarded_samples_total_across_all_block_calculations'],
                      'max_relative_difference':max_relative_difference,'notes':notes,
                      'figure':str(OUT/f'{STEM}.png')},indent=2))


if __name__=='__main__':
    main()
