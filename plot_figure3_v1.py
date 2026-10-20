"""Compare copper MD expansion with the revised NIST SRM 736 certificate."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/paper_figures'
STEM='figure3_reference_composite_v1'
TEMPS=[300,400,500,600,700,800]
COLORS=['#477BA8','#4C91B0','#55A5B3','#62B5AC','#79BEA5','#95C7A4']
TCRIT=4.3026527299
CERTIFICATE_URL='https://tsapps.nist.gov/srmext/certificates/archives/736.pdf'
CERTIFICATE_SHA256='4a76360a7bde1f58a53d38c6e2af7bfe8bda80a53c624e4913020666d7961bd3'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary=ROOT/'results/run_summary.csv'
    reference=ROOT/'data/reference/nist_srm736.csv'
    notebook=ROOT/'铜热膨胀研究.ipynb'
    protected={str(p.relative_to(ROOT)):digest(p) for p in [summary,reference,notebook]}
    prod=pd.read_csv(summary).query("group == 'production'")
    assert len(prod)==18 and not prod.duplicated(['temperature','replicate']).any()
    assert prod.groupby('temperature').size().eq(3).all()
    assert sorted(prod.temperature.unique())==TEMPS
    provenance=[]
    cols=['step','time_ps','temperature_K','pressure_bar','potential_eV','lx_A','ly_A','lz_A','volume_A3']
    for row in prod.itertuples(index=False):
        folder=ROOT/'outputs/production'/row.id
        assert json.loads((folder/'metadata.json').read_text())['status']=='completed'
        raw=folder/'thermo.dat'
        frame=pd.read_csv(raw,sep=r'\s+',comment='#',names=cols)
        sample=frame.loc[(frame.time_ps>200)&(frame.time_ps<=700)]
        assert row.equil_ps==200 and row.sample_ps==500 and row.cells==10
        assert len(sample)==5000 and np.isfinite(sample.to_numpy()).all()
        assert np.isclose((sample.lx_A/10).mean(),row.a_A,rtol=1e-12)
        provenance.append({'path':str(raw.relative_to(ROOT)),'sha256':digest(raw)})
    certificate=ROOT/'.tmp/nist_srm736.pdf'
    assert digest(certificate)==CERTIFICATE_SHA256
    ref=pd.read_csv(reference).set_index('temperature_K').loc[TEMPS].copy()
    # Values independently read from p.1 of the archived 7 October 1990 certificate.
    assert ref.expansion_from_293_ppm.tolist()==[117,1831,3620,5480,7411,9420]
    e293=ref.expansion_from_293_ppm*1e-6
    ref_eps=(1+e293)/(1+e293.loc[300])-1
    a=prod.pivot(index='temperature',columns='replicate',values='a_A').sort_index()
    eps=a.div(a.loc[300],axis=1)-1
    eps_mean=eps.mean(axis=1)
    eps_ci=TCRIT*eps.std(axis=1,ddof=1)/np.sqrt(3)
    residual=eps.sub(ref_eps,axis=0)*1e6
    res_mean=residual.mean(axis=1)
    res_ci=TCRIT*residual.std(axis=1,ddof=1)/np.sqrt(3)
    coefficients=eps.loc[800]/500*1e6
    estimate=float(coefficients.mean())
    halfwidth=float(TCRIT*coefficients.std(ddof=1)/np.sqrt(3))
    refalpha=float(ref_eps.loc[800]/500*1e6)
    bias=100*(estimate/refalpha-1)
    assert eps.loc[300].eq(0).all() and residual.loc[300].eq(0).all()
    assert np.allclose(res_ci,eps_ci*1e6)
    assert np.isclose(estimate,17.888,atol=.0005,rtol=0)
    assert np.isclose(halfwidth,.015,atol=.0005,rtol=0)
    assert np.isclose(refalpha,18.6038,atol=.00005,rtol=0)
    assert np.isclose(bias,-3.85,atol=.005,rtol=0)

    available={f.name for f in font_manager.fontManager.ttflist}
    font='Times New Roman' if 'Times New Roman' in available else 'STIXGeneral'
    plt.rcParams.update({
        'font.family':font,'font.size':8.5,'axes.labelsize':9.5,
        'xtick.labelsize':8.5,'ytick.labelsize':8.5,'text.color':'#303840',
        'axes.labelcolor':'#303840','axes.edgecolor':'#303840',
        'xtick.color':'#303840','ytick.color':'#303840','axes.linewidth':.65,
        'axes.spines.top':False,'axes.spines.right':False,
        'figure.facecolor':'white','axes.facecolor':'white','mathtext.fontset':'stix',
        'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':True,
    })
    fig=plt.figure(figsize=(185/25.4,132/25.4))
    axa=fig.add_axes([.103,.14,.390,.770])
    axb=fig.add_axes([.650,.620,.330,.290])
    axc=fig.add_axes([.650,.205,.330,.245])
    for axis in [axa,axb,axc]:
        axis.set_axisbelow(True)
        axis.tick_params(direction='out',length=3,width=.6,pad=4)
    axa.grid(axis='y',color='#E8EDF0',lw=.55)
    axa.plot(TEMPS,100*ref_eps,'D--',color='#59636D',lw=1.35,ms=3.5,
             mfc='white',mec='#59636D',mew=.7,zorder=2)
    axa.errorbar(TEMPS,100*eps_mean,yerr=100*eps_ci,fmt='-',color='#365F78',
                 lw=1.5,elinewidth=.8,capsize=2.2,capthick=.8,zorder=3)
    axa.scatter(TEMPS,100*eps_mean,c=COLORS,s=24,edgecolors='white',linewidths=.55,zorder=4)
    axa.set(xlim=(278,822),ylim=(-.035,1.015),xlabel='Temperature / K',
            ylabel='Linear expansion relative to 300 K / %')
    axa.set_xticks(TEMPS)
    axa.set_yticks([0,.2,.4,.6,.8,1.0])
    axa.xaxis.labelpad=8
    axa.yaxis.labelpad=8
    handles=[Line2D([],[],color='#365F78',marker='o',ms=4,lw=1.5,label='MD'),
             Line2D([],[],color='#59636D',marker='D',mfc='white',ms=3.5,
                    lw=1.35,ls='--',label='NIST SRM 736')]
    axa.legend(handles=handles,frameon=False,loc='lower right',fontsize=8.5,
               handlelength=2.2,borderaxespad=.65,labelspacing=.6)
    axa.text(.065,.89,'300–800 K relative deviation',transform=axa.transAxes,
             fontsize=8.5,color='#B77D55',va='top')
    axa.text(.065,.835,f'{bias:.2f}%'.replace('-','−'),transform=axa.transAxes,
             fontsize=13,color='#B77D55',va='top')

    # Treat the reference values as fixed here; do not invent reference uncertainty.
    axb.axvline(0,color='#89949B',lw=.75,zorder=1)
    axb.grid(axis='y',color='#E8EDF0',lw=.5)
    for t,color in zip(TEMPS,COLORS):
        axb.errorbar(res_mean.loc[t],t,xerr=res_ci.loc[t],fmt='o',color=color,
                     ms=4.3,mec='white',mew=.5,elinewidth=.85,capsize=2,capthick=.8,zorder=3)
    axb.set(xlim=(-400,30),ylim=(273,827),xlabel='MD − NIST / ppm',ylabel='Temperature / K')
    axb.set_xticks([-400,-200,0])
    axb.set_yticks(TEMPS)
    axb.xaxis.labelpad=6
    axb.yaxis.labelpad=6

    # Individual coefficients keep their actual x values; only category y is offset.
    repeat_offsets=[-.14,0,.14]
    axc.grid(axis='x',color='#E8EDF0',lw=.55)
    axc.errorbar(estimate,1,xerr=halfwidth,fmt='o',color='#62A99A',ms=5,
                 elinewidth=1,capsize=3,capthick=.8,mec='white',mew=.5,zorder=4)
    for rep,offset,marker in zip([1,2,3],repeat_offsets,['s','D','^']):
        axc.scatter(coefficients.loc[rep],1+offset,s=18,marker=marker,
                    facecolors='none',edgecolors='#536C68',linewidths=.65,zorder=5)
    axc.scatter(refalpha,0,marker='D',s=25,c='#59636D',edgecolors='white',linewidths=.5,zorder=4)
    axc.set(xlim=(17.76,18.72),ylim=(-.65,1.65))
    axc.set_yticks([1,0],labels=['MD','NIST'])
    axc.set_xticks([17.8,18.1,18.4,18.7])
    axc.set_xlabel('Average linear expansion coefficient\n/ $10^{-6}$ K$^{-1}$',labelpad=7,fontsize=9)
    axc.text(estimate+.07,1.27,f'{estimate:.3f} ± {halfwidth:.3f}',fontsize=8,color='#438C7E',ha='left')
    axc.text(refalpha,-.40,f'{refalpha:.4f}',fontsize=8,color='#59636D',ha='center')
    fig.text(.025,.940,'a',fontsize=12,fontweight='bold')
    fig.text(.565,.940,'b',fontsize=12,fontweight='bold')
    fig.text(.565,.485,'c',fontsize=12,fontweight='bold')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    essential=list(fig.texts)
    for axis in [axa,axb,axc]:
        essential += [axis.xaxis.label,axis.yaxis.label,*axis.texts,
                      *axis.get_xticklabels(),*axis.get_yticklabels()]
    for text in essential:
        if not text.get_visible() or not text.get_text():
            continue
        box=text.get_window_extent(renderer)
        assert box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.width and box.y1<=fig.bbox.height,text.get_text()
    assert ((eps_mean-eps_ci)*100).min()>axa.get_ylim()[0]
    assert ((eps_mean+eps_ci)*100).max()<axa.get_ylim()[1]
    assert (res_mean-res_ci).min()>axb.get_xlim()[0]
    assert (res_mean+res_ci).max()<axb.get_xlim()[1]
    assert estimate-halfwidth>axc.get_xlim()[0] and refalpha<axc.get_xlim()[1]
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f'{STEM}.pdf')
    fig.savefig(OUT/f'{STEM}.png',dpi=600)
    fig.savefig(OUT/f'{STEM}_preview.png',dpi=150)
    plt.close(fig)

    table=pd.DataFrame({'temperature_K':TEMPS,'md_expansion_mean_percent':100*eps_mean,
        'md_expansion_ci95_halfwidth_percent':100*eps_ci,'nist_expansion_percent':100*ref_eps,
        'residual_mean_ppm':res_mean,'residual_ci95_halfwidth_ppm':res_ci})
    for rep in [1,2,3]:
        table[f'md_expansion_r{rep}_percent']=100*eps[rep]
        table[f'residual_r{rep}_ppm']=residual[rep]
    table.to_csv(OUT/f'{STEM}_data.csv',index=False)
    coefficient_table=pd.DataFrame({'replicate':[1,2,3],'alpha_micro_Kinv':coefficients.to_numpy()})
    coefficient_table.to_csv(OUT/f'{STEM}_replicates.csv',index=False)
    pd.DataFrame({'source':['MD','NIST SRM 736'],
        'alpha_micro_Kinv':[estimate,refalpha],
        'ci95_halfwidth_micro_Kinv':[halfwidth,np.nan],
        'reference_uncertainty_propagated':[False,False]}).to_csv(OUT/f'{STEM}_coefficient_summary.csv',index=False)
    caption='''**Figure 3. Comparison of copper thermal expansion predicted by the Mishin EAM1 potential with the revised NIST SRM 736 reference.** (a) Linear expansion relative to 300 K, with MD means obtained from three independently initialized production trajectories per temperature. (b) Mean MD-minus-NIST residuals in ppm, with temperature increasing upward. (c) Engineering interval-average coefficients, calculated for each MD repeat as [ā(800 K)/ā(300 K) − 1]/500 K and then averaged. Error bars in all panels are two-sided 95% Student-t confidence intervals with two degrees of freedom and describe MD sampling variability only; reference uncertainty has not been propagated, and no reference interval is implied to be zero. MD strains are normalized separately to each repeat's own 300 K mean. NIST fractional length changes relative to 293 K are converted to the 300 K reference using the revised certificate dated 7 October 1990, p. 1. The coincident curves and zero residual with zero interval at 300 K follow from normalization, not independent validation. In (c), open squares, diamonds, and triangles represent r1, r2, and r3, with vertical category offsets of −0.14, 0, and +0.14 for visibility only; coefficient values are not shifted. Temperature colors are consistent between (a) and (b). Connecting lines are guides to the eye, and some intervals may be smaller than the markers. All MD averages exclude 200 ps of equilibration and use the subsequent 500 ps production segment.\n'''
    (ROOT/'docs/figure3_caption_v1.md').write_text(caption,encoding='utf8')
    for path,sha in protected.items():
        assert digest(ROOT/path)==sha
    for record in provenance:
        assert digest(ROOT/record['path'])==record['sha256']
    report={'protected_inputs_sha256':protected,'raw_production_files':provenance,
        'reference':{'url':CERTIFICATE_URL,'revision':'1990-10-07','page':1,
                     'archived_pdf_sha256':CERTIFICATE_SHA256,
                     'table_ppm':[117,1831,3620,5480,7411,9420],
                     'verification':'Certificate p.1 visually checked against CSV; official URL accessed during generation.'},
        'normalization_md':'a_mean_r(T)/a_mean_r(300)-1',
        'normalization_reference':'(1+ppm(T)*1e-6)/(1+ppm(300)*1e-6)-1',
        'residual':'1e6*(strain_r(T)-reference_strain(T))',
        'coefficient':'strain_r(800)/500; arithmetic mean of three independent estimates',
        'uncertainty':'Two-sided 95% Student-t, df=2; t=4.3026527299; MD sampling only',
        'reference_uncertainty':'Not propagated; omitted intervals do not imply zero reference uncertainty',
        'font':font,'font_file':font_manager.findfont(font),'size_mm':[185,132],
        'png_dpi':600,'preview_dpi':150,'category_offsets':repeat_offsets,
        'results':{'md_coefficient_micro_Kinv':estimate,'ci95_halfwidth_micro_Kinv':halfwidth,
                   'nist_coefficient_micro_Kinv':refalpha,'relative_deviation_percent':bias,
                   'residual_800_ppm':float(res_mean.loc[800])},
        'endpoint_checks':'All calculated quantities agree with expected rounded checks',
        'notes':[] if font=='Times New Roman' else ['Times New Roman unavailable; STIXGeneral substituted.'],
        'verification':['18 production means checked against raw samples','18 raw files unchanged',
                        'Summary, reference CSV and Notebook unchanged','Text extents inside canvas','Intervals inside axes limits']}
    (OUT/f'{STEM}_audit.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps({'font':font,'results':report['results'],'notes':report['notes'],
                      'figure':str(OUT/f'{STEM}.png')},indent=2))


if __name__=='__main__':
    main()
