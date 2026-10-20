"""Sensitivity figure: verified settings matrix and aligned coefficient changes."""
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
from matplotlib.colors import to_rgb
from matplotlib.patches import Rectangle

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/paper_figures'
STEM='figure4_settings_sensitivity_v1'
CONDITIONS=[('baseline','Baseline','#68717A'),('size12','Larger cell','#9A8AB4'),
            ('dt05','Smaller timestep','#D1A16D'),('long','Longer duration','#839BA8')]
EXPECTED={'baseline':(4000,1.,200,500),'size12':(6912,1.,200,500),
          'dt05':(4000,.5,200,500),'long':(4000,1.,1000,1000)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tint(color,amount=.26):
    return tuple(1-amount*(1-c) for c in to_rgb(color))


def main():
    protected_files=['results/run_summary.csv','results/validation_comparison.csv',
                     'outputs/study_manifest.json','铜热膨胀研究.ipynb']
    protected={name:digest(ROOT/name) for name in protected_files}
    summary=pd.read_csv(ROOT/'results/run_summary.csv').set_index('id')
    validation=pd.read_csv(ROOT/'results/validation_comparison.csv').set_index('check')
    manifest=json.loads((ROOT/'outputs/study_manifest.json').read_text(encoding='utf8'))
    by_id={j['id']:j for j in manifest}
    raw_cols=['step','time_ps','temperature_K','pressure_bar','potential_eV',
              'lx_A','ly_A','lz_A','volume_A3']
    rows,provenance=[],[]
    for key,label,color in CONDITIONS:
        means=[]
        parameters=[]
        for temp in [300,800]:
            job_id=f'prod_T{temp}_r1' if key=='baseline' else f'val_{key}_T{temp}'
            job=by_id[job_id]
            group='production' if key=='baseline' else 'validation'
            folder=ROOT/'outputs'/group/job_id
            meta_path=folder/'metadata.json'
            input_path=folder/'input.in'
            raw_path=folder/'thermo.dat'
            meta=json.loads(meta_path.read_text(encoding='utf8'))
            assert meta['status']=='completed'
            for field in ['cells','dt','equil_ps','sample_ps','temperature']:
                assert meta[field]==job[field]
            atoms=4*job['cells']**3
            assert meta['final_atoms']==atoms
            settings=(atoms,job['dt']*1000,job['equil_ps'],job['sample_ps'])
            assert settings==EXPECTED[key],(key,settings)
            parameters.append(settings)
            source=input_path.read_text(encoding='ascii')
            assert re.search(r'^units\s+metal\s*$',source,re.M)
            timestep=float(re.search(r'^timestep\s+(\S+)',source,re.M).group(1))
            assert timestep==job['dt']
            frame=pd.read_csv(raw_path,sep=r'\s+',comment='#',names=raw_cols)
            sample=frame.loc[(frame.time_ps>job['equil_ps']) &
                             (frame.time_ps<=job['equil_ps']+job['sample_ps'])]
            expected_count=round(job['sample_ps']/.1)
            assert len(sample)==expected_count and np.isfinite(sample.to_numpy()).all()
            mean_a=float((sample.lx_A/job['cells']).mean())
            assert np.isclose(mean_a,summary.loc[job_id,'a_A'],rtol=1e-12)
            means.append(mean_a)
            provenance.append({'id':job_id,'temperature_K':temp,
                'production_window_ps':[job['equil_ps'],job['equil_ps']+job['sample_ps']],
                'window_convention':'left boundary excluded; right boundary included',
                'samples':expected_count,'lattice_mean_A':mean_a,
                'files_sha256':{str(p.relative_to(ROOT)):digest(p) for p in [raw_path,meta_path,input_path]}})
        assert parameters[0]==parameters[1]
        atoms,dt,eq,production=parameters[0]
        alpha=(means[1]/means[0]-1)/500
        rows.append({'check':key,'condition':label,'atoms':atoms,'timestep_fs':dt,
                     'equilibration_ps':eq,'production_ps':production,
                     'a300_A':means[0],'a800_A':means[1],'alpha_Kinv':alpha,'color':color})
    baseline=rows[0]['alpha_Kinv']
    for row in rows:
        row['baseline_r1_alpha_Kinv']=baseline
        row['change_percent']=100*(row['alpha_Kinv']/baseline-1)
        row['practical_limit_percent']=2.
        if row['check']!='baseline':
            check=validation.loc[row['check']]
            assert np.isclose(row['alpha_Kinv'],check.alpha_Kinv,rtol=1e-11)
            assert np.isclose(baseline,check.baseline_r1_alpha_Kinv,rtol=1e-11)
            assert np.isclose(row['change_percent'],check.change_percent,atol=1e-9,rtol=0)
    expected_changes={'size12':.0680,'dt05':-.0087,'long':.0280}
    for row in rows[1:]:
        assert np.isclose(row['change_percent'],expected_changes[row['check']],atol=.00005,rtol=0)
    fonts={f.name for f in font_manager.fontManager.ttflist}
    font='Times New Roman' if 'Times New Roman' in fonts else 'STIXGeneral'
    plt.rcParams.update({'font.family':font,'font.size':8.5,'axes.labelsize':9,
        'text.color':'#303840','axes.labelcolor':'#303840','axes.edgecolor':'#303840',
        'xtick.color':'#303840','ytick.color':'#303840','xtick.labelsize':8,
        'figure.facecolor':'white','axes.facecolor':'white','axes.linewidth':.65,
        'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':True})
    fig=plt.figure(figsize=(185/25.4,100/25.4))
    axa=fig.add_axes([.020,.23,.565,.65])
    axb=fig.add_axes([.650,.23,.330,.65])
    axa.set(xlim=(0,1),ylim=(-.6,4.1))
    axb.set(ylim=(-.6,4.1),xlim=(-.045,.110))
    axa.axis('off')
    # Panel headings and table headers are separate to retain a clean hierarchy.
    fig.text(.020,.94,'a',fontsize=12,fontweight='bold',va='baseline')
    fig.text(.054,.94,'Simulation settings',fontsize=10.5,fontweight='bold',va='baseline')
    fig.text(.650,.94,'b',fontsize=12,fontweight='bold',va='baseline')
    fig.text(.684,.94,'Coefficient change',fontsize=10.5,fontweight='bold',va='baseline')
    bounds=[0,.29,.43,.61,.82,1.]
    centers=[(bounds[i]+bounds[i+1])/2 for i in range(5)]
    headers=['Condition','Atoms','Timestep\n/ fs','Equilibration\n/ ps','Production\n/ ps']
    for i,header in enumerate(headers):
        axa.text(.008 if i==0 else centers[i],3.82,header,ha='left' if i==0 else 'center',
                 va='center',fontsize=8.5,fontweight='bold',linespacing=1.25)
    axa.hlines(3.42,0,1,color='#D5DDE2',lw=.75)
    changed={'baseline':[],'size12':[1],'dt05':[2],'long':[3,4]}
    annotations=[]
    for idx,row in enumerate(rows):
        y=3-idx
        axa.text(.008,y,row['condition'],va='center',ha='left',fontsize=8.5)
        values=[f"{row['atoms']:,}",f"{row['timestep_fs']:.1f}",
                f"{row['equilibration_ps']:,}",f"{row['production_ps']:,}"]
        for col,value in enumerate(values,1):
            face=tint(row['color']) if col in changed[row['check']] else '#F3F5F6'
            axa.add_patch(Rectangle((bounds[col]+.008,y-.31),
                bounds[col+1]-bounds[col]-.016,.62,facecolor=face,edgecolor='none',zorder=1))
            axa.text(centers[col],y,value,ha='center',va='center',fontsize=9,
                     fontweight='bold' if col in changed[row['check']] else 'normal',zorder=2)
        if idx<3:
            axa.hlines(y-.5,0,1,color='#E1E6EA',lw=.45)
        if idx==0:
            axb.scatter(0,y,s=20,color='#68717A',zorder=4)
            annotations.append(axb.text(.005,y,'Reference',ha='left',va='center',fontsize=8.5,color='#68717A'))
        else:
            value=row['change_percent']
            axb.barh(y,value,height=.36,color=row['color'],edgecolor='none',zorder=3)
            label=f'{value:+.4f}%'.replace('-','−')
            annotations.append(axb.text(value+(.003 if value>0 else -.003),y,label,
                 ha='left' if value>0 else 'right',va='center',fontsize=8.5))
    axb.axvline(0,color='#87939C',lw=.8,zorder=2)
    axb.set_axisbelow(True)
    axb.set_xticks([-.04,0,.04,.08],labels=['−0.04','0','0.04','0.08'])
    axb.grid(axis='x',color='#E8EDF0',lw=.55)
    axb.set_yticks([])
    axb.set_xlabel('Change relative to baseline r1 / %',labelpad=8,fontsize=9)
    axb.tick_params(axis='x',length=3,width=.6,pad=5)
    for name in ['top','right','left']:
        axb.spines[name].set_visible(False)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    all_text=list(fig.texts)+list(axa.texts)+list(axb.texts)+[axb.xaxis.label,*axb.get_xticklabels()]
    for artist in all_text:
        box=artist.get_window_extent(renderer)
        assert box.x0>=0 and box.y0>=0 and box.x1<=fig.bbox.width and box.y1<=fig.bbox.height,artist.get_text()
    # Bar-value labels must fit inside the result panel, not spill into the matrix.
    for artist in annotations:
        box=artist.get_window_extent(renderer)
        assert box.x0>=axb.bbox.x0 and box.x1<=axb.bbox.x1,artist.get_text()
    assert np.allclose([axa.transData.transform((0,y))[1] for y in range(4)],
                       [axb.transData.transform((0,y))[1] for y in range(4)])
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f'{STEM}.pdf')
    fig.savefig(OUT/f'{STEM}.png',dpi=600)
    fig.savefig(OUT/f'{STEM}_preview.png',dpi=150)
    plt.close(fig)
    pd.DataFrame(rows).drop(columns=['color']).to_csv(OUT/f'{STEM}_data.csv',index=False)
    caption='''**Figure 4. Sensitivity to simulation settings.** (a) Baseline and modified settings, with colored cells identifying the parameters changed from the baseline. Unchanged values use neutral backgrounds; colors identify conditions and are not a numerical scale across unlike units. (b) Relative changes in the engineering 300–800 K average linear expansion coefficient, calculated as 100 × (ᾱ_condition/ᾱ_baseline,r1 − 1). Every comparison uses production repeat r1 at 300 and 800 K as the baseline, rather than the three-repeat mean. Each modified condition contains one trajectory at each endpoint, and no independent-repeat uncertainty intervals are available for these checks. The horizontal axis is expanded to resolve small changes; the prespecified ±2% practical sensitivity target lies outside the displayed range. Passing these limited checks supports robustness under the tested changes, not universal convergence or an exact bound on numerical error.\n'''
    (ROOT/'docs/figure4_caption_v1.md').write_text(caption,encoding='utf8')
    for path,sha in protected.items():
        assert digest(ROOT/path)==sha
    for record in provenance:
        for path,sha in record['files_sha256'].items():
            assert digest(ROOT/path)==sha
    report={'protected_files_sha256':protected,'trajectories':provenance,
        'baseline':'Production repeat r1; independently averaged endpoint trajectories at 300 and 800 K',
        'coefficient_definition':'(a800/a300-1)/500 K',
        'relative_change_definition':'100*(alpha_condition/alpha_baseline_r1-1)',
        'timestep_conversion':'LAMMPS metal units: ps to fs by multiplication by 1000',
        'parameter_verification':'Manifest and metadata agree; atom counts match final_atoms; input units and timestep checked',
        'endpoint_checks':'Computed means and coefficients agree with run_summary and validation_comparison; settings and rounded changes agree with expected checks',
        'independent_validation_repeats_per_endpoint':1,'uncertainty_intervals':'Not estimated or plotted',
        'practical_target_percent':2,'displayed_x_limits_percent':[-.045,.110],
        'font':font,'font_file':font_manager.findfont(font),'size_mm':[185,100],
        'png_dpi':600,'preview_dpi':150,'rows':rows,
        'notes':[] if font=='Times New Roman' else ['Times New Roman unavailable; STIXGeneral used.'],
        'checks':['Eight raw trajectory averages verified','Raw files, inputs, metadata, summaries and Notebook unchanged',
                  'Matrix and chart row centers aligned exactly','Text extents inside canvas','Result labels inside chart boundaries']}
    (OUT/f'{STEM}_audit.json').write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps({'baseline_alpha_Kinv':baseline,'changes_percent':{row['check']:row['change_percent'] for row in rows},
                      'font':font,'notes':report['notes'],'figure':str(OUT/f'{STEM}.png')},indent=2))


if __name__=='__main__':
    main()
