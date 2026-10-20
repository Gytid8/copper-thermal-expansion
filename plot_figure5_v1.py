"""Full recorded r1 trajectories: 3D volume evolution and endpoint diagnostics."""
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
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/paper_figures'
STEM='figure5_trajectory_diagnostics_v1'
TEMPS=[300,400,500,600,700,800]
COLORS=['#477BA8','#4C91B0','#55A5B3','#62B5AC','#79BEA5','#95C7A4']
COLS=['step','time_ps','temperature_K','pressure_bar','potential_eV','lx_A','ly_A','lz_A','volume_A3']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rolling_separately(time,values):
    """Trailing 50-record display mean; reset at the production boundary."""
    output=np.full(len(time),np.nan)
    for mask in [time<=200,time>200]:
        output[mask]=pd.Series(values[mask]).rolling(50,min_periods=1).mean().to_numpy()
    return output


def main():
    protected={name:digest(ROOT/name) for name in ['results/run_summary.csv','铜热膨胀研究.ipynb']}
    summary=pd.read_csv(ROOT/'results/run_summary.csv').set_index('id')
    frames,means,provenance={}, {}, []
    for temp in TEMPS:
        job_id=f'prod_T{temp}_r1'
        folder=ROOT/'outputs/production'/job_id
        meta_path=folder/'metadata.json'
        input_path=folder/'input.in'
        raw_path=folder/'thermo.dat'
        meta=json.loads(meta_path.read_text(encoding='utf8'))
        assert meta['status']=='completed' and meta['replicate']==1
        assert meta['temperature']==temp and meta['final_atoms']==4000
        assert meta['equil_ps']==200 and meta['sample_ps']==500
        assert meta['sample_interval_ps']==.1 and meta['dt']==.001
        source=input_path.read_text(encoding='ascii')
        assert re.search(r'^units\s+metal\s*$',source,re.M)
        frame=pd.read_csv(raw_path,sep=r'\s+',comment='#',names=COLS)
        assert len(frame)==7001 and np.isfinite(frame.to_numpy()).all()
        assert np.allclose(frame.time_ps,np.arange(7001)*.1,rtol=0,atol=1e-8)
        assert np.array_equal(frame.step,np.arange(7001)*100)
        production=frame.loc[(frame.time_ps>200)&(frame.time_ps<=700)]
        assert len(production)==5000
        assert np.isclose(production.volume_A3.mean(),summary.loc[job_id,'volume_A3'],rtol=1e-12)
        assert np.isclose((production.lx_A/10).mean(),summary.loc[job_id,'a_A'],rtol=1e-12)
        frames[temp]=frame
        means[temp]={'temperature_mean_K':float(production.temperature_K.mean()),
                     'pressure_mean_bar':float(production.pressure_bar.mean()),
                     'potential_mean_eV_atom':float(production.potential_eV.mean()/4000),
                     'volume_mean_A3':float(production.volume_A3.mean())}
        provenance.append({'id':job_id,'temperature_K':temp,'atoms':4000,
            'recorded_time_range_ps':[0,700],'recorded_points':7001,'production_points':5000,
            'initial_record':'Recorded step 0 retained; no samples were synthesized',
            'sha256':{str(p.relative_to(ROOT)):digest(p) for p in [raw_path,meta_path,input_path]}})
    denominator=means[300]['volume_mean_A3']
    for temp,frame in frames.items():
        frame['volume_change_percent']=100*(frame.volume_A3/denominator-1)
        frame['temperature_deviation_K']=frame.temperature_K-temp
        frame['pressure_GPa']=frame.pressure_bar*1e-4
        frame['potential_deviation_eV_atom']=frame.potential_eV/4000-means[temp]['potential_mean_eV_atom']
        assert abs(frame.loc[frame.time_ps>200,'potential_deviation_eV_atom'].mean())<1e-12
    fonts={f.name for f in font_manager.fontManager.ttflist}
    font='Times New Roman' if 'Times New Roman' in fonts else 'STIXGeneral'
    plt.rcParams.update({'font.family':font,'font.size':8.5,'axes.labelsize':9,
        'xtick.labelsize':8,'ytick.labelsize':8,'text.color':'#303840',
        'axes.labelcolor':'#303840','axes.edgecolor':'#303840',
        'xtick.color':'#303840','ytick.color':'#303840','axes.linewidth':.65,
        'figure.facecolor':'white','axes.facecolor':'white','mathtext.fontset':'stix',
        'pdf.fonttype':42,'ps.fonttype':42,'axes.unicode_minus':True,
        'path.simplify':False})
    fig=plt.figure(figsize=(185/25.4,155/25.4))
    ax=fig.add_axes([.05,.435,.90,.535],projection='3d')
    ax.set_proj_type('ortho')
    ax.set_box_aspect((3.8,2.0,1.55),zoom=1.08)
    elevation,azimuth=22,-65
    ax.view_init(elev=elevation,azim=azimuth)
    global_min=min(frame.volume_change_percent.min() for frame in frames.values())
    global_max=max(frame.volume_change_percent.max() for frame in frames.values())
    zpad=.08*(global_max-global_min)
    zlo,zhi=global_min-zpad,global_max+zpad
    for temp,color in zip(TEMPS,COLORS):
        frame=frames[temp]
        eq=frame.time_ps<=200
        post=frame.time_ps>=200  # Retain the recorded boundary point for line continuity only.
        ax.plot(frame.time_ps[eq],np.full(eq.sum(),temp),frame.volume_change_percent[eq],
                color='#C5CDD2',lw=.48,alpha=.85)
        ax.plot(frame.time_ps[post],np.full(post.sum(),temp),frame.volume_change_percent[post],
                color=color,lw=.62,alpha=.95)
    # A labeled, neutral boundary plane is a reference, not a fitted data surface.
    plane=[[(200,300,zlo),(200,800,zlo),(200,800,zhi),(200,300,zhi)]]
    ax.add_collection3d(Poly3DCollection(plane,facecolor='#89949B',alpha=.035,edgecolor='none'))
    for temp in [300,800]:
        ax.plot([200,200],[temp,temp],[zlo,zhi],ls=(0,(3,3)),lw=.5,color='#A3ADB4',alpha=.7)
    ax.set(xlim=(0,700),ylim=(270,830),zlim=(zlo,zhi),
           xlabel='Time / ps',ylabel='Target temperature / K',zlabel='Volume change / %')
    ax.set_xticks([0,200,400,700])
    ax.set_yticks(TEMPS)
    ax.zaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.labelpad=9
    ax.yaxis.labelpad=10
    ax.zaxis.labelpad=7
    ax.tick_params(axis='x',labelsize=8,pad=1)
    ax.tick_params(axis='y',labelsize=8,pad=1)
    ax.tick_params(axis='z',labelsize=8,pad=1)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.pane.set_facecolor((1,1,1,0))
        axis.pane.set_edgecolor((1,1,1,0))
        axis._axinfo['grid'].update({'color':'#E8EDF0','linewidth':.4})
    ax.text2D(.14,.92,'Production starts: 200 ps',transform=ax.transAxes,
               fontsize=8.5,color='#687780')
    # Three matched diagnostic panels; all recorded values remain in view.
    lefts=[.095,.425,.755]
    lower=[fig.add_axes([left,.110,.225,.235]) for left in lefts]
    definitions=[('temperature_deviation_K','Temperature deviation / K'),
                 ('pressure_GPa','Pressure / GPa'),
                 ('potential_deviation_eV_atom','Potential-energy deviation\n/ eV atom$^{-1}$')]
    for axis,(key,label) in zip(lower,definitions):
        axis.set_axisbelow(True)
        axis.grid(axis='y',color='#E8EDF0',lw=.45)
        axis.axhline(0,color='#8B969E',lw=.65,zorder=1)
        axis.axvline(200,color='#8B969E',ls=(0,(3,3)),lw=.65,zorder=2)
        extrema=[]
        for temp,color,style in [(300,COLORS[0],'-'),(800,COLORS[-1],'--')]:
            frame=frames[temp]
            time=frame.time_ps.to_numpy()
            values=frame[key].to_numpy()
            axis.plot(time,values,color=color,alpha=.28,lw=.35,zorder=3)
            rolling=rolling_separately(time,values)
            # Draw separately to avoid a graphical bridge across the averaging reset.
            for mask in [time<=200,time>200]:
                axis.plot(time[mask],rolling[mask],color=color,ls=style,lw=.85,zorder=4)
            extrema.extend([values.min(),values.max()])
        lo,hi=min(extrema),max(extrema)
        padding=.08*(hi-lo)
        axis.set(xlim=(0,700),ylim=(lo-padding,hi+padding),xlabel='Time / ps',ylabel=label)
        axis.set_xticks([0,200,400,700])
        axis.yaxis.set_major_locator(MaxNLocator(4))
        axis.tick_params(direction='out',length=2.5,width=.6,pad=3)
        axis.xaxis.labelpad=6
        axis.yaxis.labelpad=5
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
    for letter,left in zip(['b','c','d'],lefts):
        fig.text(left-.065,.365,letter,fontsize=12,fontweight='bold')
    fig.text(.025,.952,'a',fontsize=12,fontweight='bold')
    handles=[Line2D([],[],color=COLORS[0],lw=1,ls='-',label='300 K'),
             Line2D([],[],color=COLORS[-1],lw=1,ls='--',label='800 K')]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.51,.402),
               frameon=False,ncol=2,fontsize=8.5,handlelength=2.4,columnspacing=2)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    essential=list(fig.texts)
    for axis in lower:
        essential.extend([axis.xaxis.label,axis.yaxis.label,*axis.get_xticklabels(),*axis.get_yticklabels()])
    # 3D tick-label transforms are finalized by the renderer; inspect them in the
    # exported proof rather than using their 2D get_window_extent coordinates.
    essential.extend([ax.xaxis.label,ax.yaxis.label,ax.zaxis.label,*ax.texts])
    clipping=[]
    for artist in essential:
        if not artist.get_visible() or not artist.get_text():
            continue
        bounds=artist.get_window_extent(renderer)
        if bounds.x0<0 or bounds.y0<0 or bounds.x1>fig.bbox.width or bounds.y1>fig.bbox.height:
            clipping.append(artist.get_text())
    assert not clipping,clipping
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f'{STEM}.pdf')
    fig.savefig(OUT/f'{STEM}.png',dpi=600)
    fig.savefig(OUT/f'{STEM}_preview.png',dpi=150)
    plt.close(fig)
    rows=[]
    for temp,frame in frames.items():
        row={'temperature_K':temp,'replicate':1,'reference_volume_A3':denominator,
             'recorded_rows':len(frame),'production_rows':5000,**means[temp]}
        for key in ['volume_change_percent','temperature_deviation_K','pressure_GPa','potential_deviation_eV_atom']:
            values=frame[key]
            production=values[frame.time_ps>200]
            row.update({f'{key}_full_min':float(values.min()),f'{key}_full_max':float(values.max()),
                        f'{key}_production_mean':float(production.mean()),
                        f'{key}_production_sd_descriptive':float(production.std(ddof=1))})
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT/f'{STEM}_summary.csv',index=False)
    caption='''**Figure 5. Recorded trajectory diagnostics for periodic single-crystal copper.** (a) Separate 3D volume trajectories for production repeat r1 at six target temperatures. Volume changes are normalized to the single common denominator V̄_300,r1 = '''+f'{denominator:.4f}'+''' Å³, averaged over 200 < time ≤ 700 ps. Target temperature labels separate simulation conditions rather than instantaneous temperatures. Gray traces indicate the first 200 ps; colored traces indicate the following production segment, with the recorded 200 ps boundary point retained for visual line continuity only. The faint plane marks 200 ps and is not an interpolated data surface. (b–d) Recorded endpoint r1 trajectories at 300 and 800 K: actual temperature minus target temperature, instantaneous pressure (1 bar = 10⁻⁴ GPa), and potential energy per atom minus that trajectory's own production mean. Light traces show every recorded point, including step 0; darker solid (300 K) and dashed (800 K) traces show trailing 50-record means (5 ps at the 0.1 ps recording interval), reset at the production boundary with shorter startup windows. These means are display aids only. No temporal decimation, detrending, variance rescaling, or temperature interpolation is used. All curves represent individual r1 trajectories, not averages over three repeats. Energy centering sets the production mean to zero by construction. The plots provide diagnostic evidence but do not alone prove complete equilibration or independent sampling; fluctuations are not confidence intervals.\n'''
    (ROOT/'docs/figure5_caption_v1.md').write_text(caption,encoding='utf8')
    for name,sha in protected.items():
        assert digest(ROOT/name)==sha
    for item in provenance:
        for name,sha in item['sha256'].items():
            assert digest(ROOT/name)==sha
    audit={'protected_inputs_sha256':protected,'sources':provenance,'reference_volume_A3':denominator,
        'volume_transform':'100*(V_T,r1(t)/mean(V_300,r1, 200<t<=700)-1)',
        'temperature_transform':'recorded temperature minus target temperature',
        'pressure_transform':'recorded bar * 1e-4 -> GPa',
        'energy_transform':'recorded total potential energy / 4000 minus own production mean',
        'display_rolling_mean':{'window_records':50,'nominal_window_ps':5,'alignment':'trailing',
                                'min_periods':1,'segments':'time<=200 and time>200 separately'},
        'decimation':'None; all 7001 recorded points per trajectory displayed',
        '3D_view':{'projection':'orthographic','elevation_deg':elevation,'azimuth_deg':azimuth,
                   'box_aspect':[3.8,2.,1.55],'zoom':1.08,'z_limits_percent':[zlo,zhi]},
        'font':font,'font_file':font_manager.findfont(font),'size_mm':[185,155],
        'png_dpi':600,'preview_dpi':150,'uncertainty':'No confidence intervals inferred from displayed traces',
        'notes':[] if font=='Times New Roman' else ['Times New Roman unavailable; STIXGeneral used.'],
        'checks':['Six completed r1 trajectories with 4000 atoms','Exact recorded steps and time range checked',
                  'Production means agree with run summary','All source and Notebook hashes unchanged',
                  'Figure text, 2D labels and 3D axis labels inside export canvas',
                  '3D tick-label placement requires exported visual review',
                  'Full recorded ranges retained in axes']}
    (OUT/f'{STEM}_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8')
    print(json.dumps({'font':font,'reference_volume_A3':denominator,'3D_view':audit['3D_view'],
                      'figure':str(OUT/f'{STEM}.png'),'notes':audit['notes']},indent=2))


if __name__=='__main__':
    main()
