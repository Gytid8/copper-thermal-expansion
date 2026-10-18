"""Analyze completed MD files; fail rather than invent missing results."""
from __future__ import annotations
import json,os
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
from md_project import plotting,save_json
from run_study import manifest,folder
COLS=['step','time_ps','temperature_K','pressure_bar','potential_eV','lx_A','ly_A','lz_A','volume_A3']
RESULT=Path('results'); FIG=RESULT/'figures'

def read_job(job):
    dest=folder(job);meta=json.loads((dest/'metadata.json').read_text())
    if meta['status']!='completed':raise RuntimeError(f'Incomplete: {job["id"]}')
    frame=pd.read_csv(dest/'thermo.dat',sep=r'\s+',comment='#',names=COLS)
    steps=round((job['equil_ps']+job['sample_ps'])/job['dt'])
    interval=round(.1/job['dt'])
    assert np.array_equal(frame.step,np.arange(0,steps+1,interval))
    assert np.isfinite(frame.to_numpy()).all()
    assert np.allclose(frame.lx_A,frame.ly_A,rtol=1e-9)
    assert np.allclose(frame.lx_A,frame.lz_A,rtol=1e-9)
    assert np.allclose(frame.volume_A3,frame.lx_A*frame.ly_A*frame.lz_A,rtol=1e-9)
    assert meta['final_atoms']==4*job['cells']**3
    log=(dest/'log.lammps').read_text(errors='replace')
    assert 'ERROR' not in log
    sample=frame[frame.time_ps>job['equil_ps']].copy()
    assert len(sample)==round(job['sample_ps']/.1)
    sample['a_A']=sample.lx_A/job['cells']
    values=sample.a_A.to_numpy();half=len(values)//2
    row={**job,'atoms':4*job['cells']**3,'samples':len(sample),
         'a_A':values.mean(),'volume_A3':sample.volume_A3.mean(),
         'T_mean_K':sample.temperature_K.mean(),'P_mean_bar':sample.pressure_bar.mean(),
         'first_half_a_A':values[:half].mean(),'second_half_a_A':values[half:].mean(),
         'halves_delta_A':values[half:].mean()-values[:half].mean(),
         'wall_seconds':meta['elapsed_seconds']}
    blocks=[]
    for length in [1,2,5,10,20,50]:
        n=round(length/.1);count=len(sample)//n
        means=values[:n*count].reshape(count,n).mean(axis=1)
        sem=means.std(ddof=1)/np.sqrt(count)
        blocks.append({'job':job['id'],'block_ps':length,'blocks':count,'a_sem_A':sem})
        if length==20:
            row['a_block_sem_A']=sem
            # Match the two physical half-windows; discard each half's short tail.
            left=values[:half];right=values[half:]
            kl=len(left)//n;kr=len(right)//n
            ml=left[:kl*n].reshape(kl,n).mean(axis=1)
            mr=right[:kr*n].reshape(kr,n).mean(axis=1)
            row['halves_delta_sem_A']=np.sqrt(ml.var(ddof=1)/kl+mr.var(ddof=1)/kr)
            for key in ['temperature_K','pressure_bar']:
                m=sample[key].to_numpy()[:n*count].reshape(count,n).mean(axis=1)
                row[key+'_sem']=m.std(ddof=1)/np.sqrt(count)
    # Window-to-window agreement is a diagnostic, not a proof of equilibrium.
    row['halves_z']=abs(row['halves_delta_A'])/row['halves_delta_sem_A']
    row['temperature_offset_K']=row['T_mean_K']-job['temperature']
    return row,blocks,frame

def reference():
    data=pd.read_csv('data/reference/nist_srm736.csv')
    e=data.expansion_from_293_ppm.to_numpy()*1e-6
    e300=float(data.loc[data.temperature_K==300,'expansion_from_293_ppm'].iloc[0])*1e-6
    data['strain_from_300']=(1+e)/(1+e300)-1
    data['alpha_fixed300_Kinv']=data.expansivity_fixed293_micro_Kinv*1e-6/(1+e300)
    mean=float(data.loc[data.temperature_K==800,'strain_from_300'].iloc[0])/500
    return data,mean

def main():
    RESULT.mkdir(exist_ok=True);FIG.mkdir(exist_ok=True)
    rows=[];blockrows=[];frames={}
    for job in manifest():
        row,blocks,frame=read_job(job);rows.append(row);blockrows+=blocks;frames[job['id']]=frame
    runs=pd.DataFrame(rows);runs.to_csv(RESULT/'run_summary.csv',index=False)
    pd.DataFrame(blockrows).to_csv(RESULT/'block_statistics.csv',index=False)
    prod=runs[runs.group=='production'].copy()
    summary=prod.groupby('temperature').agg(a_A=('a_A','mean'),a_sd_A=('a_A','std'),
        volume_A3=('volume_A3','mean'),volume_sd_A3=('volume_A3','std'),
        temperature_K=('T_mean_K','mean'),pressure_bar=('P_mean_bar','mean'),replicates=('id','count')).reset_index()
    assert (summary.replicates==3).all()
    summary['a_sem_A']=summary.a_sd_A/np.sqrt(3)
    summary['a_ci95_A']=4.3026527299*summary.a_sem_A
    a300=float(summary.loc[summary.temperature==300,'a_A'].iloc[0]);v300=float(summary.loc[summary.temperature==300,'volume_A3'].iloc[0])
    summary['strain_from_300']=summary.a_A/a300-1
    summary['relative_volume']=summary.volume_A3/v300-1
    summary.to_csv(RESULT/'temperature_summary.csv',index=False)
    coefficients=[]
    for rep in [1,2,3]:
        sub=prod[prod.replicate==rep].set_index('temperature')
        a0=float(sub.loc[300,'a_A']);a1=float(sub.loc[800,'a_A'])
        coefficients.append({'replicate':rep,'a300_A':a0,'a800_A':a1,'alpha_Kinv':(a1/a0-1)/500})
    coefs=pd.DataFrame(coefficients);coefs.to_csv(RESULT/'alpha_replicates.csv',index=False)
    alpha=(float(summary.loc[summary.temperature==800,'a_A'].iloc[0])/a300-1)/500
    alpha_sd=coefs.alpha_Kinv.std(ddof=1);alpha_sem=alpha_sd/np.sqrt(3)
    ref,refalpha=reference();ref.to_csv(RESULT/'reference_rebased.csv',index=False)
    comparisons=[]
    for kind in ['size12','dt05','long']:
        sub=runs[runs.id.str.contains('val_'+kind+'_')].set_index('temperature')
        av=(sub.loc[800,'a_A']/sub.loc[300,'a_A']-1)/500
        base=prod[prod.replicate==1].set_index('temperature')
        av0=(base.loc[800,'a_A']/base.loc[300,'a_A']-1)/500
        row={'check':kind,'alpha_Kinv':av,'baseline_r1_alpha_Kinv':av0,
             'change_percent':(av/av0-1)*100,'practical_limit_percent':2.0,
             'within_2_percent':abs(av/av0-1)<=.02}
        for t in [300,800]:
            row[f'a_delta_{t}_A']=sub.loc[t,'a_A']-base.loc[t,'a_A']
            row[f'a_combined_sem_{t}_A']=np.hypot(sub.loc[t,'a_block_sem_A'],base.loc[t,'a_block_sem_A'])
        comparisons.append(row)
    checks=pd.DataFrame(comparisons);checks.to_csv(RESULT/'validation_comparison.csv',index=False)
    metrics={'alpha_300_800_Kinv':alpha,'alpha_replicate_sd_Kinv':alpha_sd,
             'alpha_replicate_sem_Kinv':alpha_sem,'alpha_ci95_halfwidth_Kinv':4.3026527299*alpha_sem,
             'nist_alpha_300_800_Kinv':refalpha,'relative_error_percent':(alpha/refalpha-1)*100,
             'production_runs':18,'validation_runs':6,'all_practical_checks_pass':bool(checks.within_2_percent.all()),
             'max_halves_z':float(runs.halves_z.max()),'confidence_note':'t distribution, df=2, 3 independent repeats; statistical uncertainty only',
             'reference_note':'NIST archived SRM 736 certificate, revision 1990-10-07, Table page 1; 293 K rebased to 300 K'}
    save_json(RESULT/'metrics.json',metrics)
    plot_all(summary,prod,runs,frames,ref,checks,metrics,pd.DataFrame(blockrows))
    print(json.dumps(metrics,ensure_ascii=False,indent=2))

def plot_all(summary,prod,runs,frames,ref,checks,metrics,blocks):
    plt=plotting();color='#b8662c';blue='#226b8e';x=summary.temperature.to_numpy()
    fig,ax=plt.subplots(figsize=(8,5))
    ax.errorbar(x,summary.a_A,yerr=summary.a_ci95_A,fmt='o-',color=color,capsize=4,label='MD：3次重复均值，95%统计区间')
    ax.set(xlabel='温度 (K)',ylabel='晶格常数 (Å)',title='单晶铜晶格常数随温度变化');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(FIG/'lattice_temperature.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    axes[0].plot(x,100*summary.strain_from_300,'o-',color=color,label='Mishin EAM1 / MD')
    axes[0].plot(ref.temperature_K,100*ref.strain_from_300,'s--',color=blue,label='NIST SRM 736（1990修订）')
    axes[0].set(xlabel='温度 (K)',ylabel='相对300 K的线膨胀 (%)',title='同参考温度的实验对照');axes[0].legend(fontsize=8)
    axes[1].plot(x,100*summary.relative_volume,'o-',color=color)
    axes[1].set(xlabel='温度 (K)',ylabel='相对300 K的体积变化 (%)',title='体积热膨胀')
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(FIG/'expansion_comparison.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,5))
    ax.bar(['MD温区平均','NIST温区平均'],[metrics['alpha_300_800_Kinv']*1e6,metrics['nist_alpha_300_800_Kinv']*1e6],color=[color,blue])
    ax.errorbar(0,metrics['alpha_300_800_Kinv']*1e6,yerr=metrics['alpha_ci95_halfwidth_Kinv']*1e6,fmt='none',ecolor='black',capsize=5)
    ax.set(ylabel=r'300–800 K平均线膨胀系数 ($10^{-6}$ K$^{-1}$)',title='工程温区平均；误差棒仅含MD重复统计误差')
    fig.tight_layout();fig.savefig(FIG/'alpha_comparison.png');plt.close(fig)
    for t in [300,800]:
        jobid=f'prod_T{t}_r1';frame=frames[jobid];jobrow=runs[runs.id==jobid].iloc[0]
        fig,axes=plt.subplots(4,1,figsize=(10,10),sharex=True)
        for ax,key,label,scale in zip(axes,['temperature_K','pressure_bar','potential_eV','volume_A3'],
            ['温度 (K)','压力 (bar)','势能 (eV/atom)','体积 (Å³)'],[1,1,1/4000,1]):
            ax.plot(frame.time_ps,frame[key]*scale,lw=.4,color='#aec3ce')
            ax.plot(frame.time_ps,(frame[key]*scale).rolling(50,min_periods=1).mean(),lw=1,color=color)
            ax.axvline(jobrow.equil_ps,ls='--',color='black',lw=.8);ax.set_ylabel(label);ax.grid(alpha=.2)
        axes[-1].set_xlabel('时间 (ps)；虚线后采样，深线为5 ps展示平滑')
        fig.suptitle(f'{t} K正式计算 r1：200 ps平衡＋500 ps采样')
        fig.tight_layout();fig.savefig(FIG/f'equilibration_{t}.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.6))
    for ax,t in zip(axes,[300,800]):
        sub=blocks[blocks.job==f'prod_T{t}_r1'];ax.plot(sub.block_ps,sub.a_sem_A*1e6,'o-',color=color)
        ax.set(xlabel='分块长度 (ps)',ylabel=r'晶格常数均值标准误 ($10^{-6}$ Å)',title=f'{t} K分块敏感性');ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(FIG/'block_sensitivity.png');plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.8))
    ax.axhline(0,color='gray',lw=.8);ax.axhspan(-2,2,color='#9ec5b1',alpha=.25,label='预设实用阈值±2%')
    ax.bar(['尺寸10→12晶胞','时间步1→0.5 fs','时长延长至1＋1 ns'],checks.change_percent,color=color)
    ax.set(ylabel='平均膨胀系数相对基准r1的变化 (%)',title='端点验证：每项各300 K和800 K一条轨迹');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(FIG/'validation.png');plt.close(fig)

if __name__=='__main__':main()
