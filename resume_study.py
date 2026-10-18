"""Recover an interrupted trajectory from its saved equilibrium checkpoint.

Discard the interrupted production suffix; retain the equilibrium prefix and
collect a fresh, uninterrupted 500 ps production segment from the checkpoint.
"""
import hashlib,json,time
import shutil
from pathlib import Path
from md_project import initialization,potential,new_lammps,save_json
from run_study import folder,source

def validate_and_clean_boundary(job):
    """Allow only one numerically identical checkpoint-boundary duplicate."""
    import numpy as np
    dest=folder(job);path=dest/'thermo.dat'
    text=path.read_text(encoding='ascii');lines=text.splitlines()
    rows=[line.split() for line in lines if line.strip() and not line.startswith('#')]
    data=np.asarray(rows,dtype=float)
    assert data.ndim==2 and data.shape[1]==9 and np.isfinite(data).all()
    end=round((job['equil_ps']+job['sample_ps'])/job['dt'])
    interval=round(.1/job['dt']);expected=np.arange(0,end+1,interval)
    step=data[:,0];boundary=round(job['equil_ps']/job['dt'])
    unique,counts=np.unique(step,return_counts=True)
    assert np.array_equal(unique,expected), 'Missing or unexpected output steps'
    duplicates=unique[counts>1]
    if len(duplicates):
        assert np.array_equal(duplicates,[boundary]) and counts[unique==boundary][0]==2
        indices=np.flatnonzero(step==boundary)
        assert indices[1]==indices[0]+1
        assert np.allclose(data[indices[0]],data[indices[1]],rtol=1e-10,atol=1e-7), 'Boundary values disagree'
        raw=dest/'thermo_raw_before_boundary_cleanup.dat'
        if raw.exists():assert raw.read_text(encoding='ascii')==text
        else:shutil.copy2(path,raw)
        kept=[];seen=False
        for line in lines:
            parts=line.split()
            if len(parts)==9 and parts[0].isdigit() and int(parts[0])==boundary:
                if seen:continue
                seen=True
            kept.append(line)
        path.write_text('\n'.join(kept)+'\n',encoding='ascii')
    clean=np.loadtxt(path,comments='#')
    assert np.array_equal(clean[:,0],expected)
    assert np.isfinite(clean).all()
    assert np.allclose(clean[:,1],clean[:,0]*job['dt'],rtol=1e-10,atol=1e-8)
    assert np.allclose(clean[:,8],clean[:,5]*clean[:,6]*clean[:,7],rtol=1e-9)
    assert (clean[:,1]>job['equil_ps']).sum()==round(job['sample_ps']/.1)
    return {'rows':len(clean),'production_rows':round(job['sample_ps']/.1),
            'removed_boundary_duplicates':len(duplicates),
            'clean_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def continuation_source(job):
    dest=folder(job).as_posix();dt=job['dt']
    text=source(job);tail=text[text.index(f'timestep {dt}\n'):]
    tail='\n'.join(line for line in tail.split('\n')
                   if not line.startswith(('reset_timestep ','velocity ')))
    tail=tail.replace(f"run {round(job['equil_ps']/dt)}\nwrite_restart {dest}/equilibrated.restart\n",'')
    first=round(job['equil_ps']/dt)+round(.1/dt)
    tail=tail.replace(f'file {dest}/thermo.dat ',f'append {dest}/thermo.dat start {first} ')
    return initialization()+f'read_restart {dest}/equilibrated.restart\n'+potential()+tail

def resume_one(job):
    dest=folder(job);plan=json.loads((dest/'recovery.json').read_text(encoding='utf8'))
    meta_file=dest/'metadata.json'
    if meta_file.exists():
        previous=json.loads(meta_file.read_text(encoding='utf8'))
        if previous.get('status')=='completed':return
        raise RuntimeError('Recovery already attempted; inspect before restarting again')
    text=continuation_source(job)
    (dest/'input.in').write_text(text,encoding='ascii')
    (dest/'fresh_run.in').write_text(source(job),encoding='ascii')
    prefix=(dest/'equilibration_prefix.dat').read_text(encoding='ascii')
    (dest/'thermo.dat').write_text(prefix,encoding='ascii')
    meta={**job,'status':'running','started_unix':time.time(),
          'backend':'LAMMPS 22Jul2025 update4 / eam/alloy/omp double precision',
          'sample_interval_ps':.1,'recovery':plan,
          'potential_sha256':hashlib.sha256(Path('potentials/Cu01.eam.alloy').read_bytes()).hexdigest(),
          'input_sha256':hashlib.sha256(text.encode('ascii')).hexdigest()}
    save_json(meta_file,meta);start=time.perf_counter()
    lmp=new_lammps((dest/'log.lammps').as_posix(),job['threads'])
    try:
        lmp.file((dest/'input.in').as_posix())
        end=round((job['equil_ps']+job['sample_ps'])/job['dt'])
        assert int(lmp.get_thermo('step'))==end
        assert int(lmp.get_natoms())==4*job['cells']**3
        lmp.command('unfix scalars')
        meta['boundary_check']=validate_and_clean_boundary(job)
        meta.update(status='completed',final_atoms=int(lmp.get_natoms()),lammps_version=lmp.version())
    except BaseException as exc:
        meta.update(status='failed',error=str(exc));raise
    finally:
        lmp.close();meta['elapsed_seconds']=time.perf_counter()-start
        meta['elapsed_note']='Recovery segment only; interrupted wall time is recorded in recovery.json'
        save_json(meta_file,meta)
    print('COMPLETED RECOVERY',job['id'],meta['elapsed_seconds'],flush=True)
