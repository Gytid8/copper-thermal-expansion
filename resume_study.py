"""Recover an interrupted trajectory from its saved equilibrium checkpoint.

Discard the interrupted production suffix; retain the equilibrium prefix and
collect a fresh, uninterrupted 500 ps production segment from the checkpoint.
"""
import hashlib,json,time
from pathlib import Path
from md_project import initialization,potential,new_lammps,save_json
from run_study import folder,source

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
        steps=[]
        for line in (dest/'thermo.dat').read_text(encoding='ascii').splitlines():
            parts=line.split()
            if len(parts)==9 and parts[0].isdigit():steps.append(int(parts[0]))
        assert steps==list(range(0,end+1,round(.1/job['dt'])))
        meta.update(status='completed',final_atoms=int(lmp.get_natoms()),lammps_version=lmp.version())
    except BaseException as exc:
        meta.update(status='failed',error=str(exc));raise
    finally:
        lmp.close();meta['elapsed_seconds']=time.perf_counter()-start
        meta['elapsed_note']='Recovery segment only; interrupted wall time is recorded in recovery.json'
        save_json(meta_file,meta)
    print('COMPLETED RECOVERY',job['id'],meta['elapsed_seconds'],flush=True)
