"""Reproducible production and validation jobs. No fabricated/cached results."""
from __future__ import annotations
import argparse,concurrent.futures,hashlib,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OMP_WAIT_POLICY','PASSIVE')
from md_project import new_lammps,initialization,potential,save_json,CONFIG

def manifest():
    jobs=[]
    for rep,base in enumerate([13579,24681,97531],1):
        for t in [300,400,500,600,700,800]:
            jobs.append(dict(id=f'prod_T{t}_r{rep}',group='production',temperature=t,
                             seed=base+t*19,cells=10,dt=0.001,equil_ps=200,sample_ps=500,
                             replicate=rep,threads=4))
    for t in [300,800]:
        seed=13579+t*19
        for label,changes in [('size12',dict(cells=12)),('dt05',dict(dt=0.0005)),
                               ('long',dict(equil_ps=1000,sample_ps=1000))]:
            job=dict(id=f'val_{label}_T{t}',group='validation',temperature=t,seed=seed,
                     cells=10,dt=0.001,equil_ps=200,sample_ps=500,replicate=1,threads=4)
            job.update(changes);jobs.append(job)
    return jobs

def folder(job):return Path('outputs')/job['group']/job['id']

def source(job):
    n=job['cells'];dt=job['dt'];dest=folder(job).as_posix()
    if n==10:
        build='read_data outputs/model/relaxed.data\n'+potential()
    else:
        build=f'''lattice fcc {CONFIG['initial_lattice_A']}
region box block 0 {n} 0 {n} 0 {n} units lattice
create_box 1 box
create_atoms 1 box
'''+potential()+'''fix relax all box/relax iso 0 vmax 0.001
min_style cg
minimize 1e-12 1e-10 10000 100000
unfix relax
'''
    every=round(0.1/dt)
    return initialization()+build+f'''timestep {dt}
reset_timestep 0
velocity all create {job['temperature']} {job['seed']} mom yes rot no dist gaussian
fix integrator all npt temp {job['temperature']} {job['temperature']} 0.1 iso 0 0 1
variable clock equal time
variable temperature equal temp
variable pressure equal press
variable potential_energy equal pe
variable length_x equal lx
variable length_y equal ly
variable length_z equal lz
variable volume equal vol
thermo 1000
thermo_style custom step time atoms temp press pe etotal lx ly lz vol
thermo_modify format float %.12g flush yes
fix scalars all ave/time {every} 1 {every} v_clock v_temperature v_pressure v_potential_energy v_length_x v_length_y v_length_z v_volume file {dest}/thermo.dat format " %.12g"
dump frames all custom {round(100/dt)} {dest}/trajectory.lammpstrj id type x y z
dump_modify frames sort id
run {round(job['equil_ps']/dt)}
write_restart {dest}/equilibrated.restart
run {round(job['sample_ps']/dt)}
write_data {dest}/final.data
write_restart {dest}/final.restart
'''

def run_one(job):
    dest=folder(job);meta_file=dest/'metadata.json'
    if meta_file.exists():
        previous=json.loads(meta_file.read_text(encoding='utf8'))
        if previous.get('status')=='completed':return
        raise RuntimeError(f'Existing incomplete job protected: {dest}')
    dest.mkdir(parents=True,exist_ok=True)
    text=source(job);(dest/'input.in').write_text(text,encoding='ascii')
    meta={**job,'status':'running','potential_sha256':hashlib.sha256(Path(CONFIG['potential']).read_bytes()).hexdigest(),
          'input_sha256':hashlib.sha256(text.encode('ascii')).hexdigest(),
          'backend':'LAMMPS 22Jul2025 update4 / eam/alloy/omp double precision',
          'sample_interval_ps':0.1,'started_unix':time.time()}
    save_json(meta_file,meta);start=time.perf_counter()
    lmp=new_lammps((dest/'log.lammps').as_posix(),job['threads'])
    try:
        lmp.file((dest/'input.in').as_posix())
        assert int(lmp.get_natoms())==4*job['cells']**3
        meta.update(status='completed',final_atoms=int(lmp.get_natoms()),lammps_version=lmp.version())
    except BaseException as exc:
        meta.update(status='failed',error=str(exc));raise
    finally:
        lmp.close();meta['elapsed_seconds']=time.perf_counter()-start
        save_json(meta_file,meta)
    print('COMPLETED',job['id'],round(meta['elapsed_seconds'],1),flush=True)

def batch(workers):
    jobs=manifest();Path('outputs').mkdir(exist_ok=True)
    Path('.tmp').mkdir(exist_ok=True)
    save_json('outputs/study_manifest.json',jobs)
    # Long checks start early; each task is an isolated Python/LAMMPS process.
    jobs.sort(key=lambda j:(0 if '_long_' in j['id'] else 1,j['group'],j['id']))
    def launch(job):
        dest=folder(job)
        if (dest/'metadata.json').exists() and json.loads((dest/'metadata.json').read_text()).get('status')=='completed':
            return (job['id'],0)
        log=ROOT/'.tmp'/f"worker_{job['id']}.txt"
        with log.open('w',encoding='utf8') as stream:
            result=subprocess.run([sys.executable,str(Path(__file__).resolve()),'one','--job',job['id']],
                                  stdout=stream,stderr=subprocess.STDOUT,env={**os.environ,'PYTHONUTF8':'1','OMP_NUM_THREADS':'1'})
        return job['id'],result.returncode
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(launch,j) for j in jobs]
        for future in concurrent.futures.as_completed(futures):
            name,code=future.result();print('JOB',name,'exit',code,flush=True)
            if code:raise RuntimeError(f'Job failed: {name}; see .tmp worker log')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['batch','one','status'])
    parser.add_argument('--job');parser.add_argument('--workers',type=int,default=4);args=parser.parse_args()
    if args.action=='batch':batch(args.workers)
    elif args.action=='one':run_one(next(j for j in manifest() if j['id']==args.job))
    else:
        for j in manifest():
            p=folder(j)/'metadata.json'
            state=json.loads(p.read_text()).get('status') if p.exists() else 'pending'
            print(j['id'],state)
