"""Verify finished recovery outputs and repair only a duplicate boundary record."""
import json,shutil,time
from md_project import new_lammps,initialization,potential,save_json
from run_study import manifest,folder
from resume_study import validate_and_clean_boundary

for job in manifest():
    if job['id'] not in ['val_dt05_T300','val_dt05_T800']:continue
    dest=folder(job);file=dest/'metadata.json'
    meta=json.loads(file.read_text(encoding='utf8'))
    assert meta['status']=='failed' and meta.get('recovery')
    log=(dest/'log.lammps').read_text(errors='replace')
    assert 'ERROR' not in log and 'for 1000000 steps with 4000 atoms' in log
    assert 'write_restart '+dest.as_posix()+'/final.restart' in log
    check=validate_and_clean_boundary(job)
    lmp=new_lammps('.tmp/final_checkpoint_check_'+str(job['temperature'])+'.log',4)
    try:
        lmp.commands_string(initialization()+'read_restart '+(dest/'final.restart').as_posix()+'\n'+potential()+
            f"fix integrator all npt temp {job['temperature']} {job['temperature']} 0.1 iso 0 0 1\nrun 0\n")
        assert int(lmp.get_thermo('step'))==1400000
        assert abs(lmp.get_thermo('time')-700)<1e-8
        assert lmp.get_natoms()==4000
        version=lmp.version()
    finally:lmp.close()
    backup=dest/'metadata_before_boundary_cleanup.json'
    assert not backup.exists();shutil.copy2(file,backup)
    meta.update(status='completed',final_atoms=4000,lammps_version=version,boundary_check=check,
                repaired_unix=time.time(),repair_note='MD completed; only duplicated 200 ps checkpoint row removed after strict validation; original raw data and metadata preserved')
    meta.pop('error',None)
    save_json(file,meta)
    print('VERIFIED COMPLETED',job['id'],check,flush=True)
