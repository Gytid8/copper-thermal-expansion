"""Archive only the two confirmed interrupted jobs; prepare checkpoint recovery."""
import hashlib,json,shutil
from datetime import datetime
from pathlib import Path
from md_project import ROOT,save_json
from run_study import manifest,folder
from resume_study import continuation_source

archive_root=ROOT/'outputs'/'interrupted'/datetime.now().strftime('%Y%m%d_%H%M%S')
targets=[j for j in manifest() if j['id'] in ['val_dt05_T300','val_dt05_T800']]
for job in targets:
    dest=folder(job).resolve();archive=(archive_root/job['id']).resolve()
    assert dest.is_relative_to(ROOT/'outputs') and archive.is_relative_to(ROOT/'outputs')
    assert not archive.exists()
    old=json.loads((dest/'metadata.json').read_text(encoding='utf8'))
    assert old['status']=='running'
    checkpoint=dest/'equilibrated.restart';assert checkpoint.is_file()
    eqstep=round(job['equil_ps']/job['dt']);prefix=[];steps=[]
    for line in (dest/'thermo.dat').read_text(encoding='ascii').splitlines():
        parts=line.split()
        if line.startswith('#'):prefix.append(line)
        elif len(parts)==9 and parts[0].isdigit() and int(parts[0])<=eqstep:
            prefix.append(line);steps.append(int(parts[0]))
    assert steps==list(range(0,eqstep+1,round(.1/job['dt'])))
    # Validate continuation before making the reversible archive move.
    text=continuation_source(job)
    assert 'velocity all create' not in text and 'reset_timestep' not in text
    assert text.count('\nrun ')==1
    archive.parent.mkdir(parents=True,exist_ok=True)
    shutil.move(str(dest),str(archive))
    dest.mkdir()
    shutil.copy2(archive/'equilibrated.restart',dest/'equilibrated.restart')
    (dest/'equilibration_prefix.dat').write_text('\n'.join(prefix)+'\n',encoding='ascii')
    save_json(dest/'recovery.json',{
        'reason':'Worker process absent; scalar output stopped around 2026-10-17 20:17; no LAMMPS error recorded',
        'archive':archive.relative_to(ROOT).as_posix(),
        'checkpoint_step':eqstep,'checkpoint_time_ps':job['equil_ps'],
        'checkpoint_sha256':hashlib.sha256((dest/'equilibrated.restart').read_bytes()).hexdigest(),
        'discarded_production_note':'Interrupted production suffix is archived, not mixed into final averages; collect full 500 ps again',
        'original_metadata':old})
    print('PREPARED',job['id'],'checkpoint at',job['equil_ps'],'ps',flush=True)
