"""Wait for the authorized study, then verify/report/publish. No fake completion."""
from pathlib import Path
import argparse,json,subprocess,sys,time
from datetime import datetime,timezone
from run_study import manifest,folder
ROOT=Path(__file__).resolve().parent

def command(args):
    subprocess.run(args,cwd=ROOT,check=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--publish',action='store_true');a=p.parse_args()
    last=None
    while True:
        statuses=[]
        for j in manifest():
            path=folder(j)/'metadata.json'
            try:
                statuses.append(json.loads(path.read_text(encoding='utf8')).get('status') if path.exists() else 'pending')
            except (json.JSONDecodeError,PermissionError):
                statuses.append('updating')
        if 'failed' in statuses:
            raise RuntimeError('A simulation failed. Inspect metadata/log; no final report published.')
        complete=statuses.count('completed')
        if complete!=last:
            print(datetime.now(timezone.utc).isoformat(),f'completed {complete}/24',flush=True)
            command([sys.executable,'write_readme.py']);last=complete
        if complete==24:break
        time.sleep(30)
    for script in ['check_methods.py','analyze_study.py','write_paper.py','execute_notebook.py','write_readme.py']:
        command([sys.executable,script])
    m=json.loads((ROOT/'results/metrics.json').read_text(encoding='utf8'))
    if not m['all_practical_checks_pass']:
        raise RuntimeError('Validation exceeds preset 2% threshold. Results retained for review, no final push.')
    if m['max_halves_z']>4:
        raise RuntimeError('Strong half-window drift diagnostic: manual scientific review required before final push.')
    if a.publish:
        command(['git','add','.'])
        command(['git','commit','-m','Complete copper thermal expansion study and executed Notebook'])
        command(['git','push','origin','main'])
        print('Remote main:',flush=True);command(['git','ls-remote','origin','refs/heads/main'])
    print('Completed report pipeline',flush=True)

if __name__=='__main__':main()
