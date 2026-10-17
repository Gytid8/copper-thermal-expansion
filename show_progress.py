"""Read-only live progress monitor; closing it does not stop simulations."""
import ctypes,json,os,time
from datetime import datetime
from pathlib import Path

ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
from run_study import manifest,folder

if os.name=='nt':
    kernel=ctypes.windll.kernel32
    handle=kernel.GetStdHandle(-11)
    mode=ctypes.c_ulong()
    if kernel.GetConsoleMode(handle,ctypes.byref(mode)):
        kernel.SetConsoleMode(handle,mode.value|4)

def last_step(path):
    try:
        with path.open('rb') as stream:
            stream.seek(max(0,path.stat().st_size-8192))
            lines=stream.read().decode('ascii',errors='replace').splitlines()
        for line in reversed(lines):
            parts=line.split()
            if len(parts)==9 and parts[0].isdigit():return int(parts[0])
    except OSError:pass
    return 0

def duration(seconds):
    if seconds is None:return '正在估算'
    minutes=max(0,round(seconds/60))
    return f'{minutes//60}小时{minutes%60:02d}分钟'

def snapshot():
    now=time.time();rows=[];rates=[];done=0;total=0;completed=0;weighted_remaining=0
    for j in manifest():
        count=round((j['equil_ps']+j['sample_ps'])/j['dt']);total+=count
        path=folder(j);meta={}
        try:meta=json.loads((path/'metadata.json').read_text(encoding='utf8'))
        except (OSError,json.JSONDecodeError):pass
        state=meta.get('status','pending');step=last_step(path/'thermo.dat') if meta else 0
        if state=='completed':step=count;completed+=1
        step=min(step,count);done+=step
        weighted_remaining+=(count-step)*(j['cells']/10)**3
        elapsed=now-meta.get('started_unix',now)
        restart_step=meta.get('recovery',{}).get('checkpoint_step',0)
        rate=max(0,step-restart_step)/elapsed if elapsed>5 and step else 0
        if state=='running' and rate:rates.append(rate*(j['cells']/10)**3)
        label={'pending':'等待中','running':'运行中','completed':'已完成','failed':'失败'}.get(state,state)
        remaining=(count-step)/rate if state=='running' and rate else None
        rows.append((j['id'],label,100*step/count,step*j['dt'],(j['equil_ps']+j['sample_ps']),remaining))
    percent=100*done/total;bar='#'*int(percent/2)+'-'*(50-int(percent/2))
    lines=['铜热膨胀模拟 · 实时进度（每5秒刷新）',datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
           '',f'总步数进度 [{bar}] {percent:5.1f}%',
           f'已完成任务：{completed}/24    当前运行：{sum(r[1]=="运行中" for r in rows)}',
           f'总步数：{done:,} / {total:,}',
           '预计计算剩余：'+duration(weighted_remaining/sum(rates) if rates else None)+'（粗估，尺寸检查可能更慢）',
           '',f'{"任务":<24} {"状态":<7} {"进度":>7}    模拟时间 / ps       当前任务剩余', '-'*94]
    for name,label,pct,ps,allps,remaining in rows:
        tail=duration(remaining) if label=='运行中' else ''
        lines.append(f'{name:<24} {label:<7} {pct:6.1f}%    {ps:7.1f}/{allps:<7.1f}   {tail}')
    lines+=['','prod = 正式计算；val_long = 时长验证；val_size12 = 尺寸验证；val_dt05 = 步长验证。',
            '任务完成后会继续统计、生成论文和上传。若检查失败，收尾流程会暂停。',
            '这是只读查看窗口：关闭窗口或按 Ctrl+C 仅停止显示，不会停止模拟。']
    err=ROOT/'.tmp/finalize_stderr.txt'
    if err.exists() and err.stat().st_size:
        lines+=['','收尾流程有提示，请查看 .tmp/finalize_stderr.txt：',err.read_text(encoding='utf8',errors='replace')[-1200:]]
    return '\n'.join(lines)

if __name__=='__main__':
    try:
        while True:
            print('\033[2J\033[H'+snapshot(),flush=True)
            time.sleep(5)
    except KeyboardInterrupt:pass
