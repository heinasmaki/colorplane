"""Measure CPU scaling (and optionally GPU), retaining energies and provenance."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import time


def benchmark(exe,output,threads=(1,2,4),repeats=2,size=192,spacing=.05,steps=3,sweeps=10,gpu=False):
    exe,output=Path(exe).resolve(),Path(output).resolve()
    if not exe.is_file():
        raise ValueError('Executable does not exist')
    if min(repeats,size,steps,sweeps)<1 or size<2 or not threads or any(t<1 for t in threads):
        raise ValueError('Invalid benchmark dimensions or thread counts')
    output.mkdir(parents=True,exist_ok=False)
    manifest=dict(executable=str(exe),executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),
                  platform=platform.platform(),logical_cpus=os.cpu_count(),threads=list(threads),
                  repeats=repeats,size=size,spacing=spacing,steps=steps,sweeps=sweeps,gpu=gpu,
                  environment={k:v for k,v in os.environ.items() if k.startswith('OMP_')})
    buildfile=exe.with_suffix('.build.json')
    if buildfile.exists():
        manifest['build']=json.loads(buildfile.read_text(encoding='utf-8'))
        if manifest['build']['executable_sha256']!=manifest['executable_sha256']:
            raise ValueError('Executable and build manifest disagree')
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    configurations=[('blocked',None)]+[('cpu',t) for t in dict.fromkeys(threads)]
    if gpu:
        configurations.append(('gpu',None))
    results=[]
    for repeat in range(repeats):
        reference_hash=None
        for backend,count in configurations:
            label=backend if count is None else f'cpu-{count}'
            directory=output/f'repeat-{repeat+1}'/label
            directory.mkdir(parents=True)
            text='\n'.join(map(str,[3,size,spacing,steps,.9,sweeps,'run','n','n']))+'\n'
            (directory/'input.inp').write_text(text,encoding='ascii')
            command=[str(exe),'--backend',backend,'--seed',str(123+repeat)]
            if count is not None:
                command+=['--threads',str(count)]
            (directory/'command.json').write_text(json.dumps(command),encoding='utf-8')
            start=time.perf_counter()
            proc=subprocess.run(command,input=text,text=True,capture_output=True,cwd=directory)
            elapsed=time.perf_counter()-start
            (directory/'console.log').write_text(proc.stdout+proc.stderr,encoding='utf-8')
            if proc.returncode:
                raise RuntimeError(f'{label} failed; see {directory / "console.log"}')
            final_hash=hashlib.sha256((directory/'run.lat').read_bytes()).hexdigest()
            if backend=='blocked':
                reference_hash=final_hash
            if backend=='cpu' and final_hash!=reference_hash:
                raise RuntimeError('CPU thread count changed the trajectory; benchmark invalid')
            with (directory/'run.metrics.csv').open() as source:
                final=list(csv.DictReader(source))[-1]
            timing=re.search(r'Annealing wall seconds:\s*([0-9.]+)',proc.stdout)
            row=dict(repeat=repeat+1,label=label,backend=backend,threads=count,
                     process_wall_seconds=elapsed,annealing_wall_seconds=float(timing.group(1)),
                     final_energy=float(final['energy']),matches_blocked_lattice=final_hash==reference_hash)
            results.append(row)
            print(f'{label}: {row["annealing_wall_seconds"]:.6f} s annealing, '
                  f'{elapsed:.3f} s process, energy {row["final_energy"]:.9g}',flush=True)
            with (output/'results.csv').open('w',newline='',encoding='utf-8') as target:
                writer=csv.DictWriter(target,fieldnames=list(row))
                writer.writeheader()
                writer.writerows(results)
    baseline=statistics.median(r['annealing_wall_seconds'] for r in results if r['label']=='blocked')
    summary=[]
    for label in dict.fromkeys(r['label'] for r in results):
        duration=statistics.median(r['annealing_wall_seconds'] for r in results if r['label']==label)
        summary.append(dict(label=label,median_annealing_seconds=duration,
                            speedup_vs_blocked=baseline/duration if duration else None))
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return results


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exe',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--threads',nargs='+',type=int,default=[1,2,4])
    parser.add_argument('--repeats',type=int,default=2)
    parser.add_argument('--size',type=int,default=192)
    parser.add_argument('--spacing',type=float,default=.05)
    parser.add_argument('--steps',type=int,default=3)
    parser.add_argument('--sweeps',type=int,default=10)
    parser.add_argument('--gpu',action='store_true')
    args=parser.parse_args()
    benchmark(args.exe,args.output,args.threads,args.repeats,args.size,args.spacing,args.steps,args.sweeps,args.gpu)
