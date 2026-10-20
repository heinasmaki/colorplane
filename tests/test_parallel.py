"""Race-safety, independent Monte Carlo reference, replay, and optional GPU checks."""
import argparse
import csv
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from build import build,compiler_config
from benchmark_parallel import benchmark
from ensemble import run_ensemble

EXE=PROBE=None
REQUIRE_GPU=False
MASK=(1<<32)-1


def philox(counter,key):
    """Independent Python integer arithmetic, with native unbounded products."""
    c=list(counter)
    k=list(key)
    for _ in range(10):
        p0=0xD2511F53*c[0]
        p1=0xCD9E8D57*c[2]
        c=[(p1>>32)^c[1]^k[0],p1&MASK,(p0>>32)^c[3]^k[1],p0&MASK]
        k=[(k[0]+0x9E3779B9)&MASK,(k[1]+0xBB67AE85)&MASK]
    return c


def random_move(seed,site,sweep,q):
    trial=0
    limit=((1<<32)//(q-1))*(q-1)
    while True:
        r=philox([site,sweep&MASK,sweep>>32,trial],[seed&MASK,seed>>32])
        if r[0]<limit:
            return 1+r[0]%(q-1),((r[2]<<21)+(r[3]>>11))*2**-53
        trial+=1


def stencil(spacing,name):
    d=math.ceil(2+1/spacing)
    result=[]
    for dx in range(-d,d+1):
        for dy in range(-d,d+1):
            r=math.hypot(dx,dy)-1/spacing
            phi=max(0,1-abs(r)) if name=='hat' else (0.25*(1+math.cos(math.pi*r/2)) if abs(r)<2 else 0)
            w=phi/spacing
            if w>1e-14:
                result.append((dx,dy,w))
    return result


def reference_sweeps(q,nx,ny,a,name,seed,t,loops):
    offsets=stencil(a,name)
    tile=min(max(abs(dx),abs(dy)) for dx,dy,_ in offsets)
    radius=max(max(abs(dx),abs(dy)) for dx,dy,_ in offsets)
    period=math.ceil(radius/tile)+1
    lattice=[[1+(3*(x+1)+7*(y+1)+(x+1)*(y+1))%q for x in range(nx)] for y in range(ny)]
    groups=[[] for _ in range(period**2)]
    for x in range(nx):
        for y in range(ny):
            groups[(x//tile)%period+period*((y//tile)%period)].append((x,y))
    accepted=changed=0
    for sweep in range(1,loops+1):
        for group in groups:
            # This is a sequential implementation; independent updates commute.
            for x,y in group:
                old=lattice[y][x]
                candidate,u=random_move(seed,1+y+x*ny,sweep,q)
                if candidate>=old:
                    candidate+=1
                eold=enew=0.0
                for dx,dy,w in offsets:
                    if 0<=x+dx<nx and 0<=y+dy<ny:
                        neighbor=lattice[y+dy][x+dx]
                        if neighbor==old:
                            eold+=w
                        if neighbor==candidate:
                            enew+=w
                delta=enew-eold
                if delta<=0 or u<math.exp(-delta/t):
                    lattice[y][x]=candidate
                    accepted+=1
                    changed+=delta!=0
    return lattice,accepted,changed


def reference_map(lattice,a,name):
    ny,nx=len(lattice),len(lattice[0])
    offsets=stencil(a,name)
    norm=math.fsum(w for _,_,w in offsets)
    result=[]
    for y in range(ny):
        row=[]
        for x in range(nx):
            row.append(math.fsum(w for dx,dy,w in offsets if 0<=x+dx<nx and 0<=y+dy<ny
                                 and lattice[y][x]==lattice[y+dy][x+dx])/norm)
        result.append(row)
    return result


class ParallelTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='colorplane-parallel-')
        self.directory=Path(self.temp.name)
        self.env=dict(os.environ,OMP_NUM_THREADS='4',OMP_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1')

    def tearDown(self):
        self.temp.cleanup()

    def probe(self,mode,text,backend=None):
        proc=subprocess.run([str(PROBE),mode,*([backend] if backend else [])],input=text,text=True,
                            capture_output=True,env=self.env,timeout=120)
        self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
        return proc.stdout.splitlines()

    def run_model(self,root,backend,steps=3,sweeps=5,offset=0,seed=23,order='typewriter',threads=None,
                  kernel='cosine',ok=True):
        text='\n'.join(map(str,[3,31,.1,steps,.9,sweeps,root,'y','y',37,offset]))+'\n'
        command=[str(EXE),'--backend',backend,'--order',order,'--kernel',kernel,
                 '--reheat-every','2','--reheat-factor','2']
        if seed is not None:
            command+=['--seed',str(seed)]
        if threads is not None:
            command+=['--threads',str(threads)]
        proc=subprocess.run(command,input=text,text=True,cwd=self.directory,capture_output=True,
                            env=self.env,timeout=120)
        self.assertEqual(proc.returncode==0,ok,proc.stdout+proc.stderr)
        return proc

    def rows(self,root):
        with (self.directory/(root+'.metrics.csv')).open() as source:
            return list(csv.DictReader(source))

    def test_published_rng_vectors_and_integer_overflow_safety(self):
        # Random123 tests/kat_vectors: Philox4x32-10, zero/all-one/pi inputs.
        vectors=[([0]*4,[0]*2,['6627e8d5','e169c58d','bc57ac4c','9b00dbd8']),
                 ([MASK]*4,[MASK]*2,['408f276d','41c83b0e','a20bc7c6','6d5451fd']),
                 ([0x243f6a88,0x85a308d3,0x13198a2e,0x03707344],[0xa4093822,0x299f31d0],
                  ['d16cfe09','94fdcceb','5001e420','24126ea1'])]
        rng=random.Random(197)
        for _ in range(50):
            c=[rng.getrandbits(32) for _ in range(4)]
            k=[rng.getrandbits(32) for _ in range(2)]
            vectors.append((c,k,[f'{v:08x}' for v in philox(c,k)]))
        lines=self.probe('rng','\n'.join(' '.join(map(str,c+k)) for c,k,_ in vectors)+'\n')
        self.assertEqual(len(lines),len(vectors))
        for line,(_,_,expected) in zip(lines,vectors):
            self.assertEqual(line.lower().split(),expected)

    def test_move_stream_uses_site_seed_and_full_64_bit_sweep(self):
        inputs=[(s,site,sweep,q) for s in (0,2147483646) for site in (1,10001)
                # Large q is a unit-test input to force rejection attempts in
                # the generic RNG helper; physical models still allow 2..255.
                for sweep in (1,2**32+5) for q in (2,3,7,255,1500000001)]
        lines=self.probe('moves','\n'.join(' '.join(map(str,row)) for row in inputs)+'\n')
        for line,args in zip(lines,inputs):
            candidate,u=line.split()
            expected,expected_u=random_move(*args)
            self.assertEqual(int(candidate),expected)
            self.assertEqual(float(u),expected_u)

    def test_coloring_has_no_interaction_edges_inside_a_batch(self):
        for a in (.49,.3,.2,.1,.05):
            for name in ('cosine','hat'):
                nx,ny=math.ceil(4/a)+7,math.ceil(3/a)+5
                lines=self.probe('plan',f'3 {nx} {ny} {a}\n{name}\n')
                tile,period,groups,padding=map(int,lines[0].split())
                self.assertEqual(groups,period**2)
                colors={}
                for line in lines[1:]:
                    site,group,center=map(int,line.split())
                    self.assertNotIn(site,colors)
                    colors[site]=group
                    y=(site-1)%ny
                    x=(site-1)//ny
                    self.assertEqual(center,1+y+padding+(x+padding)*(ny+2*padding))
                self.assertEqual(set(colors),set(range(1,nx*ny+1)))
                offsets=stencil(a,name)
                for x in range(nx):
                    for y in range(ny):
                        color=colors[1+y+x*ny]
                        for dx,dy,_ in offsets:
                            if 0<=x+dx<nx and 0<=y+dy<ny:
                                self.assertNotEqual(color,colors[1+y+dy+(x+dx)*ny])

    def test_complete_sweeps_match_independent_python_reference(self):
        for name in ('cosine','hat'):
            for backend in ('blocked','cpu','target-host'):
                for q,a,t in ((3,.125,1.7),(2,.2,.8),(255,.1,.3),(3,.49,1.7)):
                    nx,ny,seed,loops=21,17,125,4
                    lines=self.probe('sweep',f'{q} {nx} {ny} {a}\n{name}\n{seed} {t} {loops}\n',backend)
                    energy,accepted,changed=lines[0].split()
                    actual=[list(map(int,row.split())) for row in lines[1:]]
                    expected,ea,ec=reference_sweeps(q,nx,ny,a,name,seed,t,loops)
                    self.assertEqual(actual,expected)
                    self.assertEqual((int(accepted),int(changed)),(ea,ec))
                    energy_map=reference_map(actual,a,name)
                    self.assertAlmostEqual(float(energy),math.fsum(sum(energy_map,[]))/(nx*ny),places=13)

    def test_cold_sweeps_never_increase_energy(self):
        values=[]
        for loops in (1,2,3,4):
            lines=self.probe('sweep',f'3 21 17 0.125\ncosine\n125 1e-30 {loops}\n','cpu')
            values.append(float(lines[0].split()[0]))
        self.assertTrue(all(after<=before+1e-14 for before,after in zip(values,values[1:])))

    def test_cpu_thread_counts_and_target_host_replay_exactly(self):
        self.run_model('reference','blocked',order='random')
        for backend,threads in (('cpu',1),('cpu',2),('cpu',4),('target-host',None)):
            root=f'{backend}-{threads}'
            self.run_model(root,backend,threads=threads,order='random')
            for ext in ('.lat','.ngb'):
                self.assertEqual((self.directory/(root+ext)).read_bytes(),
                                 (self.directory/('reference'+ext)).read_bytes())
            self.assertEqual(self.rows(root),self.rows('reference'))

    def test_exact_restart_and_cpu_thread_count_change(self):
        for order in ('typewriter','random'):
            self.run_model(order+'full','cpu',steps=5,sweeps=9,offset=1,threads=4,order=order)
            self.run_model(order+'split','blocked',steps=2,sweeps=9,offset=1,order=order)
            self.run_model(order+'split','cpu',steps=3,sweeps=9,offset=1,threads=2,seed=None,order=order)
            self.assertEqual((self.directory/(order+'full.lat')).read_bytes(),
                             (self.directory/(order+'split-res.lat')).read_bytes())
            self.assertEqual(self.rows(order+'full')[2:],self.rows(order+'split-res'))
        # A legacy trajectory cannot be resumed as a blocked trajectory.
        self.run_model('legacy','serial')
        p=self.run_model('legacy','cpu',seed=None,ok=False)
        self.assertIn('algorithm mismatch',p.stderr)

    def test_parallel_map_matches_independent_final_energy(self):
        self.run_model('cpu','cpu',threads=3)
        rows=self.rows('cpu')
        local=[float(v) for line in (self.directory/'cpu.ngb').read_text().splitlines() for v in line.split()]
        lattice=[list(map(int,line.split())) for line in (self.directory/'cpu.lat').read_text().splitlines()]
        expected=sum(reference_map(lattice,.1,'cosine'),[])
        for actual,wanted in zip(local,expected):
            self.assertAlmostEqual(actual,wanted,places=13)
        self.assertAlmostEqual(math.fsum(local)/len(local),float(rows[-1]['energy']),places=13)
        for row in rows:
            self.assertAlmostEqual(float(row['acceptance_rate']),float(row['energy_change_acceptance_rate'])+
                                   float(row['zero_energy_acceptance_rate']),places=14)
        self.assertIn('tile_side_sites=',(self.directory/'cpu.meta').read_text())

    def test_gpu_never_silently_runs_on_host(self):
        if REQUIRE_GPU:
            self.run_model('gpu','gpu')
        else:
            env=self.env
            self.env=dict(env,OMP_TARGET_OFFLOAD='DISABLED')
            p=self.run_model('nogpu','gpu',ok=False)
            self.assertTrue('GPU' in p.stderr or 'device' in p.stderr)

    def test_ensemble_selects_backend_and_threads(self):
        # Limit OpenMP environment for the child processes launched by the runner.
        config={'cases':[dict(q=3,nx=30,ny=29,spacing=.1,steps=2,cooling=.9,sweeps=2,
                              backend='cpu',threads=2,order='random')]}
        rows=run_ensemble(EXE,config,[8,9],self.directory/'ensemble')
        self.assertTrue(all(r['status']=='ok' for r in rows))
        for row in rows:
            text=(Path(row['directory'])/'run.meta').read_text()
            self.assertIn('backend=cpu',text)
            self.assertIn('openmp_max_threads=2',text)

    def test_benchmark_validates_trajectories_and_keeps_timings(self):
        rows=benchmark(EXE,self.directory/'benchmark',threads=[1,2],repeats=1,size=24,
                       spacing=.1,steps=2,sweeps=2)
        self.assertEqual(len(rows),3)
        self.assertTrue(all(row['matches_blocked_lattice'] for row in rows))
        self.assertTrue(all(row['annealing_wall_seconds']>=0 for row in rows))

    def test_real_gpu_matches_reference_and_restarts(self):
        if not REQUIRE_GPU:
            self.skipTest('Real GPU test requires --require-gpu and a GPU-offload build')
        self.run_model('reference','blocked')
        self.run_model('gpu','gpu')
        self.assertEqual((self.directory/'gpu.lat').read_bytes(),(self.directory/'reference.lat').read_bytes())
        self.run_model('whole','gpu',steps=5)
        self.run_model('split','gpu',steps=2)
        self.run_model('split','gpu',steps=3,seed=None)
        self.assertEqual((self.directory/'whole.lat').read_bytes(),(self.directory/'split-res.lat').read_bytes())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler',default=os.environ.get('FC','gfortran'))
    parser.add_argument('--flags')
    parser.add_argument('--compiler-kind',choices=['auto','gnu','nvhpc'],default='auto')
    parser.add_argument('--require-gpu',action='store_true',help='Fail unless a real offload device is exercised')
    args,remaining=parser.parse_known_args()
    config=compiler_config(args.compiler,args.compiler_kind)
    if args.flags is None:
        args.flags=('-O1 -g -Mbounds' if config['kind']=='nvhpc' else
                    '-O1 -g -std=f2008 -fcheck=all -ftrapv -ffpe-trap=invalid,zero,overflow')
    REQUIRE_GPU=args.require_gpu
    EXE=build(args.compiler,args.flags,ROOT/'build'/'parallel-tests'/
              ('colorplane.exe' if os.name=='nt' else 'colorplane'),openmp=True,compiler_kind=config['kind'])
    PROBE=EXE.parent/('parallel_probe.exe' if os.name=='nt' else 'parallel_probe')
    subprocess.run([args.compiler,*args.flags.split(),*config['openmp_flags'],'-I',str(EXE.parent),
                    config['module_flag'],str(EXE.parent),
                    str(EXE.parent/'colorplane_build_info.f90'),
                    *[str(ROOT/'src'/(name+'.f90')) for name in
                      ('colorplane_consts','colorplane_utils','colorplane_reliability','colorplane_rng','colorplane_parallel')],
                    str(ROOT/'tests'/'parallel_probe.f90'),'-o',str(PROBE)],cwd=EXE.parent,check=True)
    unittest.main(argv=[sys.argv[0],*remaining],verbosity=2)
