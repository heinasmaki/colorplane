"""Scientific regression tests; no third-party Python dependencies required."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from build import build
from ensemble import initial_pattern, run_ensemble

EXE = None
PROBE = None


def kernel(distance, spacing, name):
    r = (distance - 1) / spacing
    value = (max(0, 1-abs(r)) if name == "hat" else
             (0.25*(1+math.cos(math.pi*r/2)) if abs(r)<2 else 0)) / spacing
    return value if value > 1e-14 else 0.0


def brute_force(lattice, spacing, name):
    """Enumerate unordered pairs directly, without production stencil arrays."""
    ny, nx = len(lattice), len(lattice[0])
    d = math.ceil(1 / spacing + 2)
    norm = math.fsum(kernel(math.hypot(x,y)*spacing, spacing, name)
                     for x in range(-d,d+1) for y in range(-d,d+1))
    values = [[0.0]*nx for _ in range(ny)]
    points = [(y,x) for y in range(ny) for x in range(nx)]
    for index, (y,x) in enumerate(points):
        for v,u in points[index+1:]:
            if lattice[y][x] == lattice[v][u]:
                weight = kernel(math.hypot(x-u,y-v)*spacing, spacing, name) / norm
                values[y][x] += weight
                values[v][u] += weight
    return values


class ReliabilityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="colorplane-test-")
        self.directory = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def run_model(self, root="run", steps=3, q=3, nx=9, ny=7, spacing=0.2, sweeps=4,
                  offset=0, args=(), ok=True):
        text = "\n".join(map(str, [q,ny,spacing,steps,0.9,sweeps,root,"y","y",nx,offset])) + "\n"
        proc = subprocess.run([str(EXE), *args], input=text, text=True, cwd=self.directory,
                              capture_output=True, timeout=30)
        if ok:
            self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
        else:
            self.assertNotEqual(proc.returncode,0,proc.stdout)
        return proc

    def rows(self, root):
        with (self.directory / f"{root}.metrics.csv").open() as source:
            return list(csv.DictReader(source))

    def lattice(self, root):
        return [[int(x) for x in row.split()] for row in
                (self.directory / f"{root}.lat").read_text().splitlines()]

    def test_seeded_replay_and_final_measurements(self):
        for root in ("one","two"):
            self.run_model(root,args=["--seed","123"])
        self.assertEqual(self.lattice("one"),self.lattice("two"))
        self.assertEqual(self.rows("one"),self.rows("two"))
        expected=brute_force(self.lattice("one"),0.2,"cosine")
        actual=[[float(x) for x in row.split()] for row in
                (self.directory / "one.ngb").read_text().splitlines()]
        for a,b in zip(sum(expected,[]),sum(actual,[])):
            self.assertAlmostEqual(a,b,places=13)
        last=self.rows("one")[-1]
        self.assertAlmostEqual(float(last["energy"]),math.fsum(sum(expected,[]))/63,places=13)
        self.assertEqual(int(last["interior_sites"]),0)
        self.assertTrue(math.isnan(float(last["interior_energy"])))
        for row in self.rows("one"):
            self.assertAlmostEqual(float(row["acceptance_rate"]),
                                   float(row["energy_change_acceptance_rate"])+
                                   float(row["zero_energy_acceptance_rate"]),places=14)

    def test_exact_restart_including_random_order_reheating_and_offset(self):
        for order in ("typewriter","random"):
            args=["--order",order,"--reheat-every","2","--reheat-factor","3"]
            self.run_model(order+"full",steps=5,sweeps=8,offset=1,args=args+["--seed","77"])
            self.run_model(order+"split",steps=2,sweeps=8,offset=1,args=args+["--seed","77"])
            self.run_model(order+"split",steps=3,sweeps=8,offset=1,args=args)
            self.assertEqual(self.lattice(order+"full"),self.lattice(order+"split-res"))
            self.assertEqual(self.rows(order+"full")[2:],self.rows(order+"split-res"))

    def test_different_seed_changes_initial_configuration(self):
        self.run_model("one",args=["--seed","1"])
        self.run_model("two",args=["--seed","2"])
        self.assertNotEqual((self.directory/"one.initial.lat").read_bytes(),
                            (self.directory/"two.initial.lat").read_bytes())

    def test_255_colors(self):
        self.run_model(q=255,args=["--seed","90"])
        values=sum(self.lattice("run"),[])
        self.assertTrue(all(1<=v<=255 for v in values))
        self.assertTrue(any(v>127 for v in values))

    def test_invalid_schedule_spacing_and_missing_input(self):
        self.run_model(steps=5,sweeps=4,offset=1,ok=False)
        self.run_model(spacing=0.5,ok=False)
        self.run_model(spacing="nan",ok=False)
        proc=subprocess.run([str(EXE)],input="",text=True,cwd=self.directory,capture_output=True,timeout=5)
        self.assertNotEqual(proc.returncode,0)

    def test_restart_parameter_and_checksum_validation(self):
        self.run_model()
        self.run_model(nx=10,ok=False)
        self.run_model(args=["--seed","3"],ok=False)
        self.run_model(args=["--kernel","hat"],ok=False)
        path=self.directory/"run.lat"
        text=path.read_text()
        path.write_text(text.replace("1","2",1) if "1" in text else text.replace("2","1",1))
        proc=self.run_model(ok=False)
        self.assertIn("checksum",proc.stderr.lower())

    def test_existing_outputs_are_preserved(self):
        self.run_model()
        self.run_model()  # creates run-res
        before=(self.directory/"run-res.lat").read_bytes()
        self.run_model(ok=False)
        self.assertEqual(before,(self.directory/"run-res.lat").read_bytes())

    def probe(self,lattice,spacing,name):
        path=self.directory/"probe.lat"
        path.write_text("\n".join(" ".join(map(str,row)) for row in lattice)+"\n")
        text=f"3 {len(lattice[0])} {len(lattice)} {spacing}\n{name}\n{path}\n"
        proc=subprocess.run([str(PROBE)],input=text,text=True,capture_output=True,timeout=30)
        self.assertEqual(proc.returncode,0,proc.stdout+proc.stderr)
        rows=proc.stdout.splitlines()
        return list(map(float,rows[0].split())), [list(map(float,row.split())) for row in rows[1:]]

    def test_independent_pair_energy_and_single_spin_delta(self):
        rng=random.Random(54)
        for name in ("cosine","hat"):
            for spacing in (0.2,0.125):
                lattice=[[rng.randint(1,3) for _ in range(9)] for _ in range(7)]
                summary,actual=self.probe(lattice,spacing,name)
                expected=brute_force(lattice,spacing,name)
                for a,b in zip(sum(actual,[]),sum(expected,[])):
                    self.assertAlmostEqual(a,b,places=13)
                self.assertAlmostEqual(summary[5],summary[6],places=13)

    def test_stripe_refinement_and_interior(self):
        errors=[]
        for width_sites in (8,16,32):
            a=math.sqrt(3)/2/width_sites
            n=8*width_sites
            lattice=[[1+(y//width_sites)%2 for _ in range(n)] for y in range(n)]
            summary,local=self.probe(lattice,a,"cosine")
            # A generic finite interior window contains partial stripe periods;
            # compare the infinite periodic reference using complete periods.
            margin=math.ceil(1/a+2)
            count=((n-2*margin)//(2*width_sites))*(2*width_sites)
            self.assertGreater(count,0)
            periodic_mean=math.fsum(local[y][n//2] for y in range(margin,margin+count))/count
            errors.append(abs(periodic_mean-1/3))
            self.assertLess(summary[0],summary[1])  # boundary bias is observable
            self.assertGreater(summary[4],0)
        self.assertLess(errors[-1],errors[0])
        self.assertLess(errors[-1],0.003)

    def test_uniform_interior_normalization_and_malformed_initial_lattices(self):
        for name in ("cosine","hat"):
            summary,local=self.probe([[1]*24 for _ in range(22)],0.2,name)
            self.assertAlmostEqual(summary[1],1.0,places=14)
            self.assertAlmostEqual(summary[3],1.0,places=14)
        for index,text in enumerate(("1 1\n", "0 1\n1 1\n", "1 1 1\n1 1\n",
                                     "1 1\n1 1\n1 1\n", "1 x\n1 1\n")):
            path=self.directory/f"bad-{index}.lat"
            path.write_text(text)
            self.run_model(root=f"bad-{index}",q=2,nx=2,ny=2,args=["--initial",str(path)],ok=False)

    def test_ensemble_and_reference_patterns(self):
        cases=[dict(q=2,nx=14,ny=12,spacing=0.2,steps=2,cooling=0.9,sweeps=2,
                    pattern="stripes",angle=30),
               dict(q=3,nx=8,ny=7,spacing=0.2,steps=2,cooling=0.9,sweeps=2,
                    order="random",kernel="hat")]
        output=self.directory/"ensemble"
        results=run_ensemble(EXE,{"cases":cases},[4,5],output)
        self.assertEqual(len(results),4)
        self.assertTrue(all(row["status"]=="ok" for row in results))
        summary=json.loads((output/"summary.json").read_text())
        self.assertTrue(all(case["successful_runs"]==2 for case in summary["cases"]))
        with self.assertRaises(FileExistsError):
            run_ensemble(EXE,{"cases":cases},[4,5],output)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--compiler",default=os.environ.get("FC","gfortran"))
    args,remaining=parser.parse_known_args()
    flags="-O1 -g -std=f2008 -fcheck=all -ffpe-trap=invalid,zero,overflow"
    EXE=build(args.compiler,flags,ROOT/"build"/"tests"/("colorplane.exe" if os.name=="nt" else "colorplane"))
    PROBE=EXE.parent/("energy_probe.exe" if os.name=="nt" else "energy_probe")
    subprocess.run([args.compiler,*flags.split(),"-I",str(EXE.parent),"-J",str(EXE.parent),
                    str(EXE.parent/"colorplane_build_info.f90"),
                    *[str(ROOT/"src"/(n+".f90")) for n in
                      ("colorplane_consts","colorplane_utils","colorplane_reliability")],
                    str(ROOT/"tests"/"energy_probe.f90"),"-o",str(PROBE)],cwd=EXE.parent,check=True)
    unittest.main(argv=[sys.argv[0],*remaining],verbosity=2)
