"""Build Colorplane with source provenance. Python standard library only."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def compiler_config(compiler,kind="auto"):
    if kind=="auto":
        kind="nvhpc" if "nvfortran" in Path(compiler).name.lower() else "gnu"
    if kind not in ("gnu","nvhpc"):
        raise ValueError("Compiler kind must be gnu or nvhpc")
    return {"kind":kind,"module_flag":"-module" if kind=="nvhpc" else "-J",
            "openmp_flags":["-mp=gpu"] if kind=="nvhpc" else ["-fopenmp"]}


def build(compiler="gfortran", flags="-O3", output=None, openmp=False,compiler_kind="auto",spin_bits=None):
    config=compiler_config(compiler,compiler_kind)
    if spin_bits is None:
        spin_bits=32 if config["kind"]=="nvhpc" else 16
    if spin_bits not in (16,32):
        raise ValueError("Spin storage must be 16 or 32 bits")
    if config["kind"]=="nvhpc" and openmp and spin_bits==16:
        raise ValueError("NVIDIA OpenMP GPU builds need --spin-bits 32 to avoid the INT2 device runtime")
    output = Path(output or ROOT / "build" / ("colorplane.exe" if os.name == "nt" else "colorplane")).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    sources = [ROOT / "src" / (name + ".f90") for name in
               ("colorplane_consts", "colorplane_utils", "colorplane_reliability", "colorplane_rng",
                "colorplane_parallel", "colorplane")]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = "unknown"
    digest = hashlib.sha256(json.dumps({"sources":hashes,"spin_bits":spin_bits},sort_keys=True).encode()).hexdigest()
    identity = revision + ":" + digest
    generated = output.parent / "colorplane_build_info.f90"
    generated.write_text("module colorplane_build_info\n  implicit none\n"
                         f"  integer, parameter :: spin_kind = selected_int_kind({9 if spin_bits==32 else 4})\n"
                         "  character(len=*), parameter :: build_identity = &\n"
                         f"    '{identity}'\nend module\n", encoding="utf-8")
    command = [compiler, *shlex.split(flags), *(config["openmp_flags"] if openmp else []),
               config["module_flag"], str(output.parent), "-I", str(output.parent),
               str(generated), *map(str, sources), "-o", str(output)]
    subprocess.run(command, cwd=output.parent, check=True)
    manifest = {"revision": revision, "source_sha256": hashes, "build_identity": identity,
                "compiler_kind":config["kind"],
                "spin_bits":spin_bits,
                "command": command, "executable_sha256": hashlib.sha256(output.read_bytes()).hexdigest()}
    output.with_suffix(".build.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--flags", default="-O3")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--openmp", action="store_true", help="Enable CPU parallelism and compiled GPU targets")
    parser.add_argument("--compiler-kind",choices=["auto","gnu","nvhpc"],default="auto")
    parser.add_argument("--spin-bits",type=int,choices=[16,32],help="Default: GNU 16, NVIDIA 32")
    args = parser.parse_args()
    print(build(args.compiler, args.flags, args.output, args.openmp,args.compiler_kind,args.spin_bits))
