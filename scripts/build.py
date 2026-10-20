"""Build Colorplane with source provenance. Python standard library only."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def build(compiler="gfortran", flags="-O3", output=None):
    output = Path(output or ROOT / "build" / ("colorplane.exe" if os.name == "nt" else "colorplane")).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    sources = [ROOT / "src" / (name + ".f90") for name in
               ("colorplane_consts", "colorplane_utils", "colorplane_reliability", "colorplane")]
    hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        revision = "unknown"
    digest = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    identity = revision + ":" + digest
    generated = output.parent / "colorplane_build_info.f90"
    generated.write_text("module colorplane_build_info\n  implicit none\n"
                         "  character(len=*), parameter :: build_identity = &\n"
                         f"    '{identity}'\nend module\n", encoding="utf-8")
    command = [compiler, *shlex.split(flags), "-J", str(output.parent), "-I", str(output.parent),
               str(generated), *map(str, sources), "-o", str(output)]
    subprocess.run(command, cwd=output.parent, check=True)
    manifest = {"revision": revision, "source_sha256": hashes, "build_identity": identity,
                "command": command, "executable_sha256": hashlib.sha256(output.read_bytes()).hexdigest()}
    output.with_suffix(".build.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--flags", default="-O3")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    print(build(args.compiler, args.flags, args.output))
