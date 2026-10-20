"""Run independent seeded ensembles and summarize finite-lattice candidate energies."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import time


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def initial_pattern(case, path):
    """Reference patterns in physical coordinates; angle is in degrees."""
    pattern = case.get("pattern", "random")
    if pattern == "random":
        return None
    if pattern not in ("stripes", "checkerboard") or case["q"] != 2:
        raise ValueError("stripes/checkerboard patterns require q=2")
    width = case.get("width", math.sqrt(3) / 2 if pattern == "stripes" else 1.0)
    angle = math.radians(case.get("angle", 0.0))
    if width <= 0 or not math.isfinite(width) or not math.isfinite(angle):
        raise ValueError("Invalid reference pattern width/angle")
    ca, sa, a = math.cos(angle), math.sin(angle), case["spacing"]
    rows = []
    for y in range(case["ny"]):
        row = []
        for x in range(case["nx"]):
            u, v = a * (ca*x + sa*y), a * (-sa*x + ca*y)
            index = math.floor(v / width)
            if pattern == "checkerboard":
                index += math.floor(u / width)
            row.append(str(1 + index % 2))
        rows.append(" ".join(row))
    path.write_text("\n".join(rows) + "\n", encoding="ascii")
    return path


def run_ensemble(exe, config, seeds, output):
    exe, output = Path(exe).resolve(), Path(output).resolve()
    if not exe.is_file():
        raise ValueError("Executable does not exist")
    if len(set(seeds)) != len(seeds) or any(s < 0 or s > 2147483646 for s in seeds):
        raise ValueError("Seeds must be unique integers in 0..2147483646")
    cases = config["cases"]
    if not cases or not seeds:
        raise ValueError("At least one case and seed are required")
    # Never reuse an output directory: stale files would silently trigger restarts.
    output.mkdir(parents=True, exist_ok=False)
    manifest = {"config": config, "seeds": seeds, "executable": str(exe),
                "executable_sha256": sha256(exe)}
    build_manifest = exe.with_suffix(".build.json")
    if build_manifest.exists():
        manifest["build"] = json.loads(build_manifest.read_text(encoding="utf-8"))
        if manifest["build"]["executable_sha256"] != manifest["executable_sha256"]:
            raise ValueError("Executable does not match build manifest")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    results, summaries = [], []
    for index, case in enumerate(cases, 1):
        group = []
        for seed in seeds:
            directory = output / f"case-{index:03d}" / f"seed-{seed}"
            directory.mkdir(parents=True)
            nx, ny = case["nx"], case["ny"]
            offset = case.get("offset", 0)
            text = "\n".join(map(str, [case["q"], ny, case["spacing"], case["steps"],
                               case["cooling"], case["sweeps"], "run", "y", "y", nx, offset])) + "\n"
            (directory / "input.inp").write_text(text, encoding="ascii")
            command = [str(exe), "--seed", str(seed), "--order", case.get("order", "typewriter"),
                       "--kernel", case.get("kernel", "cosine"), "--reheat-every",
                       str(case.get("reheat_every", 0)), "--reheat-factor", str(case.get("reheat_factor", 2.0))]
            command += ["--backend", case.get("backend", "serial")]
            if "threads" in case:
                command += ["--threads", str(case["threads"])]
            initial = initial_pattern(case, directory / "pattern.lat")
            if "initial" in case:
                if initial:
                    raise ValueError("Choose either initial file or pattern")
                initial = directory / "pattern.lat"
                shutil.copyfile(Path(case["initial"]).resolve(), initial)
            if initial:
                command += ["--initial", str(initial)]
            (directory / "command.json").write_text(json.dumps(command), encoding="utf-8")
            started = time.perf_counter()
            with (directory / "console.log").open("w", encoding="utf-8") as log:
                proc = subprocess.run(command, input=text, text=True, cwd=directory, stdout=log, stderr=log)
            result = {"case": index, "seed": seed, "status": "ok" if proc.returncode == 0 else "failed",
                      "seconds": time.perf_counter() - started, "directory": str(directory)}
            if proc.returncode == 0:
                with (directory / "run.metrics.csv").open(encoding="ascii") as source:
                    final = list(csv.DictReader(source))[-1]
                result.update({key: float(value) for key, value in final.items()})
                group.append(result)
            results.append(result)
        energies = [r["energy"] for r in group]
        summaries.append({"case": index, "parameters": case, "successful_runs": len(group),
                          "failed_runs": len(seeds)-len(group), "best_energy": min(energies) if energies else None,
                          "median_energy": statistics.median(energies) if energies else None,
                          "mean_energy": statistics.mean(energies) if energies else None,
                          "sample_stdev": statistics.stdev(energies) if len(energies)>1 else None,
                          "zero_energy_runs": sum(e == 0 for e in energies)})
        # Save partial progress after each case; failures remain visible.
        fields = list(dict.fromkeys(key for row in results for key in row))
        with (output / "results.csv").open("w", newline="", encoding="utf-8") as target:
            writer = csv.DictWriter(target, fields)
            writer.writeheader()
            writer.writerows(results)
        (output / "summary.json").write_text(json.dumps({
            "interpretation": "Best energies are candidate upper bounds for these finite lattice models; "
                              "failure to find zero does not prove a positive lower bound.",
            "cases": summaries}, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--seeds", required=True, type=int, nargs="+")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = run_ensemble(args.exe, json.loads(args.config.read_text(encoding="utf-8")), args.seeds, args.output)
    print(f"Saved {len(rows)} runs to {args.output}")
    raise SystemExit(any(row["status"] != "ok" for row in rows))
