# Colorplane

Spin model simulator for plane colorings. Computes the ground state of a multicomponent spin system in a lattice, where the number of equal spin values unit distance apart is minimized. When interpreting spin values as colors, the output is an approximation to a plane coloring. The code thus attempts to produce approximate solutions to the (probabilistic) [Hadwiger-Nelson problem](https://en.wikipedia.org/wiki/Hadwiger%E2%80%93Nelson_problem).

The principles and methods used by the code are described in the associated paper [Lattice approach to plane colorings](https://arxiv.org/abs/1908.03880).

The code is interactive, and you can look at the input file examples in the "materials/inputs" subfolder for reasonable parameters. 

## Building and checking

Requires Python 3.9+ (standard library only) and a Fortran 2008 compiler; gfortran
is the default. From the repository root:

```sh
python scripts/build.py
python tests/test_reliability.py
```

For multicore execution:

```sh
python scripts/build.py --openmp
build/colorplane --backend cpu --threads 8 --seed 123 < materials/inputs/q2-10-1.inp
python tests/test_parallel.py
```

For NVIDIA/Linux GPU execution, independent-tile scheduling, correctness tests,
benchmarks, and cluster submission, see [Parallel spin updates](docs/parallelism.md).
The default `serial` mode preserves the original update schedule. `blocked`
provides a single-thread reference for `cpu` and `gpu` modes. GPU execution
requires an offload-enabled compiler/runtime and refuses silent CPU fallback.

The executable is `build/colorplane` (`build/colorplane.exe` on Windows).
Alternatively, `make -C src` retains the old `src/colorplane` location. Override
the compiler with `--compiler /path/to/gfortran` or `FC=...`, and compiler flags
with `--flags="-O3"`. The build records the Git revision, exact source hashes,
compiler options, and executable hash. Do not use unsafe floating-point flags
such as `-ffast-math` for scientific comparisons.

## Reproducible simulations

Existing fresh-run `.inp` files still work:

```sh
build/colorplane --seed 123 < materials/inputs/q2-10-1.inp
```

The default seed is **1**; choose different seeds explicitly for independent
runs. Repeatability requires the same executable/compiler, options, input,
and starting configuration. Compiler RNG implementations can differ, so seeds
are not a promise of bitwise agreement across compilers or versions.

Optional controls are `--order random` (a new permutation of sites in serial
mode, or of independent tile colors in parallel modes), `--kernel hat`,
`--reheat-every N`, and `--reheat-factor F` (default 2).
Defaults retain typewriter order, the cosine kernel, and geometric cooling.
Reheating multiplies the scheduled temperature before global stages N+1,
2N+1, etc.; it does not restore a discarded configuration. The hat kernel is
`max(0,1-|x|)`, an experimental sensitivity check, not a validated substitute
for all the smoothing properties discussed in the paper.

`--initial FILE` starts from a supplied lattice instead of random colors.
Use a new output root and specify its dimensions, spacing, and color count
in the input. Lattice files contain only colors and cannot verify their own
physical spacing. Starting configurations and the RNG state immediately
before the first Monte Carlo move are saved for every run.

### Measurements and output

| Extension | Contents |
| --- | --- |
| `.lat` | Final spin lattice |
| `.initial.lat` | Configuration before the first sweep of this run |
| `.rng-initial` | RNG state immediately before the first sweep |
| `.metrics.csv` | Measurements of the fixed configuration after each temperature stage |
| `.sum` | Human-readable parameters and exact stage energies |
| `.dat` | Serial sweep-time estimates; parallel modes use `.metrics.csv` instead |
| `.ngb` | Optional final local-energy map, normalized by the full interaction stencil |
| `.meta` | Seed, build/compiler identity, parameters, and measurement conventions |
| `.state` | Full-precision restart parameters, RNG state, and lattice consistency checksum |

The CSV records whole-lattice energy, interior energy and site count, exact-zero
site fraction, maximum local energy, and accepted moves per attempted move.
Accepted energy-changing and zero-energy moves are reported separately, with
the same denominator; their rates add to the total acceptance rate.
The interior excludes a conservative border wide enough to contain the full
nonzero interaction stencil. If no interior exists, its energy is `NaN` and
its count is zero. The zero fraction and maximum refer to the entire finite
lattice. Zero means exact floating-point zero after interaction weights at
or below `1e-14` have been discarded, not a proof about the continuum model.

Free boundaries still lower the whole-lattice energy. Interior measurements
help diagnose this bias but do not establish an infinite-volume limit.
These extra measurements require one additional full energy traversal per
temperature stage. No automatic convergence claim is made: sweep counts
remain user-specified.

### Restarting

As before, supplying a root with existing `.lat` and `.sum` files requests a
restart and writes a new root with `-res` appended. A matching `.state` is now
required. The code validates lattice parameters, annealing settings, compiler,
build identity, and a checksum of the saved lattice before restoring the RNG.
Do not pass `--seed` or `--initial` on an exact restart. Keep the original
cooling, sweep count, offset, order, kernel, and reheating parameters; the input
temperature-step count specifies **additional** stages. Offsets and reheating
use global stage indices, so a split run follows the same schedule as an
uninterrupted run. Restart validation deliberately rejects changed builds.

Legacy lattice files can be used through `--initial FILE` with a fresh seed
and output root. Their old rounded summaries cannot support an exact restart.
Existing destination files are never silently overwritten. Checkpoints are
written at successful completion, not after every stage; interruption recovery
within a run is not provided.

Spacing must be in `(0,0.5)` to exclude self-interaction. Every temperature
must have at least one sweep, and underflowing/overflowing temperature
schedules are rejected. Signed 16-bit storage supports all 2–255 colors.

## Repeated runs and comparison experiments

```sh
python scripts/ensemble.py --exe build/colorplane --config experiments.example.json --seeds 1 2 3 --output experiments/comparison
```

On Windows use the `.exe` suffix. The destination must be new. Each case/seed
gets its own input, command, console log, lattice, and diagnostics. The runner
writes a provenance manifest, `results.csv`, and `summary.json` with best,
mean, median, sample standard deviation, failure count, and zero-energy success
count. Failed simulations remain visible and cause a nonzero exit status.

The example is a small workflow demonstration, **not a converged scientific
study**. It varies physical size independently of spacing, and demonstrates
alternative update order, kernel, reheating, and rotated stripe initialization.
Case parameters are `q`, `nx`, `ny`, `spacing`, `steps`, `cooling`, and `sweeps`;
optional keys include `offset`, `order`, `kernel`, `reheat_every`,
`reheat_factor`, and `initial` (an explicit lattice path). Alternatively use
`pattern: "stripes"` or `"checkerboard"` for two colors, with optional `width`
in physical units and `angle` in degrees. Relative initial-file paths are
resolved from the directory in which the runner is invoked.

For scientific comparisons, increase both optimization budgets and the number
of independent runs. Compare interior observables across spacing and domain
size; compare algorithms at similar measured computation cost. Best energies
are candidate upper bounds for the finite lattice objective. A positive
residual energy is not evidence of a rigorous positive lower bound.

The test suite checks direct pair-enumerated energies, single-spin energy
differences (including the factor of two from ordered pairs), final maps,
seeded replay, exact restart continuity, 255-color storage, malformed input,
output preservation, stripe refinement, and the experiment runner. It uses
bounds checks and floating-point exception traps. It does not replace a full
scientific convergence study.

You should regard this software as a reference implementation: I would actually like you to write your own simulator which could use different choices of the simulated annealing algorithm. This would provide more information and could even give some important insight into the problem.

The source code in the folder "src" is subject to the GPL 3.0 License. See the license file in the main folder for the full text. The associated files in the folder "materials" and its subfolders are licensed under a Creative Commons Attribution-ShareAlike 4.0 International License. The material includes the input files used in the paper, as well as the lattices shown there. You can access [the full text](materials/License_CC) of the CC license.
