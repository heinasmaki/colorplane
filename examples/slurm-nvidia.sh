#!/usr/bin/env bash
# Adapt account, partition, GPU request syntax, and module setup to your cluster.
# Submit from the repository root, after compiling:
# sbatch examples/slurm-nvidia.sh build/nvidia/colorplane materials/inputs/q2-10-1.inp 123
#SBATCH --job-name=colorplane
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --mem=4G

set -euo pipefail
if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo 'Usage: sbatch slurm-nvidia.sh EXECUTABLE INPUT [SEED]' >&2
    exit 2
fi
cd "${SLURM_SUBMIT_DIR:?Submit this script with sbatch}"
executable=$(realpath "$1")
input=$(realpath "$2")
seed=${3:-123}
work="experiments/slurm-${SLURM_JOB_ID:?Missing job ID}"
mkdir -p experiments
mkdir "$work"
cp "$input" "$work/input.inp"
cd "$work"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_TARGET_OFFLOAD=MANDATORY
# Load compiler/runtime modules before sbatch, or add module commands here.
# Respect scheduler CUDA_VISIBLE_DEVICES; 0 is the first visible OpenMP device.
export OMP_DEFAULT_DEVICE=0
srun "$executable" --backend gpu --seed "$seed" < input.inp
