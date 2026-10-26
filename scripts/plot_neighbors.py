"""Plot weighted equal-color neighbor probabilities. GPL-3.0-or-later."""
import argparse
from pathlib import Path

import numpy as np

from lattice_visualization import (crop_curve, image_axes, load_lattice,
                                   load_probability_map, neighbor_probabilities, pyplot, settings)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lattice", type=Path)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--map", type=Path, help="Read this .ngb map instead of ROOT.ngb")
    source.add_argument("--recompute", action="store_true", help="Ignore existing maps; calculate from lattice")
    parser.add_argument("--spacing", type=float, help="Override sibling .meta/.sum; required for recomputation")
    parser.add_argument("--kernel", choices=["cosine", "hat"], help="Default: sibling .meta, then cosine")
    parser.add_argument("--normalization", choices=["full", "available"], default="full",
                        help="full matches Colorplane; available conditions on neighbors inside domain")
    parser.add_argument("--output", type=Path, help="PNG/PDF/SVG; default: ROOT.neighbors.png")
    parser.add_argument("--data-output", type=Path, help="Also save the displayed probability matrix as text")
    parser.add_argument("--title")
    parser.add_argument("--vmax", type=float, default=1, help="Heatmap upper limit (default 1); lower limit is 0")
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    try:
        if not 0 < args.vmax <= 1 or args.dpi < 1:
            raise ValueError("vmax must be in (0,1] and DPI must be positive")
        lattice = load_lattice(args.lattice)
        spacing, kernel, _ = settings(args.lattice, args.spacing, args.kernel)
        map_path = args.map or args.lattice.with_suffix(".ngb")
        if args.data_output and args.data_output.resolve() in {args.lattice.resolve(), map_path.resolve()}:
            raise ValueError("Use a separate --data-output path to preserve the input lattice and neighbor map")
        if args.map and args.normalization != "full":
            raise ValueError("--map uses full normalization; use --recompute for available normalization")
        if not args.recompute and args.normalization == "full" and (args.map or map_path.exists()):
            probabilities = load_probability_map(map_path, lattice.shape)
            source_description = f"read {map_path} (full-stencil normalization)"
        else:
            if spacing is None:
                raise ValueError("No spacing found: supply --spacing or a sibling .meta/.sum file")
            probabilities = neighbor_probabilities(lattice, spacing, kernel, args.normalization)
            source_description = f"computed from lattice ({kernel}, {args.normalization} normalization)"
        plt = pyplot(args.show)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
        cmap = plt.get_cmap("magma").copy()
        cmap.set_bad("lightgray")
        image = image_axes(axes[0], probabilities, spacing, cmap=cmap, vmin=0, vmax=args.vmax)
        fig.colorbar(image, ax=axes[0], label="P(equal color)", extend="max" if args.vmax < 1 else "neither")
        axes[0].set_title("Equal-color neighbors near distance one")
        crop, means = crop_curve(probabilities)
        axes[1].plot(crop * (spacing or 1), means)
        axes[1].set_xlabel("Boundary width excluded (distance units)" if spacing else "Boundary layers excluded (sites)")
        axes[1].set_ylabel("Mean weighted probability")
        finite_means = means[np.isfinite(means)]
        axes[1].set_ylim(0, min(1, max(.01, float(finite_means.max()) * 1.12)) if finite_means.size else 1)
        axes[1].grid(alpha=.25)
        axes[1].set_title("Mean over the remaining central rectangle")
        fig.suptitle((args.title or args.lattice.name) + f" | {args.normalization} normalization")
        output = args.output or args.lattice.with_suffix(".neighbors.png")
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=args.dpi)
        if args.data_output:
            args.data_output.parent.mkdir(parents=True, exist_ok=True)
            np.savetxt(args.data_output, probabilities, fmt="%.17e",
                       header=f"Equal-color probability; {source_description}; spacing={spacing}")
        valid = probabilities[np.isfinite(probabilities)]
        print(source_description)
        print(f"Defined sites: {valid.size}/{probabilities.size}; mean: {valid.mean() if valid.size else float('nan'):.12g}")
        print(f"Saved {output}")
        if args.show:
            plt.show()
        plt.close(fig)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
