"""Render a Colorplane .lat file with discrete colors. GPL-3.0-or-later."""
import argparse
from pathlib import Path

from lattice_visualization import lattice_image, load_lattice, pyplot, settings


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lattice", type=Path)
    parser.add_argument("--output", type=Path, help="PNG/PDF/SVG; default: ROOT.lattice.png")
    parser.add_argument("--spacing", type=float, help="Override sibling .meta/.sum; otherwise use site axes")
    parser.add_argument("--colors", type=int, help="Total color count, including absent colors")
    parser.add_argument("--title")
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--show", action="store_true", help="Also open an interactive plot window")
    args = parser.parse_args()
    try:
        lattice = load_lattice(args.lattice)
        spacing, _, colors = settings(args.lattice, args.spacing, colors=args.colors)
        if args.dpi < 1:
            raise ValueError("DPI must be positive")
        plt = pyplot(args.show)
        fig, ax = plt.subplots(figsize=(7, 6), layout="constrained")
        image, ticks = lattice_image(ax, lattice, spacing, colors)
        fig.colorbar(image, ax=ax, ticks=ticks, label="Color label")
        ax.set_title(args.title or args.lattice.name)
        output = args.output or args.lattice.with_suffix(".lattice.png")
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=args.dpi)
        print(f"Saved {output}")
        if args.show:
            plt.show()
        plt.close(fig)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
