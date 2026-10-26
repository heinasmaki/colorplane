"""Shared readers and model-matched neighbor probabilities. GPL-3.0-or-later."""
from io import StringIO
import math
from pathlib import Path
import re

import numpy as np


def load_matrix(path):
    """Read rectangular Fortran/Python text arrays, including D exponents."""
    text = Path(path).read_text(encoding="utf-8").replace("D", "E").replace("d", "e")
    data = np.loadtxt(StringIO(text), ndmin=2)
    if not data.size or not np.isfinite(data).all():
        raise ValueError(f"{path}: expected a nonempty, finite numeric matrix")
    return data


def load_lattice(path):
    data = load_matrix(path)
    if np.any(data != np.floor(data)) or np.any((data < 1) | (data > 255)):
        raise ValueError("Lattice colors must be integers in 1..255")
    return data.astype(np.int16)


def settings(path, spacing=None, kernel=None, colors=None):
    """Prefer explicit options, then sibling .meta, then legacy .sum spacing."""
    path = Path(path)
    metadata = {}
    if path.with_suffix(".meta").exists():
        for line in path.with_suffix(".meta").read_text(encoding="utf-8").splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                metadata[key.strip()] = value.strip()
    if spacing is None and "spacing" in metadata:
        spacing = float(metadata["spacing"].replace("D", "E").replace("d", "e"))
    if spacing is None and path.with_suffix(".sum").exists():
        summary = path.with_suffix(".sum").read_text(encoding="utf-8", errors="replace")
        match = re.search(r"Lattice spacing:\s*([0-9.eEdD+-]+)", summary)
        if match:
            spacing = float(match[1].replace("D", "E").replace("d", "e"))
    if spacing is not None and (not math.isfinite(spacing) or not 0 < spacing < .5):
        raise ValueError("Spacing must be finite and in (0, 0.5), as in Colorplane")
    kernel = kernel or metadata.get("kernel", "cosine")
    if kernel not in ("cosine", "hat"):
        raise ValueError("Kernel must be cosine or hat")
    if colors is None and "q" in metadata:
        colors = int(metadata["q"])
    if colors is not None and not 2 <= colors <= 255:
        raise ValueError("Number of colors must be in 2..255")
    return spacing, kernel, colors


def interaction_kernel(spacing, name="cosine"):
    """Full retained stencil, with exactly the Fortran support and cutoff."""
    if not math.isfinite(spacing) or not 0 < spacing < .5:
        raise ValueError("Computing probabilities requires spacing in (0, 0.5)")
    radius = math.ceil(2 + 1 / spacing)
    y, x = np.ogrid[-radius:radius+1, -radius:radius+1]
    r = np.sqrt(x.astype(float)**2 + y.astype(float)**2) - 1 / spacing
    if name == "cosine":
        weights = np.where(np.abs(r) < 2, .25 * (1 + np.cos(np.pi * r / 2)), 0)
    elif name == "hat":
        weights = np.maximum(0, 1 - np.abs(r))
    else:
        raise ValueError("Kernel must be cosine or hat")
    weights /= spacing
    weights[weights <= 1e-14] = 0
    return weights


def fft_convolver(shape, weights):
    """Linear convolution with zero padding; never wrap opposite boundaries."""
    padded = tuple(1 << (n + k - 2).bit_length() for n, k in zip(shape, weights.shape))
    transform = np.fft.rfft2(weights, s=padded)
    dy, dx = (n // 2 for n in weights.shape)

    def convolve(values):
        result = np.fft.irfft2(np.fft.rfft2(values, s=padded) * transform, s=padded)
        return result[dy:dy+shape[0], dx:dx+shape[1]]

    return convolve


def neighbor_probabilities(lattice, spacing, kernel="cosine", normalization="full"):
    """Weighted P(same color) near distance one, using FFTs per present color."""
    if normalization not in ("full", "available"):
        raise ValueError("Normalization must be full or available")
    weights = interaction_kernel(spacing, kernel)
    convolve = fft_convolver(lattice.shape, weights)
    numerator = np.zeros(lattice.shape, dtype=float)
    for color in np.unique(lattice):
        mask = lattice == color
        neighbors = convolve(mask)
        numerator[mask] = neighbors[mask]
    if normalization == "full":
        result = numerator / weights.sum()
    else:
        # Missing exterior sites do not enter this conditional denominator.
        available = convolve(np.ones(lattice.shape))
        # Remove FFT noise where no retained offset reaches the lattice.
        minimum_weight = weights[weights > 0].min()
        result = np.full(lattice.shape, np.nan)
        np.divide(numerator, available, out=result, where=available > minimum_weight / 2)
    # FFT roundoff can create values a few ulps outside [0,1].
    return np.clip(result, 0, 1)


def load_probability_map(path, shape):
    data = load_matrix(path)
    if data.shape != shape:
        raise ValueError(f"Neighbor map shape {data.shape} does not match lattice {shape}")
    if np.any((data < -1e-13) | (data > 1 + 1e-13)):
        raise ValueError("Neighbor map values must be normalized probabilities in [0,1]")
    return np.clip(data, 0, 1)


def crop_curve(probabilities):
    """Mean over progressively cropped rectangles, ignoring undefined sites."""
    ny, nx = probabilities.shape
    y = np.minimum(np.arange(ny), np.arange(ny)[::-1])
    x = np.minimum(np.arange(nx), np.arange(nx)[::-1])
    layer = np.minimum(y[:, None], x[None, :])
    valid = np.isfinite(probabilities)
    sums = np.bincount(layer[valid], weights=probabilities[valid], minlength=layer.max()+1)
    counts = np.bincount(layer[valid], minlength=layer.max()+1)
    sums = np.cumsum(sums[::-1])[::-1]
    counts = np.cumsum(counts[::-1])[::-1]
    means = np.full(sums.shape, np.nan)
    np.divide(sums, counts, out=means, where=counts > 0)
    return np.arange(len(means)), means


def pyplot(show=False):
    import matplotlib
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def image_axes(ax, values, spacing, **kwargs):
    ny, nx = values.shape
    extent = (0, nx*spacing, 0, ny*spacing) if spacing else (0, nx, 0, ny)
    image = ax.imshow(values, origin="lower", interpolation="nearest", extent=extent,
                      aspect="equal", **kwargs)
    units = "distance units" if spacing else "sites"
    ax.set_xlabel(f"x ({units})")
    ax.set_ylabel(f"y ({units})")
    return image


def lattice_image(ax, lattice, spacing=None, colors=None):
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib import colormaps
    q = colors or max(2, int(lattice.max()))
    if lattice.max() > q:
        raise ValueError(f"Lattice contains color {lattice.max()}, exceeding --colors {q}")
    # Discrete colors are indexed by label, including absent colors.
    base = colormaps["tab10"] if q <= 10 else colormaps["tab20"] if q <= 20 else colormaps["hsv"]
    palette = base(np.arange(q) if q <= 20 else np.arange(q) / q)
    cmap = ListedColormap(palette)
    norm = BoundaryNorm(np.arange(.5, q+1.5), q)
    image = image_axes(ax, lattice, spacing, cmap=cmap, norm=norm)
    ticks = np.arange(1, q+1) if q <= 20 else np.unique(np.linspace(1, q, 10, dtype=int))
    return image, ticks
