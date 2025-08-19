from scipy.interpolate import RegularGridInterpolator
import numpy as np

import numpy as np
import warnings
from scipy.interpolate import RegularGridInterpolator

def _reconstruct_centers_from_model(model, axis="ra", grid_len=None):
    """
    Rebuild centers to match the (possibly padded) grid length using the model's
    ra_min/ra_max/dec_min/dec_max, bins, and padding.

    The model saved edges with length `bins`, which yields B = bins-1 actual bins
    before padding. After padding, the grid dimension is B + 2*padding.
    Centers are midpoints of each (possibly padded) bin.

    Parameters
    ----------
    model : dict
    axis : str
        "ra" or "dec"
    grid_len : int
        Target number of centers (should equal grid.shape[dim])

    Returns
    -------
    np.ndarray
    """
    assert axis in ("ra", "dec")
    amin = model[f"{axis}_min"]
    amax = model[f"{axis}_max"]
    bins = int(model.get("bins", 15))
    padding = int(model.get("padding", 1))

    # B is the number of *actual* bins before padding (because edges length = bins)
    B = bins - 1
    if B <= 0:
        raise ValueError(f"Invalid bins in model: bins={bins}")

    # native bin width (pre-padding)
    dA = (amax - amin) / B

    # original centers (length B)
    centers_core = amin + (np.arange(B) + 0.5) * dA

    # padded centers: extend by `padding` bins on both sides at same spacing
    if padding > 0:
        left = centers_core[0] - dA * np.arange(padding, 0, -1)
        right = centers_core[-1] + dA * np.arange(1, padding + 1)
        centers = np.concatenate([left, centers_core, right])
    else:
        centers = centers_core

    if grid_len is not None and len(centers) != grid_len:
        # If still mismatched, resample linearly across the padded span as a last resort.
        warnings.warn(
            f"{axis.upper()} centers length ({len(centers)}) != grid axis length ({grid_len}). "
            "Resampling centers to match grid shape."
        )
        # Rebuild centers uniformly across the total span of the padded grid:
        total_span = (amax - amin) + 2 * padding * dA
        a_min_padded = amin - padding * dA
        centers = a_min_padded + (np.arange(grid_len) + 0.5) * (total_span / grid_len)

    return centers

def zp_at_coord(model, ra, dec, margin=0.1):
    """
    Get zero-point correction for a given coordinate.

    Parameters
    ----------
    model : dict
        Zero-point calibration model loaded from JSON.
    ra : float
        Right Ascension in degrees.
    dec : float
        Declination in degrees.
    margin : float, optional
        Allowed margin (in degrees) outside the grid before warning.
        Default = 0.1 deg.

    Returns
    -------
    float
        Zero-point correction value (mag).
    """
    global_median = model.get("global_median", 0.0)

    
    if "grid" in model and "ra_centers" in model and "dec_centers" in model:
        ra_centers = np.array(model["ra_centers"])
        dec_centers = np.array(model["dec_centers"])
        grid = np.array(model["grid"])

        need_rebuild = (
            ra_centers.size == 0 or
            dec_centers.size == 0 or
            grid.shape[0] != ra_centers.size or
            grid.shape[1] != dec_centers.size
        )
        if need_rebuild:
            ra_centers = _reconstruct_centers_from_model(model, "ra", grid_len=grid.shape[0])
            dec_centers = _reconstruct_centers_from_model(model, "dec", grid_len=grid.shape[1])
            
        
        ra_min, ra_max = ra_centers.min(), ra_centers.max()
        dec_min, dec_max = dec_centers.min(), dec_centers.max()

        # Check bounds with margin
        if not (ra_min - margin <= ra <= ra_max + margin and
                dec_min - margin <= dec <= dec_max + margin):
            warnings.warn(
                f"Coordinate (RA={ra:.3f}, Dec={dec:.3f}) is outside "
                f"the grid range RA=[{ra_min:.3f}, {ra_max:.3f}], "
                f"Dec=[{dec_min:.3f}, {dec_max:.3f}]. "
                "Falling back to global median."
            )
            raise Exception(
                f"Coordinate (RA={ra}, Dec={dec}) is outside the grid range."
            )

        # Interpolator
        interpolator = RegularGridInterpolator(
            (ra_centers, dec_centers), grid,
            bounds_error=False,
            fill_value=np.nan
        )
        zp_value = interpolator([[ra, dec]])[0]

        # If interpolator returns NaN, fallback
        if np.isnan(zp_value):
            warnings.warn(
                f"Interpolation failed at (RA={ra:.3f}, Dec={dec:.3f}). "
                "Returning global median."
            )
            raise Exception(
                f"Interpolation failed at (RA={ra}, Dec={dec}). "
                "Falling back to global median."
            )

        return float(zp_value) + global_median

    # Fallback if no grid in model
    return global_median