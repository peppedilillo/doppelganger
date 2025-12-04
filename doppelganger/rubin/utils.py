from typing import Sequence

from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.visualization as vis
from astropy.wcs import WCS
from lsst import geom as geom
from lsst.afw.image import ExposureF
from lsst.geom import Point2I
from lsst.geom import SpherePoint
import matplotlib.pyplot as plt
import numpy as np


def remove_figure(fig):
    """Remove a figure and free associated memory.

    Args:
        fig: Matplotlib Figure object or (Figure, Axes) tuple to remove

    Note:
        From Rubin DP0.2 Tutorial notebooks by the Community Science Team.
    """
    import gc

    plt.show()

    # so that we can pass either a fig or a (fig, ax) tuple.
    if isinstance(fig, tuple):
        fig = fig[0]

    # Get the axes and clear their images
    for ax in fig.get_axes():
        for im in ax.get_images():
            im.remove()
    fig.clf()  # Clear the figure
    plt.close(fig)  # Close the figure
    gc.collect()  # Call the garbage collector


def table_summary(service, table_name: str):
    """Generate a summary of a TAP table schema including column descriptions.

    Args:
        service: PyVO TAP service object
        table_name: Name of the table in dp02_dc2_catalogs schema

    Returns:
        Formatted string with table description and column details
    """
    s = ""
    description = (
        service.search(
            "SELECT description "
            "FROM tap_schema.tables "
            "WHERE tap_schema.tables.schema_name = 'dp02_dc2_catalogs'"
            f"AND table_name = 'dp02_dc2_catalogs.{table_name}'"
        )
        .to_table()
        .to_pandas()
        .iloc[0]["description"]
    )
    s += f"{table_name}: {description}\n\n"
    for i, row in (
        service.search(
            "SELECT column_name, datatype, description, unit "
            "FROM tap_schema.columns "
            f"WHERE table_name = 'dp02_dc2_catalogs.{table_name}'"
        )
        .to_table()
        .to_pandas()
        .iterrows()
    ):
        s += f"{row.column_name}:  {row.description} ({row.datatype})\n"
    return s


def str_ivoa2adql(ivoa_str: str, frame: str = "ICRS") -> str:
    """Convert IVOA STC-S region string to ADQL geometric function.

    Args:
        ivoa_str: IVOA STC-S region string (e.g., "CIRCLE 150.0 2.0 0.5")
        frame: Coordinate reference frame. Defaults to "ICRS"

    Returns:
        ADQL geometric function string (e.g., "CIRCLE('ICRS', 150.0, 2.0, 0.5)")
    """
    figure, *coords = ivoa_str.strip().split()
    coord_string = ", ".join(f"{c:.7f}" for c in map(float, coords))
    return f"{figure}('{frame}', {coord_string})"


def get_mask(image: ExposureF, mask_names: str | list[str]) -> np.ndarray:
    """Extract a boolean mask array from an LSST image for specified mask planes.

    Args:
        image: LSST image object with mask information
        mask_names: Single mask plane name or list of mask plane names

    Returns:
        Boolean array. True indicates pixels flagged in any of the specified mask planes
    """
    mask = image.getMask()
    if isinstance(mask_names, str):
        mask_names = [mask_names]

    mask_array = mask.getArray()
    out = np.zeros_like(mask_array)
    for mask_name in mask_names:
        target_bit = mask.getMaskPlane(mask_name)
        out |= (mask_array & (2**target_bit)) != 0
    return out


def get_color_limits(img_data: np.array, scale: float | None = None) -> dict:
    """Calculate color scale limits for image display.

    Args:
        img_data: Image array data
        scale: Percentile value for contrast clipping (0-100). If None, uses ZScale algorithm

    Returns:
        Dictionary with 'vmin' and 'vmax' keys for matplotlib color scaling
    """
    zscale = vis.ZScaleInterval()
    if scale is None:
        return {"vmin": (_l := zscale.get_limits(img_data))[0], "vmax": _l[1]}
    elif 0.0 < scale < 100.0:
        return {
            "vmin": np.percentile(img_data, scale),
            "vmax": np.percentile(img_data, 100 - scale),
        }
    else:
        raise ValueError("Paramter `scale` is a percentile and should be comprised between 0 and 100.")


def plot_with_coords(
    image: ExposureF,
    coords: Sequence[tuple[float, float]] = (),
    figsize: tuple[int, int] = (10, 10),
    scale: float | None = None,
):
    """Plot an astronomical image, optionally with markers for specific coordinates.

    Args:
        image: Astronomical image object with getWcs() and getImage() methods.
        coords: Sequence of (ra, dec) tuples in degrees for sky coordinates to mark on first image. Optional.
        figsize: Tuple of (width, height) in inches for figure size. Defaults to (16, 9).
        scale: Float value between 0 and 100 controlling contrast via percentile clipping. If None, uses ZScale. Defaults to None.

    Returns:
        Tuple of (matplotlib Figure, list of matplotlib Axes objects).
    """
    wcs = WCS(image.getWcs().getFitsMetadata())
    image_data = image.getImage().array

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ax.imshow(image_data, cmap="gray", **get_color_limits(image_data, scale=scale))
    for ra, dec in coords:
        coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs")
        ax.plot(*wcs.world_to_pixel(coord), "ro", markerfacecolor="None", ms=20)
    ax.set(xticks=[], yticks=[], xlabel="", ylabel="")
    plt.tight_layout()
    plt.show()
    return fig, ax


def plot_side_by_side(
    image1: ExposureF,
    image2: ExposureF,
    coords1: Sequence[tuple[float, float]] = (),
    coords2: Sequence[tuple[float, float]] = (),
    figsize: tuple[int, int] = (16, 9),
    scale: float | None = None,
):
    """Plot two astronomical images side by side, optionally with markers for specific coordinates.

    Args:
        image1: Astronomical image object with getWcs() and getImage() methods.
        image2: Astronomical image object with getWcs() and getImage() methods.
        coords1: Sequence of (ra, dec) tuples in degrees for sky coordinates to mark on first image. Optional.
        coords2: Sequence of (ra, dec) tuples in degrees for sky coordinates to mark on second image. Optional.
        figsize: Tuple of (width, height) in inches for figure size. Defaults to (16, 9).
        scale: Float value between 0 and 100 controlling contrast via percentile clipping. If None, uses ZScale. Defaults to None.

    Returns:
        Tuple of (matplotlib Figure, list of matplotlib Axes objects).
    """
    fig = plt.figure(figsize=figsize)
    axs = []
    for i, (image, coords) in enumerate(zip((image1, image2), (coords1, coords2))):
        wcs = WCS(image.getWcs().getFitsMetadata())
        image_data = image.getImage().array
        ax = fig.add_subplot(1, 2, i + 1, projection=wcs)
        ax.imshow(image_data, cmap="gray", **get_color_limits(image_data, scale))
        for ra, dec in coords:
            coord = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame="icrs")
            ax.plot(*wcs.world_to_pixel(coord), "ko", markerfacecolor="None", ms=20)
        ax.set(xticks=[], yticks=[], xlabel="", ylabel="")
        axs.append(ax)
    plt.tight_layout()
    plt.show()
    return fig, axs


def plot_zoom(
    img,
    ra,
    dec,
    side_px=250,
    figsize: tuple[int, int] = (10, 10),
    scale: float | None = None,
    marker: str = "+",  # matplotlib markers
    title: str | None = None,
):
    """Plot a zoomed-in view of an astronomical image centered on specific coordinates.

    Args:
        img: Astronomical image object with getWcs(), getImage(), and getBBox() methods.
        ra: Right ascension in degrees.
        dec: Declination in degrees.
        side_px: Side length of the zoomed region in pixels. Defaults to 250.
        figsize: Tuple of (width, height) in inches for figure size. Defaults to (10, 10).
        scale: Float value between 0 and 100 controlling contrast via percentile clipping. If None, uses ZScale. Defaults to None.
        marker: Matplotlib marker style for the central coordinate. Defaults to "+".
        title: Optional title for the plot. Defaults to None.

    Returns:
        Tuple of (matplotlib Figure, matplotlib Axes object).
    """
    wcs = WCS(img.getWcs().getFitsMetadata())
    coord_center = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    img_data = img.getImage().array
    bbox_extent = (
        img.getBBox().beginX,
        img.getBBox().endX,
        img.getBBox().beginY,
        img.getBBox().endY,
    )
    center_x_pixel, center_y_pixel = wcs.world_to_pixel(coord_center)
    half_side_px = side_px // 2

    fig, ax = plt.subplots(1, figsize=figsize)
    plt.imshow(
        img_data,
        cmap="gray",
        extent=bbox_extent,
        origin="lower",
        **get_color_limits(img_data, scale),
    )
    ax.plot(*(center_x_pixel, center_y_pixel), f"k{marker}", markerfacecolor="None", ms=40)
    if title is not None:
        ax.set_title(title)
    ax.set_xlim(center_x_pixel - half_side_px, center_x_pixel + half_side_px)
    ax.set_ylim(center_y_pixel - half_side_px, center_y_pixel + half_side_px)
    return fig, ax


def image_contains(image: ExposureF, ra: float, dec: float) -> bool:
    """Check whether sky coordinates fall within image boundaries.

    Args:
        image: LSST image object with WCS and bounding box
        ra: Right ascension in degrees
        dec: Declination in degrees

    Returns:
        True if coordinates fall within image bounds, False otherwise
    """
    lsst_wcs = image.getWcs()
    sky_point = SpherePoint(ra * geom.degrees, dec * geom.degrees)
    pixel_point = Point2I(lsst_wcs.skyToPixel(sky_point))
    return image.getBBox().contains(pixel_point)
