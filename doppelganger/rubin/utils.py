import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import astropy.visualization as vis


def remove_figure(fig):
    """
    Remove a figure to reduce memory footprint.
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
    s = ""
    description = service.search(
        "SELECT description "
        "FROM tap_schema.tables "
        "WHERE tap_schema.tables.schema_name = 'dp02_dc2_catalogs'"
        f"AND table_name = 'dp02_dc2_catalogs.{table_name}'"
    ).to_table().to_pandas().iloc[0]["description"]
    s += f"{table_name}: {description}\n\n"
    for i, row in service.search(
            "SELECT column_name, datatype, description, unit "
            "FROM tap_schema.columns "
            f"WHERE table_name = 'dp02_dc2_catalogs.{table_name}'"
    ).to_table().to_pandas().iterrows():
        s += f"{row.column_name}:  {row.description} ({row.datatype})\n"
    return s


def get_color_limits(img_data: np.array, scale: float | None = None) -> dict:
    zscale = vis.ZScaleInterval()
    if scale is None:
        return {"vmin": (_l := zscale.get_limits(img_data))[0], "vmax": _l[1]}
    elif 0. < scale < 100.:
        return {"vmin": np.percentile(img_data, scale), "vmax": np.percentile(img_data, 100 - scale)}
    else:
        raise ValueError("Paramter `scale` is a percentile and should be comprised between 0 and 100.")


def plot_host(img, ra, dec, scale: float | None = None, title: str | None = None, figsize: tuple[int, int] = (10, 10)):
    img_wcs = WCS(img.getWcs().getFitsMetadata())
    img_data = img.getImage().array
    coord_galhost = SkyCoord(
        ra=ra * u.degree,
        dec=dec * u.degree,
    )
    bbox_extent = (img.getBBox().beginX, img.getBBox().endX, img.getBBox().beginY, img.getBBox().endY)

    fig, ax = plt.subplots(1, figsize=figsize)
    plt.subplot(projection=img_wcs)
    plt.plot(*img_wcs.world_to_pixel(coord_galhost), 'r*', markerfacecolor="None", ms=20, label="Host galaxy")
    plt.imshow(img_data, cmap='gray', extent=bbox_extent, **get_color_limits(img_data, scale))
    if title is not None:
        plt.title(title)
    plt.axis('on')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid(color='grey', ls='solid')
    plt.legend(loc="lower left")
    return fig, ax


def plot_zoom(
        img,
        ra,
        dec,
        side_px=250,
        scale: float | None = None,
        marker: str = "+",  # matplotlib markers
        title: str | None = None,
        figsize: tuple[int, int] = (10, 10)
):
    wcs = WCS(img.getWcs().getFitsMetadata())
    coord_center = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        frame='icrs'
    )
    img_data = img.getImage().array
    bbox_extent = (img.getBBox().beginX, img.getBBox().endX, img.getBBox().beginY, img.getBBox().endY)
    center_x_pixel, center_y_pixel = wcs.world_to_pixel(coord_center)
    half_side_px = side_px // 2

    fig, ax = plt.subplots(1, figsize=figsize)
    plt.imshow(img_data, cmap='gray', extent=bbox_extent, origin="lower", **get_color_limits(img_data, scale))
    ax.plot(*(center_x_pixel, center_y_pixel), f'k{marker}', markerfacecolor="None", ms=40)
    if title is not None:
        ax.set_title(title)
    ax.set_xlim(center_x_pixel - half_side_px, center_x_pixel + half_side_px)
    ax.set_ylim(center_y_pixel - half_side_px, center_y_pixel + half_side_px)
    return fig, ax
