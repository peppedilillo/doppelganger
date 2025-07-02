from typing import Literal
from warnings import warn

import pandas as pd

from lsst.geom import SpherePoint, Point2I
import lsst.geom as geom
from lsst.daf.butler import Butler
from lsst.afw.math import Warper, WarperConfig
from lsst.source.injection import generate_injection_catalog
from lsst.source.injection import VisitInjectConfig, VisitInjectTask
from lsst.ip.diffim.subtractImages import AlardLuptonSubtractTask, AlardLuptonSubtractConfig


def calexp_contains(calexp, ra, dec):
    lsst_wcs = calexp.getWcs()
    sky_point = SpherePoint(ra * geom.degrees, dec * geom.degrees)
    pixel_point = Point2I(lsst_wcs.skyToPixel(sky_point))
    bbox = calexp.getBBox()
    if bbox.contains(pixel_point):
        return True
    return False


def seek(
    service,
    ra: float,
    dec: float,
    dataproduct_subtype: str | None = None,
    band: Literal["u", "g", "r", "i", "z", "y"] | None = None
):
    query = f"SELECT * FROM dp02_dc2_catalogs.ObsCore WHERE CONTAINS(POINT('ICRS', {ra:.4f}, {dec:.4f}), s_region) = 1"
    df = service.search(query).to_table().to_pandas()
    if dataproduct_subtype is not None:
        df = df[df["dataproduct_subtype"] == dataproduct_subtype]
    if band is not None:
        df = df[df["lsst_band"] == band]
    return df


def fetch(service, ra: float, dec: float, band: Literal["u", "g", "r", "i", "z", "y"], max_attempts: int=5):
    def dataId(visit: pd.DataFrame) -> dict:
        """
        Takes a row from the dataframe returned by seek and parses a dataId dictionary for the butler.
        """
        return {
            "instrument": "LSSTCam-imSim",
            "detector": visit["lsst_detector"].iloc[0],
            "visit": visit["lsst_visit"].iloc[0],
            "band": band,
            "physical_filter": f"{band}_sim_1.4",
        }

    butler = Butler('dp02', collections='2.2i/runs/DP0.2')
    visits_df = seek(service, ra, dec, dataproduct_subtype="lsst.calexp", band=band)
    attempts = 0
    while attempts < max_attempts and len(visits_df) > 0:
        visit = visits_df.sample(1)
        calexp = butler.get("calexp", dataId=(_dId := dataId(visit)))
        if not calexp_contains(calexp, ra, dec):
            attempts += 1
            warn(
                "Skipping an image since it does not actually contain the candidate source, likely due to s_region padding.")
            continue
        template = butler.get("goodSeeingDiff_templateExp", dataId=_dId)
        sources = butler.get('src', dataId=_dId)
        return calexp, template, sources
    if len(visits_df) == 0:
        raise RuntimeError(f"Could not find an image containing the target, visit table is empty.")
    raise RuntimeError(f"No images containing the transient found in {max_attempts} attempts.")


def inject(calexp, ra: float, dec: float, mag: float):
    EPSILON = 10 ** -7
    injection_catalog = generate_injection_catalog(
        ra_lim=[ra, ra + EPSILON],
        dec_lim=[dec, dec + EPSILON],
        number=1,
        source_type="Star",
        mag=[mag],
    )
    inject_config = VisitInjectConfig()
    inject_task = VisitInjectTask(config=inject_config)
    injected_output = inject_task.run(
        injection_catalogs=[injection_catalog],
        input_exposure=calexp.clone(),
        psf=calexp.getPsf(),
        photo_calib=calexp.getPhotoCalib(),
        wcs=calexp.getWcs(),
    )
    return injected_output.output_exposure, injected_output.output_catalog


def warp(science, template):
    """
    Warp input template image to WCS and Bounding Box of the science image.
    From: https://github.com/LSSTDESC/dia_improvement/blob/master/notebooks/dia_kernel_exploration.ipynb
    """
    warper_config = WarperConfig()
    warper = Warper.fromConfig(warper_config)
    warped_template = warper.warpExposure(science.getWcs(), template, destBBox=science.getBBox())
    # Add PSF.  I think doing this directly without warping is wrong.
    # At least the x,y mapping should be updated
    warped_template.setPsf(template.getPsf())
    return warped_template


def subtract(science, template, sources):
    template_warped = warp(science, template)
    config = AlardLuptonSubtractConfig()
    config.sourceSelector.value.unresolved.name = 'base_ClassificationExtendedness_value'
    alTask = AlardLuptonSubtractTask(config=config)
    difference = alTask.run(template_warped, science, sources)
    return difference.difference


def sfis_pipeline(service, ra, dec, mag, band):
    print("Retrieving visit table, choosing one visit at random.")
    calexp, template, sources = fetch(service, ra, dec, band)
    print("Starting source injection.")
    calexp_injected, calexp_catalog = inject(calexp, ra, dec, mag)
    print("Starting DIA.")
    calexp_difference = subtract(calexp_injected, template, sources)
    return {
        "science": calexp_injected,
        "template": template,
        "difference": calexp_difference,
    }