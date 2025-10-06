from typing import Literal

from lsst.afw.image import ExposureF
from lsst.afw.math import Warper
from lsst.afw.math import WarperConfig
from lsst.ip.diffim.subtractImages import AlardLuptonSubtractConfig
from lsst.ip.diffim.subtractImages import AlardLuptonSubtractTask
from lsst.source.injection import generate_injection_catalog
from lsst.source.injection import VisitInjectConfig
from lsst.source.injection import VisitInjectTask


def inject(
    exposure: ExposureF,
    ra: float,
    dec: float,
    mag: float,
):
    """Inject a synthetic stellar source into an exposure.

    Args:
        exposure: LSST exposure to inject source into
        ra: Right ascension in degrees
        dec: Declination in degrees
        mag: Source magnitude

    Returns:
        Tuple of (injected exposure, injection catalog)
    """
    EPSILON = 10**-7
    injection_catalog = generate_injection_catalog(
        ra_lim=[ra - EPSILON, ra + EPSILON],
        dec_lim=[dec - EPSILON, dec + EPSILON],
        number=1,
        source_type="Star",
        mag=[mag],
    )
    inject_config = VisitInjectConfig()
    inject_task = VisitInjectTask(config=inject_config)
    injected_output = inject_task.run(
        injection_catalogs=[injection_catalog],
        input_exposure=exposure.clone(),
        psf=exposure.getPsf(),
        photo_calib=exposure.getPhotoCalib(),
        wcs=exposure.getWcs(),
    )
    return injected_output.output_exposure, injected_output.output_catalog


def warp(
    science: ExposureF,
    template: ExposureF,
) -> ExposureF:
    """Warp template image to match science image WCS and bounding box.

    Args:
        science: Science exposure defining target WCS and bounds
        template: Template exposure to be warped

    Returns:
        Warped template exposure

    Note:
        From https://github.com/LSSTDESC/dia_improvement/blob/master/notebooks/dia_kernel_exploration.ipynb
    """
    warper_config = WarperConfig()
    warper = Warper.fromConfig(warper_config)
    warped_template = warper.warpExposure(science.getWcs(), template, destBBox=science.getBBox())
    # Add PSF.  I think doing this directly without warping is wrong.
    # At least the x,y mapping should be updated
    warped_template.setPsf(template.getPsf())
    return warped_template


def subtract(
    science,
    template,
    sources,
) -> ExposureF:
    """Perform Alard-Lupton difference imaging between science and template.

    Args:
        science: Science exposure
        template: Template exposure (will be warped to match science)
        sources: Source catalog for PSF matching

    Returns:
        Difference image exposure
    """
    template_warped = warp(science, template)
    config = AlardLuptonSubtractConfig()
    config.sourceSelector.value.unresolved.name = "base_ClassificationExtendedness_value"
    alTask = AlardLuptonSubtractTask(config=config)
    difference = alTask.run(template_warped, science, sources)
    return difference.difference
