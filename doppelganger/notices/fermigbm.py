from datetime import datetime, timezone
import re

from doppelganger.types import CelestialCoords
from ._notices import get_param

__all__ = [
    "get_trigger_id",
    "get_issue_time",
    "get_obs_time",
    "get_localization_coords",
    "get_localization_error",
]

def _get_coords(notice: dict) -> dict:
    return notice["WhereWhen"]["ObsDataLocation"]["ObservationLocation"]["AstroCoords"]


_PATTERN_COORDS_SYSTEM_ID = r"(.*?)-(.*?)-(.*?)"


def _parse_coordsys(coords: dict) -> tuple[str, str, str]:
    """Helper returning timezone, frame and geographic reference from a notice."""
    if match := re.match(_PATTERN_COORDS_SYSTEM_ID, (_cs := coords["@coord_system_id"]), ):
        tz, frame, geo = match.groups()
    else:
        raise ValueError(f"Unknown coordinate system '{_cs}'.")
    return tz, frame, geo


def get_trigger_id(notice: dict) -> str | None:
    """Returns a string triggerId tag value from Fermi-GBM parsed notice.
    Returns None if no triggerId tag is found."""
    return get_param(notice["What"], "TrigID").get("@value", None)


def get_issue_time(notice: dict) -> datetime:
    """Returns a datetime object representing the time at which the notice itself was written
    to then be disseminated through GCN.."""
    dt = datetime.fromisoformat(notice["Who"]["Date"]["#text"])
    return dt


def get_obs_time(notice: dict) -> datetime:
    """Returns a datetime object representing the time of the observation."""
    coords = _get_coords(notice)
    tz_str, _, _ = _parse_coordsys(coords)

    if tz_str == "UTC":
        tz = timezone.utc
    else:
        raise ValueError(f"Unknown time zone {tz_str}")

    dt = datetime.fromisoformat(coords["Time"]["TimeInstant"]["ISOTime"]["#text"])
    return dt.replace(tzinfo=tz)


def get_localization_coords(notice: dict) -> CelestialCoords:
    """Returns a CelestialCoords object representing the most-likely sky position \
    for the celestial source.

    Usage:
        Convert to astropy's SkyCoord with:
        > from astropy.coordinates import SkyCoord
        > SkyCoord(**obsloc(notice)._asdict())
    """
    coords = _get_coords(notice)
    _, frame_str, _ = _parse_coordsys(coords)

    return CelestialCoords(
        ra=float(coords["Position2D"]["Value2"]["C1"]["#text"]),
        dec=float(coords["Position2D"]["Value2"]["C2"]["#text"]),
        unit=coords["Position2D"]["@unit"],
        frame=frame_str.lower(),
    )


def get_localization_error(notice: dict) -> float:
    """Returns the error of the celestial source."""
    coords = _get_coords(notice)
    return float(coords["Position2D"]["Error2Radius"]["#text"])