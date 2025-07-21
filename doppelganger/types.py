from typing import NamedTuple


class CelestialCoords(NamedTuple):
    ra: float
    dec: float
    unit: str
    frame: str

    def __str__(self) -> str:
        return f"RA={self.ra:.2f} {self.unit}, Dec={self.dec:.2f} {self.unit} ({self.frame})"
