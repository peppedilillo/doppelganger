from datetime import datetime
import unittest
import zipfile as zf

from doppelganger.notices import fermigbm as gbm
from doppelganger.notices import parse_notice
from doppelganger.types import CelestialCoords

from . import TEST_DIR


class TestNoticesGBM(unittest.TestCase):
    def setUp(self):
        with zf.ZipFile(TEST_DIR / "data_notices_gbm.zip", "r") as fzip:
            self.notices = [
                (fname, fzip.read(fname).decode()) for fname in fzip.namelist()
            ]

    def test_parse_gbm_notice(self):
        for fname, content in self.notices:
            _ = parse_notice(content)

    def test_get_trigger_id(self):
        for fname, content in self.notices:
            notice = parse_notice(content)
            trig_id = gbm.get_trigger_id(notice)
            self.assertTrue(len(trig_id) > 0)
            self.assertTrue(isinstance(trig_id, str))

    def test_get_issue_time(self):
        for fname, content in self.notices:
            notice = parse_notice(content)
            dt = gbm.get_issue_time(notice)
            self.assertTrue(isinstance(dt, datetime))

    def test_get_obs_time(self):
        for fname, content in self.notices:
            notice = parse_notice(content)
            dt = gbm.get_obs_time(notice)
            self.assertTrue(isinstance(dt, datetime))

    def test_get_localization_coords(self):
        for fname, content in self.notices:
            notice = parse_notice(content)
            coords = gbm.get_localization_coords(notice)
            self.assertTrue(isinstance(coords, CelestialCoords))

    def test_get_localization_error(self):
        for fname, content in self.notices:
            notice = parse_notice(content)
            err = gbm.get_localization_error(notice)
            self.assertTrue(isinstance(err, float))
