#!/usr/bin/env python
from __future__ import division, absolute_import, print_function

import unittest
from pint.observatory import Observatory, get_observatory
from pint.observatory.topo_obs import TopoObs
from pint import pulsar_mjd
import numpy as np
from astropy.time import Time
import os
from pint.config import datapath
from pinttestdata import testdir, datadir
from nose.tools import *


os.chdir(datadir)


class TestObservatory(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.test_obs = ['aro', 'ao', 'chime', 'drao']
        self.test_time = Time(np.linspace(55000, 58000, num=100),
                              scale='utc', format='pulsar_mjd')

    def test_get_obs(self):
        for tobs in self.test_obs:
            site = get_observatory(tobs, include_gps=False, include_bipm=True,
                                   bipm_version='BIPM2015')
            assert site, \
                "Observatory {} did not initilized correctly".format(tobs)

    def test_clock_corr(self):
        for tobs in self.test_obs:
            site = get_observatory(tobs, include_gps=True, include_bipm=True,
                                   bipm_version='BIPM2015')
            clock_corr = site.clock_corrections(self.test_time)
            assert len(clock_corr) == len(self.test_time)
            # Test one time
            clock_corr1 = site.clock_corrections(self.test_time[0])
            assert clock_corr1.shape == ()

    def test_get_TDBs(self):
        for tobs in self.test_obs:
            site = get_observatory(tobs, include_gps=True, include_bipm=True,
                                   bipm_version='BIPM2015')
            # Test default TDB calculation
            tdbs = site.get_TDBs(self.test_time)
            assert len(tdbs) == len(self.test_time)
            tdb1 = site.get_TDBs(self.test_time[0])
            assert tdb1.shape == (1,)

            # Test TDB calculation from ephemeris
            tdbs = site.get_TDBs(self.test_time, method='ephemeris',
                                 ephem='de430t')
            assert len(tdbs) == len(self.test_time)
            tdb1 = site.get_TDBs(self.test_time[0], method='ephemeris',
                                 ephem='de430t')
            assert tdb1.shape == (1,)


    def test_positions(self):
        for tobs in self.test_obs:
            site = get_observatory(tobs, include_gps=True, include_bipm=True,
                                   bipm_version='BIPM2015')
            posvel = site.posvel(self.test_time, ephem='de436')
            assert posvel.pos.shape == (3, len(self.test_time))
            assert posvel.vel.shape == (3, len(self.test_time))

    @raises(KeyError)
    def test_wrong_name(self):
        _ = get_observatory('Wrong_name')

    def test_wrong_path(self):
        # observatory clock correction path expections.
        fake_obs = TopoObs('Fake1', tempo_code='?', itoa_code='FK',
                            clock_fmt='tempo2', clock_file='fake2gps.clk',
                            clock_dir='TEMPO2', itrf_xyz=[0.00, 0.0, 0.0])
        site = get_observatory('Fake1', include_gps=True, include_bipm=True,
                               bipm_version='BIPM2015')
        with self.assertRaises(RuntimeError):
            _ = site.clock_corrections(self.test_time)

    def test_wrong_TDB_method(self):
        site = get_observatory('ao', include_gps=True, include_bipm=True,
                               bipm_version='BIPM2015')
        with self.assertRaises(ValueError):
            _ = site.get_TDBs(self.test_time, method='ephemeris')
        with self.assertRaises(ValueError):
            _ = site.get_TDBs(self.test_time, method='Unknown_method')
        with self.assertRaises(ValueError):
            _ = site.get_TDBs(self.test_time, method='ephemeris', ephem='de436')
