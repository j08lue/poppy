import unittest
import glob
import numpy as np
from poppy import metrics

class TestLoad(unittest.TestCase):

    def test_get_amoc(self):
        """Test if get_amoc returns a proper time series"""
        ncfiles = sorted(glob.glob('./data/x3_0801-??.nc'))
        assert(len(ncfiles)>1)
        # get amoc timeseries
        df = metrics.get_amoc(ncfiles, latlim=(30,60), zlim=(500,9999), window_size=0)
        if metrics.use_pandas:
            self.assertIsInstance(df, metrics.pd.Series)
        self.assertEqual(len(df), len(ncfiles))

    def test_get_timeseries(self):
        """Test whether get_timeseries returns reasonable results"""
        ncfiles = sorted(glob.glob('./data/x3_0801-??.nc'))
        # get sea surface temperature timeseries
        df = metrics.get_timeseries(ncfiles, 
                varn='TEMP', grid='T', 
                reducefunc=np.nanmean, 
                latlim=(60,90), lonlim=None, k=0)
        # temperature below 100 degC
        self.assertTrue(np.mean(df) < 100)

if __name__ == '__main__':
    unittest.main()
