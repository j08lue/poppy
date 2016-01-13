import unittest
import netCDF4
from poppy import grid

class TestLoad(unittest.TestCase):

    def test_find_k_depth(self):
        fname = './data/x3_0801-01.nc'
        with netCDF4.Dataset(fname) as ds:
            kz, z = grid.find_k_depth(ds, 0)
            self.assertEqual(kz, 0)
            self.assertEqual(kz, ds.variables['z_w'][kz]*1e-2)

if __name__ == '__main__':
    unittest.main()
