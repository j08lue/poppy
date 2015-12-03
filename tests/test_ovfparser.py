import unittest
from poppy.ovfparser import parse_ovf_file

class TestLoad(unittest.TestCase):

    def test_parse_ovf_file(self):
        ovf_file = './data/gx1v6_overflow'
        overflows = parse_ovf_file(ovf_file)
        self.assertIsInstance(overflows, dict)

if __name__ == '__main__':
    unittest.main()
