import unittest
from utils import metautils

#python -m unittest test.test_meta
class TypeLookup(unittest.TestCase):
    def test_find(self):
        illtype=metautils.gettype('SRR5009521')
        self.assertEqual(illtype, 'ILLUMINA')
    def test_nofind(self):
        try:
            metautils.gettype('smegma')
        except ValueError:
            pass
    def test_memory(self):
        try:
            metautils.getECS('SRR5009459','bytes')
        except Exception:
            pass
    def test_manifest(self):
        try:
            metautils.twoSampleComparison('Untreated HCT116','0.1 uM T3 treated HCT116','unvs01.txt')
def main():
    unittest.main()

if __name__ == "__main__":
    main()
