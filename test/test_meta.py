import unittest
from utils import metautils

#python -m unittest test.test_meta
class TypeLookup(unittest.TestCase):
    def test_find(self):
        illtype=metautils.getType('SRR5009521')
        self.assertEqual(illtype, 'ILLUMINA')
    def test_nofind(self):
        try:
            metautils.getType('smegma')
        except ValueError:
            pass
    def test_memory(self):
        #SRR5009459,2017-01-09 12:28:26,2016-11-09 18:37:32,113210576,16981586400,113210576,150,11506,,https://sra-download.ncbi.nlm.nih.gov/traces/sra8/SRR/004892/SRR5009459,SRX2340655,SA467_RNA-Seq_unstranded_HiSeq,RNA-Seq,cDNA,TRANSCRIPTOMIC,PAIRED,0,0,ILLUMINA,Illumina HiSeq 2000,SRP091981,PRJNA321560,,321560,SRS1757302,SAMN05919702,simple,9606,Homo sapiens,0.5 uM T3 treated HCT116,,,,,male,,no,,,,,,SRA492450,,public,FF65F819A7DD5AF3FA7B0978E56A25B5,91E47A9CDF6C97BE0843FE8B3AE3F368
        self.assertEqual(metautils.getECS('SRR5009459','bytes','rmats'),16777216000)
    def test_manifest(self):
        metautils.twoSampleComparisonManifest('Untreated HCT116','0.1 uM T3 treated HCT116','unvs01.txt')
def main():
    unittest.main()

if __name__ == "__main__":
    main()
