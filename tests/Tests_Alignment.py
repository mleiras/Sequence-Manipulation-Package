import unittest
from Alignment import Alignment

class TestAlignmentsLG(unittest.TestCase):
    
    def setUp(self):
        self.t1 = Alignment("ATGAAGGT", "AGAGAGGC", match_pairs = True, eqs = 2, diffs = 0, g = -1, alphabet = "ACGT")
        self.t2 = Alignment("LGPSGCASGIWTKSA", "TGPSGGSRIWTKSG", blosum62 = True, g = -8)
        self.t3 = Alignment("HGWAG", "PHSWG", blosum62 = True, g = -8, local = True)
        self.t4 = Alignment("GKYESVI", "KYVSSWI", blosum62 = True, g = -1)
        self.t5 = Alignment("GKYESVI", "KYVSSWI", blosum62 = True, local = True, g = -8)
        self.t6 = Alignment("ATGAAGGT", "AGAGAGGC", eqs = 2, diffs = -1, g = -1, alphabet = "ACGT")
        self.t7 = Alignment("HGWAG", "PHSWG", file = "blosum80.txt", g = -8, local = True)
    
    def test_alignScore(self):
        self.assertEqual(self.t1._alignScore(), ('global', 10))
        self.assertEqual(self.t2._alignScore(), ('global', 45))
        self.assertEqual(self.t3._alignScore(), ('local', 19))
        self.assertEqual(self.t4._alignScore(), ('global', 16))
        self.assertEqual(self.t5._alignScore(), ('local', 14))
    
    def test_getScore(self):
        self.assertEqual(self.t6._getScore('A', 'G'), -1)
        self.assertEqual(self.t6._getScore('A', 'A'), 2)
        self.assertEqual(self.t3._getScore('A', 'A'), 4)
        self.assertEqual(self.t3._getScore('A', 'Z'), -1)
        self.assertEqual(self.t7._getScore('A', 'A'), 5)
        self.assertEqual(self.t7._getScore('A', 'F'), -3)
        self.assertEqual(self.t7._getScore('*', '*'), 1)
    
    def test_alignment(self):
        self.assertEqual(self.t1.align_seq(), "Alignment:\nATGA-AGGT\nA-GAGAGGC")
        self.assertEqual(self.t1.align_seq(ties = True), "Alignment 1:\nATGA-AGGT\nA-GAGAGGC")
        self.assertEqual(self.t2.align_seq(), "Alignment:\nLGPSGCASGIWTKSA\nTGPSG-GSRIWTKSG")
        self.assertEqual(self.t2.align_seq(ties = True), "Alignment 1:\nLGPSGCASGIWTKSA\nTGPSG-GSRIWTKSG")
        self.assertEqual(self.t3.align_seq(), "Alignment:\nHGW\nHSW")
        self.assertEqual(self.t4.align_seq(ties = True), "Alignment 1:\nGKY-ESV-I\n-KYVSS-WI\n\nAlignment 2:\nGKY-ES-VI\n-KYVSSW-I")
        self.assertEqual(self.t5.align_seq(), "Alignment:\nKYES\nKYVS")
        
if __name__ == '__main__':
    # begin the unittest.main()
    unittest.main(argv=['first-arg-is-ignored'], exit = False)