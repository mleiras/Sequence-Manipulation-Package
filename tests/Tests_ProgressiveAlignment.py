import unittest
from ProgressiveAlignment import ProgressiveAlignment

class TestProgressive(unittest.TestCase):
    
    def setUp(self):
        self.t1 = ProgressiveAlignment(['TGACTA','TACGTA','TGGTA','GAT'])
        self.t2 = ProgressiveAlignment(['GTTGCACCA', 'GTCAGCA','TTCCCA','GCAGA'])
        self.t3 = ProgressiveAlignment(['ACTCAT','AGTCAT','ACGTCCT'])
        self.t4 = ProgressiveAlignment(['TGACTA','TACGTA','TGGTA','GAT'], g = -8, eqs = 2, diffs = -1)
        self.t5 = ProgressiveAlignment(['GTTGCACCA', 'GTCAGCA','TTCCCA','GCAGA'], g = 0, eqs = 2, diffs = -1)
        self.t6 = ProgressiveAlignment(['ACTCAT','AGTCAT','ACGTCCT'], g = 0, eqs = 4, diffs = -2)
        self.t7 = ProgressiveAlignment(['TGACTA','TACGTA','TGGTA','GAT'], profile = 'pssm', pseudo = 0.5)
        self.t8 = ProgressiveAlignment(['GTTGCACCA', 'GTCAGCA','TTCCCA','GCAGA'], profile = 'pssm', pseudo = 1)
        self.t9 = ProgressiveAlignment(['ACTCAT','AGTCAT','ACGTCCT'], profile = 'pssm', pseudo = 1)
        self.t10 = ProgressiveAlignment(['TGACTA','TACGTA','TGGTA','GAT'], g = -8, eqs = 2, diffs = -1, profile = 'pssm', pseudo = 0.5)
        self.t11 = ProgressiveAlignment(['GTTGCACCA', 'GTCAGCA','TTCCCA','GCAGA'], g = 0, eqs = 2, diffs = -1, profile = 'pssm', pseudo = 0.5)
        self.t12 = ProgressiveAlignment(['ACTCAT','AGTCAT','ACGTCCT'], g = 0, eqs = 4, diffs = -2, profile = 'pssm', pseudo = 0.5)
        self.t13 = ProgressiveAlignment(['ATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG','CTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAG','ACAAAAGCAACATCAAGGCTGCCTGGGGGAAGATTGGTGGCCATG'])
        self.t14 = ProgressiveAlignment(['LGPSGCASGIWTKSA', 'TGPSGGSRIWTKSG', 'GRIWTSGAVTKSA'], g = 0, eqs = 4, diffs = -2)
        self.t15 = ProgressiveAlignment(['SGIWTKSA', 'TGPSGGS', 'GRIWA'])
        
        
    def test_consensus(self):
        self.assertEqual(self.t1.consensus(self.t1.list_seqs[0],self.t1.list_seqs[1]), ('TAACTA'))        
        self.assertEqual(self.t2.consensus(self.t2.list_seqs[0],self.t2.list_seqs[1]), ('GTCACAA'))        
        self.assertEqual(self.t3.consensus(self.t3.list_seqs[0],self.t3.list_seqs[1]), ('ACTCAT'))        
        self.assertEqual(self.t4.consensus(self.t4.list_seqs[0],self.t4.list_seqs[1]), ('TAACTA'))        
        self.assertEqual(self.t5.consensus(self.t5.list_seqs[0],self.t5.list_seqs[1]), ('GTCACAA'))        
        self.assertEqual(self.t6.consensus(self.t6.list_seqs[0],self.t6.list_seqs[1]), ('ACTCAT'))
        self.assertEqual(self.t7.consensus(self.t7.list_seqs[0],self.t7.list_seqs[1]), ('TAACTA'))
        self.assertEqual(self.t8.consensus(self.t8.list_seqs[0],self.t8.list_seqs[1]), ('GTCACAA'))
        self.assertEqual(self.t9.consensus(self.t9.list_seqs[0],self.t9.list_seqs[1]), ('ACTCAT'))
        self.assertEqual(self.t10.consensus(self.t10.list_seqs[0],self.t10.list_seqs[1]), ('TAACTA'))
        self.assertEqual(self.t11.consensus(self.t11.list_seqs[0],self.t11.list_seqs[1]), ('GTCACAA'))
        self.assertEqual(self.t12.consensus(self.t12.list_seqs[0],self.t12.list_seqs[1]), ('ACTCAT'))
        self.assertEqual(self.t13.consensus(self.t13.list_seqs[0],self.t13.list_seqs[1]), ('ATGACTCCCTCCGACAAGAACAACGCCAATGCCAAATCCCTAAAG'))
        self.assertEqual(self.t14.consensus(self.t14.list_seqs[0],self.t14.list_seqs[1]), ('LGPSGCSSIITTSS'))
        self.assertEqual(self.t15.consensus(self.t15.list_seqs[0],self.t15.list_seqs[1]), ('SGISTKS'))
        
        
    def test_multiple_align(self):
        self.assertEqual(self.t1.progressive(), ['TGAC-TA', 'T-ACGTA', 'TG--GTA', '-GA--T-'])
        self.assertEqual(self.t2.progressive(), ['GTTGCACCA', 'G-T-CAGCA', '-TT-C-CCA', '---GCA-GA'])
        self.assertEqual(self.t3.progressive(), ['AC-TCAT', 'A-GTCAT', 'ACGTCCT'])
        self.assertEqual(self.t4.progressive(), ['TGACTA', 'TACGTA', 'TGG-TA', '-GA-T-'])
        self.assertEqual(self.t5.progressive(), ['GTTGCAC-CA', 'G-T-CA-GCA', '-TT-C-C-CA', '---GCA-G-A'])
        self.assertEqual(self.t6.progressive(), ['AC-TCAT', 'A-GTCAT', 'ACGTCCT'])
        self.assertEqual(self.t8.progressive(), ['GTTGCACCA', 'G-T-CAGCA', '-TT-C-CCA', '---GCA-GA'])
        self.assertEqual(self.t9.progressive(), ['AC-TCAT', 'A-GTCAT', 'ACGTCCT'])
        self.assertEqual(self.t10.progressive(), ['TGACTA', 'TACGTA', 'TGG-TA', '-GA-T-'])
        self.assertEqual(self.t11.progressive(), ['GTTGCAC-CA', 'G-T-CA-GCA', '-TT-C-C-CA', '---GCA-G-A'])
        self.assertEqual(self.t12.progressive(), ['AC-TCAT', 'A-GTCAT', 'ACGTCCT'])
        self.assertEqual(self.t13.progressive(), ['ATGAGT-CTCTCTGAT-AAGGACAAGG-CTGCTGTGAAAG-C-CCT-----AT-GG', '-CT--GTCTC-CTG--C-C-GACAAGACCAAC-GTCAAGGCCGCCTGG-GGT-AAG', 'ACAAAAGC-AAC--ATCAAGG-C-TG--C--CTG--GGGGAAGATTGGTGGCCATG'])
        self.assertEqual(self.t14.progressive(), ['L-GPSGCASG--IW-----TKSA', '-TGPSG---GSRIWTKSG------', '---------G-RIWT-SGAVTKSA'])
        self.assertEqual(self.t15.progressive(), ['SGIWTKSA', 'TGP-SGGS', 'GRIW---A'])
        
if __name__ == '__main__':
    # begin the unittest.main()
    unittest.main(argv=['first-arg-is-ignored'], exit = False)
    