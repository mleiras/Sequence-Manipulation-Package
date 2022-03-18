import unittest
from Blast import Blast
from Blast import SimpleBlast

class TestBlast(unittest.TestCase):
    
    def setUp(self):
        self.t1 = Blast('AATATGTTATATAATAATATTT', 'AATATAT', 3)
        self.t2 = Blast('aataTGTTAtatAATAATATTT', 'AATAtaT', 3)
        self.t3 = Blast('GGGCCGATAAGTGATCAGGTAGCTA', 'GGT', 2)
        self.t4 = Blast('TGGCAACGACCATG', 'CCCA', 3)
        self.t5 = Blast('MNOKLYLPBMMRWYYHJFSLPPNYBB', 'MMRWYY', 4)
        
    def test_query_map(self):
        self.assertEqual(self.t1.query_map(), {'AAT': [0], 'ATA': [1, 3], 'TAT': [2, 4]})
        self.assertEqual(self.t2.query_map(), {'AAT': [0], 'ATA': [1, 3], 'TAT': [2, 4]})
        self.assertEqual(self.t3.query_map(), {'GG': [0], 'GT': [1]})
        self.assertEqual(self.t4.query_map(), {'CCC': [0], 'CCA': [1]})
        self.assertEqual(self.t5.query_map(), {'MMRW': [0], 'MRWY': [1], 'RWYY': [2]})
    
    def test_hits(self):
        self.assertEqual(self.t1.hits(), [(0, 0), (0, 12), (0, 15), (1, 1), (1, 8), (1, 10), (1, 13), (1, 16), (3, 1), (3, 8), (3, 10), (3, 13), (3, 16), (2, 2), (2, 7), (2, 9), (2, 17), (4, 2), (4, 7), (4, 9), (4, 17)])
        self.assertEqual(self.t2.hits(), [(0, 0), (0, 12), (0, 15), (1, 1), (1, 8), (1, 10), (1, 13), (1, 16), (3, 1), (3, 8), (3, 10), (3, 13), (3, 16), (2, 2), (2, 7), (2, 9), (2, 17), (4, 2), (4, 7), (4, 9), (4, 17)])
        self.assertEqual(self.t3.hits(), [(0, 0), (0, 1), (0, 17), (1, 10), (1, 18)])
        self.assertEqual(self.t4.hits(), [(1, 9)])
        self.assertEqual(self.t5.hits(), [(0, 9), (1, 10), (2, 11)])
    
    def test_extend_hit(self):
        self.assertEqual(self.t1.extend_hit((1, 16)), (0, 15, 7, 6))
        self.assertEqual(self.t1.extend_hit((0, 0)), (0, 0, 7, 6))
        self.assertEqual(self.t1.extend_hit((4, 2)), (2, 0, 5, 4))
        self.assertEqual(self.t3.extend_hit((1, 10)), (0, 9, 3, 2))
        self.assertEqual(self.t4.extend_hit((1, 9)), (0, 8, 4, 3))
        self.assertEqual(self.t5.extend_hit((1, 10)), (0, 9, 6, 6))
        
    def test_best_hit(self):
        self.assertEqual(self.t1.best_hit(), (0, 0, 7, 6))
        self.assertEqual(self.t2.best_hit(), (0, 0, 7, 6))
        self.assertEqual(self.t3.best_hit(), (0, 17, 3, 3))
        self.assertEqual(self.t4.best_hit(), (1, 9, 4, 3))
        self.assertEqual(self.t5.best_hit(), (0, 9, 6, 6))

        
class TestSimpleBlast(unittest.TestCase):
    
    def setUp(self):
        self.t1 = SimpleBlast(['AATATGTTATATAATAATATTT', 'GGGTTATATAATAATATTTA', 'AATATGTTATATTTCAGGATCTATTT', 'AATATATAATATATAATATATAATATATAATATAT'], 3)
        self.t2 = SimpleBlast(['AATATGTTATAtaatAATATTT', 'gggttatataataatattta', 'aatatGTTATATTTCAGGATCTattt', 'aatatATAATATATAATATATAATATATaatatat'], 3)
        self.t3 = SimpleBlast(['TGGGTAGATGATGACCCGACC', 'CCCGCCATGCACCGA', 'GAGAGATCGACAGCT', 'TTTTAAACAAAT'], 3)
        self.t4 = SimpleBlast(['ACTGATCACACGTACTAGCTTCGTATACGTcgtacttgaacaaaagtaaacagtacagtact', 'gggggtacatgacgatacgatcagt', 'agactacgatcagttttttcagtacccccgagccca', 'acgacgacgatttac', 'CCCGCCATGCACCGA', 'GAGAGATCGACAGCT', 'TTTTAAACAAAT'], 3)
        self.t5 = SimpleBlast(['ACTGATCACACGTACTAGCTTCGTATACGTcgtacttgaacaaaagtaaacagtacagtact', 'gggggtacatgacgatacgatcagt', 'agactacgatcagttttttcagtacccccgagccca', 'acgacgacgatttac', 'CCCGCCATGCACCGA', 'GAGAGATCGACAGCT', 'TTTTAAACAAAT'], 6)
        self.t6 = SimpleBlast(["HSGVYLDDRDKPHFQDLNSFILDMSNNQRD", "RFIIQTPYVRWKQRPEFWYPYEFEAVQTNS", "MTNQNHEVCHFYMCIFMLFIHQKKSSPYLS", "GYLWHCKTPNNMFSMPQLENCSYPEFMMQP", "ERKMIDFCCCAIHYGKFDAQEWHSEVSHVP", "WEWSGWDWANIVRSGQHWPGDCRFVNMQDC", "QTQDATVFWCCMQKIETYLYQWFLENYFRL", "DRPHEHSDSDTEAPYACITHSKQVDANDSQ", "TMCCVWDYFSVNEDRRLIQQNPYEVEHELC", "KQEDMQDLNAMLVSWLACIYEAIRWALGWI"], 7)
        self.t7 = SimpleBlast(["HSGVYLDDRDKPHFQDLNSFILDMSNNQRD", "RFIIQTPYVRWKQRPEFWYPYEFEAVQTNS", "MTNQNHEVCHFYMCIFMLFIHQKKSSPYLS", "GYLWHCKTPNNMFSMPQLENCSYPEFMMQP", "ERKMIDFCCCAIHYGKFDAQEWHSEVSHVP", "WEWSGWDWANIVRSGQHWPGDCRFVNMQDC", "QTQDATVFWCCMQKIETYLYQWFLENYFRL", "DRPHEHSDSDTEAPYACITHSKQVDANDSQ", "TMCCVWDYFSVNEDRRLIQQNPYEVEHELC", "KQEDMQDLNAMLVSWLACIYEAIRWALGWI"], 3)
        
    def test_get_hits(self):
        self.assertEqual(self.t1.get_hits('AATATGTTATATAATAATATTT', 'AATATATTTATA'), [(0, 0), (0, 12), (0, 15), (1, 1), (1, 8), (1, 10), (1, 13), (1, 16), (3, 1), (3, 8), (3, 10), (3, 13), (3, 16), (9, 1), (9, 8), (9, 10), (9, 13), (9, 16), (6, 19), (2, 2), (2, 7), (2, 9), (2, 17), (4, 2), (4, 7), (4, 9), (4, 17), (8, 2), (8, 7), (8, 9), (8, 17), (7, 6), (5, 18)])
        self.assertEqual(self.t2.get_hits('AATATGTTATATAATAATATTT', 'AATATATTTATA'), [(0, 0), (0, 12), (0, 15), (1, 1), (1, 8), (1, 10), (1, 13), (1, 16), (3, 1), (3, 8), (3, 10), (3, 13), (3, 16), (9, 1), (9, 8), (9, 10), (9, 13), (9, 16), (6, 19), (2, 2), (2, 7), (2, 9), (2, 17), (4, 2), (4, 7), (4, 9), (4, 17), (8, 2), (8, 7), (8, 9), (8, 17), (7, 6), (5, 18)])
        self.assertEqual(self.t3.get_hits('CCCGCCATGCACCGA', 'CCCGCC'), [(0, 0), (3, 3), (1, 1), (1, 11), (2, 2)])
        self.assertEqual(self.t4.get_hits('agactacgatcagttttttcagtacccccgagccca', 'acaatcagt'), [(3, 8), (6, 11), (6, 20), (4, 9), (4, 18), (5, 10), (5, 19)])
        self.assertEqual(self.t5.get_hits('agactacgatcagttttttcagtacccccgagccca', 'acaatcagtccca'), [(3, 8)])
        self.assertEqual(self.t7.get_hits('HSGVYLDDRDKPHFQDLNSFILDMSNNQRD', 'YLDCTQPWIAKVMPA'), [(0, 4)])
    
    def test_get_matches(self):
        self.assertEqual(self.t1.get_matches('AATATGTTATATAATAATATTT', 'AATATATTTATA'), {(0, 0): 7, (0, 12): 6, (0, 15): 6, (1, 1): 7, (1, 8): 9, (1, 10): 3, (1, 13): 6, (1, 16): 6, (3, 1): 5, (3, 8): 5, (3, 10): 9, (3, 13): 5, (3, 16): 7, (9, 1): 3, (9, 8): 3, (9, 10): 6, (9, 13): 2, (9, 16): 7, (6, 19): 7, (2, 2): 7, (2, 7): 1, (2, 9): 9, (2, 17): 6, (4, 2): 5, (4, 7): 6, (4, 9): 3, (4, 17): 7, (8, 2): 1, (8, 7): 2, (8, 9): 2, (8, 17): 6, (7, 6): 6, (5, 18): 3})
        self.assertEqual(self.t2.get_matches('AATATGTTATATAATAATATTT', 'AATATATTTATA'), {(0, 0): 7, (0, 12): 6, (0, 15): 6, (1, 1): 7, (1, 8): 9, (1, 10): 3, (1, 13): 6, (1, 16): 6, (3, 1): 5, (3, 8): 5, (3, 10): 9, (3, 13): 5, (3, 16): 7, (9, 1): 3, (9, 8): 3, (9, 10): 6, (9, 13): 2, (9, 16): 7, (6, 19): 7, (2, 2): 7, (2, 7): 1, (2, 9): 9, (2, 17): 6, (4, 2): 5, (4, 7): 6, (4, 9): 3, (4, 17): 7, (8, 2): 1, (8, 7): 2, (8, 9): 2, (8, 17): 6, (7, 6): 6, (5, 18): 3})
        self.assertEqual(self.t3.get_matches('CCCGCCATGCACCGA', 'CCCGCC'), {(0, 0): 6, (3, 3): 6, (1, 1): 6, (1, 11): 3, (2, 2): 6})
        self.assertEqual(self.t4.get_matches('agactacgatcagttttttcagtacccccgagccca', 'acaatcagt'), {(3, 8): 8, (6, 11): 8, (6, 20): 3, (4, 9): 8, (4, 18): 1, (5, 10): 8, (5, 19): 2})
        self.assertEqual(self.t5.get_matches('agactacgatcagttttttcagtacccccgagccca', 'acaatcagtccca'), {(3, 8): 8})
        self.assertEqual(self.t7.get_matches('HSGVYLDDRDKPHFQDLNSFILDMSNNQRD', 'YLDCTQPWIAKVMPA'), {(0, 4): 3})
        
    def test_best_alignment(self):
        self.assertEqual(self.t1.best_alignment('AATATATTTATA'), 'AATATGTTATATAATAATATTT')
        self.assertEqual(self.t1.best_alignment('aataTATTTata'), 'AATATGTTATATAATAATATTT')
        self.assertEqual(self.t2.best_alignment('AATATATTTATA'), 'AATATGTTATATAATAATATTT')
        self.assertEqual(self.t3.best_alignment('CCCGCC'), 'CCCGCCATGCACCGA')
        self.assertEqual(self.t4.best_alignment('acaatcagt'), 'GGGGGTACATGACGATACGATCAGT')
        self.assertEqual(self.t5.best_alignment('acaatcagtccca'), 'GGGGGTACATGACGATACGATCAGT')
        self.assertEqual(self.t6.best_alignment('YLDCTQPWIAKVMPA'), 'No sequence matches obtained')
        self.assertEqual(self.t7.best_alignment('YLDCTQPWIAKVMPA'), 'HSGVYLDDRDKPHFQDLNSFILDMSNNQRD')
        
if __name__ == '__main__':
    # begin the unittest.main()
    unittest.main(argv=['first-arg-is-ignored'], exit = False)