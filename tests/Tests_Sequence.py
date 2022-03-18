# -*- coding: utf-8 -*-

import unittest
from Sequence import Sequence
#from Sequence import percentage, vow_cons, transcription, comp_inverse,translation, cDNA, get_orfs, get_codons, get_all_prots


class TestSequences (unittest.TestCase):
    
    def setUp(self):
        self.t1 = Sequence("ACGT")
        self.t2 = Sequence("ACGU")
        self.t3 = Sequence("FLIMVSPTAY_HQNKDECWRG")
        self.t4 = Sequence('HSIDN')
        self.t5 = Sequence('AUGCUA')
        self.t6 = Sequence('CAgatgattt')
        self.t7 = Sequence('ACCAGAGCGG')
        self.t8 = Sequence('ATTTAATTACAAGTCTTCAGAATGCCAGAGATATACAGGATCTAACCA')
        self.t9 = Sequence('tttttttt')
        self.t10 = Sequence('CAgaugauuu')
        self.t11 = Sequence('AUUUAAUUACAAGUCUUCAGAAUGCCAGAgauau')
        self.t12 = Sequence('ATGTTGGGCATGATCAAGAACTCGCTGTTCGGAAGCGTAGAGACGTGCCTTGGCAGGTCCTA')
        self.t13 = Sequence('ATGTTAGAGTTATTAAAAAGTCTGGTATTCGCCGTAATCATGGTACCTGTCGTGATGGCCATCATCCTGGGTCTGATTTACGGTCTTGGTGAAGTATTCAACATCTTTTCTGGTGTTGGTAAAAAAGACCAGCCCGGACAAAATCATTGA')

        
    def test_percentage(self):
        self.assertEqual(self.t1.percentage(), {'A': 25, 'C':25 , 'G':25, 'T':25 } )
        self.assertEqual(self.t2.percentage(), {'A': 25, 'C': 25, 'G': 25, 'U': 25} )
        self.assertEqual(self.t3.percentage(), {'F': 4.8,'L': 4.8,'I': 4.8,'M': 4.8,'V': 4.8,'S': 4.8,'P': 4.8, 'T': 4.8,'A': 4.8,'Y': 4.8,'_': 4.8,'H': 4.8,'Q': 4.8,'N': 4.8,'K': 4.8,'D': 4.8,'E': 4.8,'C': 4.8,'W': 4.8,'R': 4.8,'G': 4.8})

        
    def test_transcription(self):
        self.assertRaises(Exception, Sequence.transcription, '')
        self.assertRaises(Exception, Sequence.transcription, 'HSIDN')
        self.assertRaises(Exception, Sequence.transcription, 'hyeod')
        self.assertRaises(Exception, Sequence.transcription, 'ATGCTX')
        self.assertRaises(Exception, Sequence.transcription, 'AUGCUA')
        self.assertRaises(Exception, Sequence.transcription, '15485')
        self.assertEqual(self.t1.transcription(), 'ACGU')
        self.assertEqual(self.t6.transcription(), 'CAGAUGAUUU')
        self.assertEqual(self.t7.transcription(), 'ACCAGAGCGG')
        self.assertEqual(self.t8.transcription(), 'AUUUAAUUACAAGUCUUCAGAAUGCCAGAGAUAUACAGGAUCUAACCA')
        self.assertEqual(self.t9.transcription(), 'UUUUUUUU')
     
    
    def test_comp_inverse (self):
        self.assertRaises(Exception, Sequence.comp_inverse, '')
        self.assertRaises(Exception, Sequence.comp_inverse, 'HSIDN')
        self.assertRaises(Exception, Sequence.comp_inverse, 'hyeod')
        self.assertRaises(Exception, Sequence.comp_inverse, 'ATGCTX')
        self.assertRaises(Exception, Sequence.comp_inverse, 'AUGCUA')
        self.assertRaises(Exception, Sequence.comp_inverse, '15485')
        self.assertTrue(self.t1.comp_inverse().isupper())
        self.assertEqual(self.t6.comp_inverse(), 'AAATCATCTG')
        self.assertEqual(self.t8.comp_inverse(), 'TGGTTAGATCCTGTATATCTCTGGCATTCTGAAGACTTGTAATTAAAT')
    
    
    def test_translation(self):
        self.assertRaises(Exception, Sequence.translation, '')
        self.assertRaises(Exception, Sequence.translation, 'HSIDN')
        self.assertRaises(Exception, Sequence.translation, 'hyeod')
        self.assertRaises(Exception, Sequence.translation, 'ATGCTX')
        self.assertRaises(Exception, Sequence.translation, 'AUGCUA')
        self.assertRaises(Exception, Sequence.translation, '15485')
        self.assertTrue(self.t1.translation().isupper())
        self.assertEqual(self.t1.translation(), 'T')
        self.assertEqual(self.t8.translation(), 'I_LQVFRMPEIYRI_P')
    
    
    def test_cDNA(self):
        self.assertEqual(self.t10.cDNA(), 'GTCTACTAAA')
        self.assertEqual(self.t11.cDNA(),'TAAATTAATGTTCAGAAGTCTTACGGTCTCTATA')

   
    def test_get_codons(self): 
        self.assertEqual(self.t8.get_codons(), ['ATT', 'TAA', 'TTA', 'CAA', 'GTC', 'TTC', 'AGA', 'ATG', 'CCA', 'GAG', 'ATA', 'TAC', 'AGG', 'ATC', 'TAA', 'CCA'])
        self.assertEqual(self.t12.get_codons(), ['ATG', 'TTG', 'GGC', 'ATG', 'ATC', 'AAG', 'AAC', 'TCG', 'CTG', 'TTC', 'GGA', 'AGC', 'GTA', 'GAG', 'ACG', 'TGC', 'CTT', 'GGC', 'AGG', 'TCC'])    
    
    
    def test_get_orfs(self):
        self.assertEqual(self.t8.get_orfs(), {'ORF +1': ['ATT', 'TAA', 'TTA', 'CAA', 'GTC', 'TTC', 'AGA', 'ATG', 'CCA', 'GAG', 'ATA', 'TAC', 'AGG', 'ATC', 'TAA', 'CCA'], 'ORF +2': ['TTT', 'AAT', 'TAC', 'AAG', 'TCT', 'TCA', 'GAA', 'TGC', 'CAG', 'AGA', 'TAT', 'ACA', 'GGA', 'TCT', 'AAC'], 'ORF +3': ['TTA', 'ATT', 'ACA', 'AGT', 'CTT', 'CAG', 'AAT', 'GCC', 'AGA', 'GAT', 'ATA', 'CAG', 'GAT', 'CTA', 'ACC'], 'ORF -1': ['TGG', 'TTA', 'GAT', 'CCT', 'GTA', 'TAT', 'CTC', 'TGG', 'CAT', 'TCT', 'GAA', 'GAC', 'TTG', 'TAA', 'TTA', 'AAT'], 'ORF -2': ['GGT', 'TAG', 'ATC', 'CTG', 'TAT', 'ATC', 'TCT', 'GGC', 'ATT', 'CTG', 'AAG', 'ACT', 'TGT', 'AAT', 'TAA'], 'ORF -3': ['GTT', 'AGA', 'TCC', 'TGT', 'ATA', 'TCT', 'CTG', 'GCA', 'TTC', 'TGA', 'AGA', 'CTT', 'GTA', 'ATT', 'AAA']})
        self.assertEqual(self.t12.get_orfs(), {'ORF +1': ['ATG', 'TTG', 'GGC', 'ATG', 'ATC', 'AAG', 'AAC', 'TCG', 'CTG', 'TTC', 'GGA', 'AGC', 'GTA', 'GAG', 'ACG', 'TGC', 'CTT', 'GGC', 'AGG', 'TCC'], 'ORF +2': ['TGT', 'TGG', 'GCA', 'TGA', 'TCA', 'AGA', 'ACT', 'CGC', 'TGT', 'TCG', 'GAA', 'GCG', 'TAG', 'AGA', 'CGT', 'GCC', 'TTG', 'GCA', 'GGT', 'CCT'], 'ORF +3': ['GTT', 'GGG', 'CAT', 'GAT', 'CAA', 'GAA', 'CTC', 'GCT', 'GTT', 'CGG', 'AAG', 'CGT', 'AGA', 'GAC', 'GTG', 'CCT', 'TGG', 'CAG', 'GTC', 'CTA'], 'ORF -1': ['TAG', 'GAC', 'CTG', 'CCA', 'AGG', 'CAC', 'GTC', 'TCT', 'ACG', 'CTT', 'CCG', 'AAC', 'AGC', 'GAG', 'TTC', 'TTG', 'ATC', 'ATG', 'CCC', 'AAC'], 'ORF -2': ['AGG', 'ACC', 'TGC', 'CAA', 'GGC', 'ACG', 'TCT', 'CTA', 'CGC', 'TTC', 'CGA', 'ACA', 'GCG', 'AGT', 'TCT', 'TGA', 'TCA', 'TGC', 'CCA', 'ACA'], 'ORF -3': ['GGA', 'CCT', 'GCC', 'AAG', 'GCA', 'CGT', 'CTC', 'TAC', 'GCT', 'TCC', 'GAA', 'CAG', 'CGA', 'GTT', 'CTT', 'GAT', 'CAT', 'GCC', 'CAA', 'CAT']})        
    
    
    def test_get_all_prots(self):
        self.assertEqual(self.t1.get_all_prots(), [])       
        self.assertEqual(self.t8.get_all_prots(), ['MPEIYRI_'])
        self.assertEqual(self.t13.get_all_prots(),['MLELLKSLVFAVIMVPVVMAIILGLIYGLGEVFNIFSGVGKKDQPGQNH_'] )

    
if __name__ == '__main__':
   unittest.main()
    













    
    


















if __name__ == '__main__':
   unittest.main()