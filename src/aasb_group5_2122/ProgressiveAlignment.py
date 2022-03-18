from Sequence import Sequence
from Alignment import Alignment
from Motifs import Motifs

class ProgressiveAlignment():
    def __init__(self, list_seqs: list, **kwargs):
        '''Multiple Alignment for a list of sequences (heuristic method?)

        Parameters
        ----------
        list_seqs : list
            List of the multiple alignments (DNA or Amino Acid)
        
        kwargs: dictionary
            May take as keys:
                profile -- pseudo -- g -- eqs -- diffs
        '''
        self.list_seqs = list_seqs
        self.profile = 'pwm'
        self.pseudo = 0
        self.dic = {'eqs': 4, 'diffs': -1, 'g': -4}
        self.abc = 'dna'
        
        if 'profile' in kwargs:
            self.profile = kwargs['profile']
        
        if 'pseudo' in kwargs:
            self.pseudo = kwargs['pseudo']
        
        if 'g' in kwargs:
            g = kwargs['g']
            self.dic['g'] = g
        
        if 'diffs' in kwargs:
            diffs = kwargs['diffs']
            self.dic['diffs'] = diffs
        
        if 'eqs' in kwargs:
            eqs = kwargs['eqs']
            self.dic['eqs'] = eqs
            
        if 'abc' in kwargs:
            self.abc = kwargs['abc']
            self.dic['abc'] = self.abc
    
    
    def consensus(self, s1: str, s2: str) -> str:
        '''Creates an instance of the Class Motif to create a probabilistic profile (PWM) and then creates a sequence consensus between two different sequences.

        Parameters
        ----------
        s1 : str
            String sequence of DNA or Amino Acids
        s2 : str
            String sequence of DNA or Amino Acids

        Returns
        -------
        str
            Sequence consensus
        '''
        seqs = [s1, s2]
        p_max = 0
        cons = ''
        motif = Motifs(seqs)
        pwm = motif.create_profile()
        for dic in pwm:
            m = max(*dic.values())
            p_max += m
            key = [k for k, v in dic.items() if v == m][0]
            cons += key
        return cons


    def progressive(self):
        '''Creates an instance of the Class Alinhament and it uses the method consensus for creating multiple alignments with the list of sequences.
        It repeats the multiple alignment with the final sequences and consensus.

        Returns
        -------
        list_align: list
            List of the multiple alignments
        '''
        s1, s2, *rest = self.list_seqs
        alinhamento = Alignment(s1, s2, g = self.dic['g'], diffs = self.dic['diffs'], eqs = self.dic['eqs'])
        if self.abc == 'dna':
            seq1, seq2 = alinhamento.simple_align(s1, s2)
        else:
            seq1, seq2 = alinhamento.simple_align(s1, s2, alfabeto = 'protein')
        cons = self.consensus(seq1,seq2)
        new_list = [seq1,seq2]
        for i in rest:
            alinhamento = Alignment(cons, i, g = self.dic['g'], diffs = self.dic['diffs'], eqs = self.dic['eqs'])
            if self.abc == 'dna':
                SEQ1, SEQ2 = alinhamento.simple_align(cons, i)
            else:
                SEQ1, SEQ2 = alinhamento.simple_align(cons, i, alfabeto = 'protein')
            cons = self.consensus(SEQ1,SEQ2)
            new_list.append(SEQ2)        
        list_align=[]
        for i in new_list:
            alinhamento = Alignment(cons, i, g = self.dic['g'], diffs = self.dic['diffs'], eqs = self.dic['eqs'])
            if self.abc == 'dna':
                cons, SEQ = alinhamento.simple_align(cons, i)
            else:
                cons, SEQ = alinhamento.simple_align(cons, i, alfabeto = 'protein')
            cons = self.consensus(cons,SEQ)
            list_align.append(SEQ)
        return list_align

