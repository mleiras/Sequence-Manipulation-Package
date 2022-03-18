from Sequence import Sequence
import math


class Motifs():
    def __init__(self, list_seq: list, **kwargs):
        '''Class that allows creating and displaying probabilistic profiles (PWM or PSSM) and
        other functions (? se ficar muito grande tira-se: calculate the probability of a sequence and
        calculate the most probable sequence of the profile created.)

        Parameters
        ----------
        list_seq : list
            Sequences of DNA, RNA or AMINO ACIS
        
        kwargs: dictionary
            May take as keys:
                profile -- pseudo
        '''
        self.list_seq = list_seq
        self.abc = self.validate_abc()
        self.profile = 'pwm'
        self.pseudo = 0
        self.n = 4
        
        if self.abc == 'DNA':
            self.abc_letters = "ACTG"
        elif self.abc == 'AMINO':
            self.abc_letters = "FLIMVSPTAYHQNKDECWRG_"
            self.n = 20
        elif self.abc == 'RNA':
            self.abc_letters = "ACGU"
        
        if 'profile' in kwargs:
            self.profile = kwargs['profile']
            
        if 'pseudo' in kwargs:
            self.pseudo = kwargs['pseudo']
        
    def validate_seq(self, seq: str) -> str:
        '''This method determines the type of the sequence: DNA, RNA or AMINO ACIDS
        It gives an error if the sequence is not recognized.

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS

        Returns
        -------
        str: str
            Type of the sequence ('DNA', 'RNA' or 'AMINO')
        '''
        test = Sequence(seq)
        return test.check
        
    
    def validate_abc(self) -> str:
        '''Through the Sequences Class, it validates if all the sequences of a given list are of the same type, be it DNA, RNA or PROTEIN (AMINO).
        If any sequence is not recognized or of the same type, an Error is displayed

        Returns
        -------
        str
            Sequence of DNA, RNA or AMINO ACIDS

        Raises
        ------
        TypeError
            None of the sequences are of the same type
        '''
        class_abc = Sequence(self.list_seq[0])
        abc = class_abc.check
        for i in self.list_seq[1:]:           
            valid = self.validate_seq(i)
            if valid != abc:
                raise TypeError("Sequences introduced are not of the same type!") 
        return abc
    
    
    def print_profile(self, profile: list) -> repr:
        '''Prints a probabilistic profile legibly

        Parameters
        ----------
        profile : list
            

        Returns
        -------
        repr
            A pretty way of showing the probabilistic profile
        '''
        bases = sorted(profile[0].keys())
        tab = [[f"{p[b]:-5.2f}" for b in bases] for p in profile]
        for p in zip(*([bases] + tab)):
            print(*p)
            
    
    def create_profile(self) -> list:
        '''Calculates the probabilistic PWM (Position Weighted Matrix) or PSSM (Position Specific Scoring Matrix) profile of a given list of **DNA** sequences.
        Default pseudocount of 0

        Returns
        -------
        matrix : list
            It returns list with a dictionary of the probabilistic profile (PWM or PSSM)

        Raises
        ------
        TypeError
            If the profile selected isn't pwm or pssm
        '''
        total = len(self.list_seq) + self.n*self.pseudo
        matrix = []
        if self.profile == 'pwm':
            for line in zip(*self.list_seq):
                matrix.append({k: (line.count(k)+self.pseudo)/total for k in self.abc_letters})
        elif self.profile == 'pssm':
            for line in zip(*self.list_seq):
                matrix.append({k: math.log2(((line.count(k)+self.pseudo)/total)/(1/self.n)) for k in self.abc_letters})
        else:
            raise TypeError("Type of profile invalid")   
        return matrix

    def prob_seq(self, seq: str, profile: list) -> float:
        '''Calculates and returns the probability of a given sequence by the associated profile
        
        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS
        profile : list
            List with a dictionary of the probabilistic profile (PWM or PSSM)

        Returns
        -------
        p : float
            A float number representative of the probability of the given sequence

        Raises
        ------
        TypeError
            If the profile selected isn't pwm or pssm
        '''
        assert len(seq) == len(profile), 'Sequence size does not match associated profile size'
        if self.profile == 'pwm':
            p = 1
            for b, c in zip(seq, profile):
                p *= c[b]
        elif self.profile == 'pssm':
            p = 0
            for b, c in zip(seq, profile):
                p += c[b]
        else:
            raise TypeError("Type of profile invalid")
        return round(p, 5)
    
    def seq_most_probable(self, seq: str, profile: list) -> str:
        '''Calculates the probability of each subsequence of a given sequence, and returns the subsequence with the highest probability according to the associated profile

        Parameters
        ----------
        seq : str
            Sequence of DNA, RNA or AMINO ACIDS
        profile : list
            List with a dictionary of the probabilistic profile (PWM or PSSM)

        Returns
        -------
        str
            Subsequence with the highest probability
        '''
        list_seq = [seq[I:I+len(profile)] for I in range(len(seq) - len(profile) + 1)]
        probs = [self.prob_seq(s, profile) for s in list_seq]
        score_max = max(probs)
        return [list_seq[I] for I,p in enumerate(probs) if p == score_max][0]

