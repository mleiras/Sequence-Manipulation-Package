class Alignment:
    def __init__(self, seq1: str, seq2: str, **kwargs):
        '''
        Creates two matrix: first with the values in each position and second with the movements to be done on the alignment
        
        Parameters
        ------------
        seq1, seq2: str
            Sequences of which the matrix will be constructed
            
        kwargs: dictionary
            May take as keys:
                blosum62 -- local -- match_pairs -- file -- g -- eqs -- diffs -- alphabet
        '''
        self.seq1 = seq1
        self.seq2 = seq2
        blosum62, self.local, match_pairs, file = (False, False, False, False)
        self.dic = {'eqs': None, 'diffs': None, 'alphabet': None, 'g': None}
        
        if 'match_pairs' in kwargs:
            match_pairs = kwargs['match_pairs']
            eqs = 1
            diffs = 0
            g = -8
            alphabet = "ACGT"
            self.dic["eqs"], self.dic["diffs"], self.dic["alphabet"], self.dic['g'] = (eqs, diffs, alphabet, g)
            
        if 'g' in kwargs:
            g = kwargs['g']
            self.dic['g'] = g
            
        if 'eqs' in kwargs:
            eqs = kwargs['eqs']
            self.dic['eqs'] = eqs
            
        if 'diffs' in kwargs:
            diffs = kwargs['diffs']
            self.dic['diffs'] = diffs
            
        if 'alphabet' in kwargs:
            alphabet = kwargs['alphabet']
            self.dic['alphabet'] = alphabet
        
        if 'blosum62' in kwargs:
            blosum62 = kwargs['blosum62']
            self.dic = self.dic_blosum62()
        
        if 'local' in kwargs:
            self.local = kwargs['local']

        if 'file' in kwargs:
            self.dic = self.dic_blosum(kwargs['file'])
            file = True
        
        if blosum62 == True and match_pairs == True:
            raise Exception("Blosum62 and Match Pairs condition are not viable if both are True!")

        elif blosum62 == True and file == True:
            raise Exception("Blosum62 and File matrix are not viable when both are used as input!")
        
        self.matriz1 = self.matriz_0(seq1, seq2, g, self.local)
        self.matriz2 = self.matriz_0(seq1, seq2, g, self.local)
        
        for i in range(len(self.matriz1)):
            for j in range(len(self.matriz1[0])):
                if i >= 2 and j >= 2:
                    if blosum62 is True or file is True:
                        x = self.dic[self.matriz1[i][0]][self.matriz1[0][j]] + self.matriz1[i - 1][j - 1]
                    else:
                        if self.matriz1[0][j] == self.matriz1[i][0]:
                            x = eqs + self.matriz1[i - 1][j - 1]
                        else:
                            x = diffs + self.matriz1[i - 1][j - 1]
                            
                    if self.local is True:
                        self.matriz1[i][j] = max(x, self.matriz1[i][j - 1] + g, self.matriz1[i - 1][j] + g, 0)
                        self.matriz2[i][j] = self.verify_side(self.matriz1[i][j], x, self.matriz1[i][j - 1] + g, self.matriz1[i - 1][j] + g)
                    else:
                        self.matriz1[i][j] = max(x, self.matriz1[i][j - 1] + g, self.matriz1[i - 1][j] + g)
                        self.matriz2[i][j] = self.verify_side(self.matriz1[i][j], x, self.matriz1[i][j - 1] + g, self.matriz1[i - 1][j] + g)
    
    def dic_blosum(self, file: str) -> dict:
        '''
        Read a matrix of values according the file given by the user
        The file must be in the same directory as the python executer and needs to follow the aspects of a matrix
        
        Parameters
        ------------
        file: str
            The file location of the matrix with the characters values
        
        Returns
        ------------
        dic: dict
            Dictionary with the all the characters and the corresponding values
        '''
        f = open(file, 'r')
        linhas = f.read()

        dic = {}
        headers, *mat = [linha.split() for linha in linhas.splitlines()]
        for linha in mat:
            letra, *resto = linha
            dic[letra] = {}
            for outra, val in zip(headers, resto):
                dic[letra][outra] = int(val)
        
        f.close()
        return dic

    def dic_blosum62(self) -> dict:
        '''
        Construction of a dictionary with the values blosum62 of the characters
        
        Returns 
        ------------
        dic: dict
            Dictionary with all the characters and the corresponding values
        '''
        linhas = '''   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
        A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
        R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
        N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
        D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
        C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
        Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
        E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
        G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
        H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
        I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
        L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
        K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
        M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
        F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
        P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
        S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
        T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
        W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
        Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
        V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
        B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
        Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
        X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
        * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1'''

        dic = {}
        headers, *mat = [linha.split() for linha in linhas.splitlines()]
        for linha in mat:
            letra, *resto = linha
            dic[letra] = {}
            for outra, val in zip(headers, resto):
                dic[letra][outra] = int(val)
        
        return dic

    def __str__(self) -> str:
        '''
        Writing of arrays obtained at object initialization
        
        Returns
        ------------
        string: str
            String with the matrices constructed at object initialization
            
        '''
        s1 = [[str(e) for e in row] for row in self.matriz1]
        lens1 = [max(map(len, col)) for col in zip(*s1)]
        fmt1 = ' '.join('{{:{}}}'.format(x) for x in lens1)
        table1 = [fmt1.format(*row) for row in s1]
        string = '\n'.join(table1)
        string += '\n\n'
    
        s2 = [[str(e) for e in row] for row in self.matriz2]
        lens2 = [max(map(len, col)) for col in zip(*s2)]
        fmt2 = ' '.join('{{:{}}}'.format(x) for x in lens2)
        table2 = [fmt2.format(*row) for row in s2]
        string += '\n'.join(table2)

        return string

    def matriz_0(self, seq1: str, seq2: str, g: int, l: bool) -> list:
        '''
        Method that obtains an 'empty' matrix in order to be altered and construct the desired ones
        Serves as a template
        Called in the initialization of the object
        
        Parameters
        ------------
        seq1, seq2: str
            Sequences of which the matrix will be constructed
        
        g: int
            Corresponds to the spacing value
        
        l: bool
            Verifys the use of local alignment in case of True argument
            
        Returns
        ------------
        matriz: list
            Matrix template to manipulate
        '''
        matriz = [[0 for c in range(len(seq2) + 2)] for r in range(len(seq1) + 2)]
        matriz[0][0] = ' '
        matriz[0][1] = '-'
        matriz[1][0] = '-'
        
        for i in range(len(seq1)):
            matriz[i + 2][0] = seq1[i]
        for i in range(len(seq2)):
            matriz[0][i + 2] = seq2[i]
        
        for i in range(len(matriz)):
            if i >= 1:
                if l is True:
                    matriz[i][1] = 0
                else:
                    matriz[i][1] = g * (i - 1)

        for i in range(len(matriz[0])):
            if i >= 1:
                if l is True:
                    matriz[1][i] = 0
                else:
                    matriz[1][i] = g * (i - 1)
                
        return matriz    

    def verify_side(self, valor: int, x: int, y: int, z: int) -> str:
        '''
        Method that returns the movement of the given arguments
        Called in the initialization of the object
        
        Parameters
        ------------
        valor: int
            Base value
        
        x, y, z: int
            Value to be compared to the value and return the specific movement
        
        Returns
        ------------
        Movement: str
            List of possible movements: (DE, DC, EC, D, E, C, DEC)
        '''
        if valor == x and valor == y: return 'DE'
        elif valor == x and valor == z: return 'DC'
        elif valor == y and valor == z: return 'EC'
        elif valor == x: return 'D'
        elif valor == y: return 'E'
        elif valor == z: return 'C'
        else: return 'DEC'

    def alignScore(self) -> str:
        '''
        Intermediate method of the "_alignScore" function
        
        Returns
        ------------
        str
            Pretty string with the type of alignment and score
        '''
        alignment, score = self._alignScore()
        return f"Best {alignment} alignment score: {score}"
    
    def _alignScore(self) -> tuple:
        '''
        Method that gives the global or local score of the constructed matrix for the given sequences
        
        Returns
        ------------
        int
            If local alignment the returned value corresponds to the max of the matrix 
            On the other hand, global alignment returns the extreme under-right score 
        '''
        max = 0
        if self.local is True:
            for i in range(len(self.matriz1)):
                for j, value in enumerate(self.matriz1[i]):
                    if type(value) == int and value > max:
                        max = value
                        x = i
                        y = j
        
        if self.local is True:
            return 'local', max
        else:
            return 'global', self.matriz1[len(self.matriz1) - 1][len(self.matriz1[0]) - 1]
    
    def getScore(self, a: str, b: str) -> str:
        '''
        Intermediate method of the "_getScore" function
        
        Returns
        ------------
        str
            Pretty string with the characters and corresponding value
        '''
        return f"The value between characters {a} and {b} is equal to {self._getScore(a, b)}"

    def _getScore(self, a: str, b: str) -> int:
        '''
        Returns the value corresponded of the individual caracteres inputted by the user
        
        Parameters
        ------------
        a, b: str
            Characters intended to verify the value
        
        Returns
        ------------
        int
            Value of the inputted characters
        '''
        if 'alphabet' in self.dic.keys():
            assert a in self.dic["alphabet"] and b in self.dic["alphabet"], "Caracters introduced are absent in alphabet!"
            if a == b:
                return self.dic["eqs"]
            elif a != b:
                return self.dic["diffs"]
        else:
            return self.dic[a][b]
        
    def align_seq(self, ties = False) -> str:
        '''
        Method to get all the possible alignments between the given sequences if parameter ties is True.
        
        Parameters
        ------------
        ties: boolean
            Takes False argument as default and indicates if we should assume the ties in the matrix or not
            If False the diagonal path is the preferred
            
        Returns
        ------------
        string: str
            Pretty string with the single alignment for the diagonal movement or all alignments in case ties equals True
        '''
        lseq1, lseq2 = self._align_seq()
        
        string = ""
        if ties is True:
            for i in range(len(lseq1)):
                if i == len(lseq1) - 1:
                    string += f"Alignment {i+1}:\n{lseq1[i][::-1]}\n{lseq2[i][::-1]}"
                else:
                    string += f"Alignment {i+1}:\n{lseq1[i][::-1]}\n{lseq2[i][::-1]}\n\n"
        elif ties is False:
            string += f"Alignment:\n{lseq1[0][::-1]}\n{lseq2[0][::-1]}"
        else:
            raise Exception("Ties can only take True or False argument!")
            
        return string
            
    def _align_seq(self) -> tuple:
        '''
        Method that constructs all the alignments and stores in the variables "self.str_seq1" and "self.str_seq2"
        
        Returns
        ------------
        self.str_seq1, self.str_seq2: list
            List of all the alignments being each position of the list the alignment of the other
            For example -> self.str_seq1[0] is the alignment of self.str_seq2[0]
        '''
        string = ""
        self.str_seq1 = [""]
        self.str_seq2 = [""]
        max = 0
        if self.local is True:
            for i in range(len(self.matriz1)):
                for j, value in enumerate(self.matriz1[i]):
                    if type(value) == int and value > max:
                        max = value
                        x = i
                        y = j
        else:
            x = len(self.matriz2) - 1
            y = len(self.matriz2[0]) - 1
        
        self.get_mov(0, x, y)
        return self.str_seq1, self.str_seq2      
    
    def get_mov(self, i: int, x: int, y: int):
        '''
        Recursive method that continuously implements in the lists "self.str_seq1" and "self.str_seq2" the originated sequences from the respective movements pathways
        Called in the "_align_seq" function
        
        Parameters
        ------------
        i: int
            Each position that will correspond to the position to be incremented on the lists
        
        x, y: int
            Coordinates to begin the recursive function in the matrix and correctly implement in the lists the alignments
        '''
        if self.local is True and self.matriz1[x][y] == 0:
            return 1
        
        if (x == 1 and y == 1):
            return 1

        elif str(self.matriz2[x][y]).isalpha() is False:
            self.str_seq1[i] += self.matriz1[x][0]
            self.str_seq2[i] += self.matriz1[0][y]
            return 1
        
        if self.matriz2[x][y] == "D":
            self.str_seq1[i] += self.matriz1[x][0]
            self.str_seq2[i] += self.matriz1[0][y]
            self.get_mov(i, x - 1, y - 1)
            
        elif self.matriz2[x][y] == "E":
            self.str_seq1[i] += "-"
            self.str_seq2[i] += self.matriz1[0][y]
            self.get_mov(i, x, y - 1)
            
        elif self.matriz2[x][y] == "C":
            self.str_seq1[i] += self.matriz1[x][0]
            self.str_seq2[i] += "-"
            self.get_mov(i, x - 1, y)
            
        elif self.matriz2[x][y] == "DE":
            self.str_seq1.append(self.str_seq1[i])
            self.str_seq2.append(self.str_seq2[i])
            # Diagonal path
            self.str_seq1[i] += self.matriz1[x][0]
            self.str_seq2[i] += self.matriz1[0][y]
            self.get_mov(i, x - 1, y - 1)
            # Left path
            self.str_seq1[i + 1] += "-"
            self.str_seq2[i + 1] += self.matriz1[0][y]
            self.get_mov(i + 1, x, y - 1)
        
        elif self.matriz2[x][y] == "DC":
            self.str_seq1.append(self.str_seq1[i])
            self.str_seq2.append(self.str_seq2[i])
            # Diagonal path
            self.str_seq1[i] += self.matriz1[x][0]
            self.str_seq2[i] += self.matriz1[0][y]
            self.get_mov(i, x - 1, y - 1)
            # Uppwards path
            self.str_seq1[i + 1] += self.matriz1[x][0]
            self.str_seq2[i + 1] += "-"
            self.get_mov(i + 1, x - 1, y)
        
        elif self.matriz2[x][y] == "EC":
            self.str_seq1.append(self.str_seq1[i])
            self.str_seq2.append(self.str_seq2[i])
            # Left path
            self.str_seq1[i] += "-"
            self.str_seq2[i] += self.matriz1[0][y]
            self.get_mov(i, x, y - 1)
            # Uppwards path
            self.str_seq1[i + 1] += self.matriz1[x][0]
            self.str_seq2[i + 1] += "-"
            self.get_mov(i + 1, x - 1, y)
            
        elif self.matriz2[x][y] == "DEC":
            self.str_seq1.append(self.str_seq1[i])
            self.str_seq2.append(self.str_seq2[i])
            self.str_seq1.append(self.str_seq1[i])
            self.str_seq2.append(self.str_seq2[i])
            # Diagonal path
            self.str_seq1[i] += self.matriz1[x][0]
            self.str_seq2[i] += self.matriz1[0][y]
            self.get_mov(i, x - 1, y - 1)
            # Left path
            self.str_seq1[i + 1] += "-"
            self.str_seq2[i + 1] += self.matriz1[0][y]
            self.get_mov(i + 1, x, y - 1)
            # Uppwards path
            self.str_seq1[i + 2] += self.matriz1[x][0]
            self.str_seq2[i + 2] += "-"
            self.get_mov(i + 2, x - 1, y)
        
        return 1
    
    def score_s(self, x1: str, x2: str) -> int:
        '''Intermediate method of the "simple_align" function when alfabeto = "dna" that returns score between two elements

        Parameters
        ----------
        x1 : str
            element of sequence for scoring
        x2 : str
            element of sequence for scoring

        Returns
        -------
        int
            Score between two elements
        '''
        if x1 == x2:
            return self.dic['eqs']
        else: return self.dic['diffs']
    
    def score_blosum62(self, x1: str, x2: str) -> int:
        '''Intermediate method of the "simple_align" function when alfabeto = "protein" that returns score between two elements based on the blosum62 matrix

        Parameters
        ----------
        x1 : str
            element of sequence for scoring
        x2 : str
            element of sequence for scoring

        Returns
        -------
        int
            Score between two elements
        '''
        if x1 == '-' or x2 == '-':
            return self.dic['g']
        else:
            dic = self.dic_blosum62()
            return dic[x1][x2]
            
    def maxMat(self, mat: list) -> tuple:
        '''Intermediate method of the "simple_align" function when tipo = "local"

        Parameters
        ----------
        mat : list
            matrix of the alignment

        Returns
        -------
        maxval, maxrow, maxcol : tuple
            Tuple of the max value of the local alignment and the coordenates of max value
        '''
        maxval = mat[0][0]
        maxrow = 0
        maxcol = 0
        for i in range(0, len(mat)):
            for j in range(0, len(mat[i])):
                if mat[i][j]>maxval:
                    maxval = mat[i][j]
                    maxrow = i
                    maxcol = j
        return (maxval, maxrow, maxcol)

    def simple_align(self, s1: str, s2: str, alfabeto: str = "protein", tipo: str = "global") -> tuple:
        '''Method that constructs the alignments between two sequences

        Parameters
        ----------
        s1 : str 
            String of DNA or AMINO ACIDS (protein)
        s2 : [type]
            String of DNA or AMINO ACIDS (protein)
        alfabeto : str, optional
            Type of alphabet: DNA or AMINO ACIDS (protein), by default "protein"
        tipo : str, optional
            type of alignment (global or local), by default "global"

        Returns
        -------
        S1, S2 : tuple
            Two alignment strings
        '''
        mat = [[0 for c in range(len(s2)+1)] for L in range(len(s1)+1)]
        tr = [[" " for c in range(len(s2)+1)] for L in range(len(s1)+1)]

        if tipo == "global":  
            for L in range(len(s1)):
                mat[L+1][0] = mat[L][0] + self.dic['g']
                tr[L+1][0] = "C"
            for C in range(len(s2)):
                mat[0][C+1] = mat[0][C] + self.dic['g']
                tr[0][C+1] = "E"

            for L,x1 in enumerate(s1):
                for C,x2 in enumerate(s2):
                    if alfabeto == "protein":
                        valor = [
                            mat[L][C] + self.score_blosum62(x1, x2),
                            mat[L+1][C] + self.dic['g'], 
                            mat[L][C+1] + self.dic['g'] ]
                    elif alfabeto == "dna":
                        valor = [
                            mat[L][C] + self.score_blosum62(x1,x2), 
                            mat[L+1][C] + self.dic['g'], 
                            mat[L][C+1] + self.dic['g'] ]
                    else: break
                    direcoes = "D E C".split()
                    mat[L+1][C+1] = max(*valor)
                    tr[L+1][C+1] = direcoes[valor.index(mat[L+1][C+1])]
                    
        elif tipo == "local":
            for L,x1 in enumerate(s1):
                for C,x2 in enumerate(s2):
                    if alfabeto == "protein":
                        valor = [
                            mat[L][C] + self.score_s(x1,x2), 
                            mat[L+1][C] + self.dic['g'], 
                            mat[L][C+1] + self.dic['g'],
                            0]
                    elif alfabeto == "dna":
                        valor = [
                            mat[L][C] + self.score_s(x1,x2), 
                            mat[L+1][C] + self.dic['g'], 
                            mat[L][C+1] + self.dic['g'],
                            0]
                    else: break
                    direcoes = "D E C".split()
                    mat[L+1][C+1] = max(*valor)
                    if mat[L+1][C+1] == 0:
                        tr[L+1][C+1] = " "
                    else:
                        tr[L+1][C+1] = direcoes[valor.index(mat[L+1][C+1])]
        else: return False

        L = len(s1)
        C = len(s2)
        S1, S2 = "",""
        
        if tipo == "global":
            while L> 0 or C> 0:
                if tr[L][C] == "D":
                    S1 = s1[L-1] + S1
                    S2 = s2[C-1] + S2
                    L -= 1
                    C -= 1
                elif tr[L][C] == "E":
                    S1 = "-" + S1
                    S2 = s2[C-1] + S2
                    C -= 1
                elif tr[L][C] == "C":
                    S2 ="-" + S2
                    S1 = s1[L-1] + S1
                    L -= 1
                else:
                    break
        
            return S1, S2
        
        elif tipo == "local":
            S1, S2 = "",""
            maxv, L, C = self.maxMat(mat)
            
            while mat[L][C] != 0:
                if tr[L][C] == "D":
                    S1 = s1[L-1] + S1
                    S2 = s2[C-1] + S2
                    L -= 1
                    C -= 1
                elif tr[L][C] == "E":
                    S1 = "-" + S1
                    S2 = s2[C-1] + S2
                    C -= 1
                elif tr[L][C] == "C":
                    S2 = "-" + S2
                    S1 = s1[L-1] + S1
                    L -= 1
                else:
                    break
            return S1, S2
        
