class Blast:
    '''
    Class to perform a Blast between a sequence query, a given sequence and the respective size of the query hits aiming to obtain in the sequence
    '''
    def __init__(self, seq: str, query: str, w: int):
        '''
        Blast Class simplified that searches good local alignments of a sequence within another sequence based on perfect matches

        Parameters
        ----------
        seq : str
            Sequence
            
        query : str
            Sequence query  
            
        w : int
            It determines the size of the subsequences to search in the seq
        '''
        assert (type(seq) == str) and (type(query) == str) and (type(w) == int) and (w > 0) and (w < len(query)), "Query and sequence must be a string and W an integer. W needs to be positive and lower than the size of the query"
        self.seq = seq.upper()
        self.query = query.upper()
        self.w = w
        
    def _get_index(self, q: str, query: str) -> list:
        '''
        Function that searches in the query the indexes of the query subsequences: "q". Auxiliar function of query_map

        Parameters
        ----------
        q: str
            Query subsequence aiming to search in the query
        
        query : str
            Query sequence serving as a template to find the indexes
            
        Returns
        ----------
        index: list
            q string positions in the query 
        '''
        assert type(q) == str, "Sequence must be a string"
        q = q.upper()
        
        index = []
        for ind in range(len(query)):
            if q == query[ind:ind + len(q)]: index.append(ind)
        return index
    
    
    def query_map(self, query = None, w = None) -> dict:
        '''
        Construction of a dictionary where the keys are the subsequences of the query with size W and the keys are the indexes 

        Parameters
        ----------
        query: str
            Query sequence serving as a template to find the indexes
        
        w: int
            Size each subsequence must have
        
        Returns
        ----------
        qm: dict
            Query subsequences with the corresponding indexes
        '''        
        if query is None: query = self.query
        if w is None: w = self.w
        
        assert (type(query) == str) and (type(w) == int) and (w > 0) and (w < len(query)), "Query must be a string and W an integer. W needs to be positive and lower than the size of the query"
        query = query.upper()
        
        qm = {}
        query_div = []
        for j in range(w):
            query_div.append([query[i: i + w] for i in range(j, len(query), w)])

        querys = []
        for i in query_div:
            for j in i:
                if j not in querys and len(j) == w: querys.append(j)

        for q in querys:
            indices = self._get_index(q, query)
            qm[q] = indices

        return qm
    

    def hits(self, qm = None, seq = None) -> list:
        '''
        Construction of a tuple list with the coordinates of the query subsequences in the query and in the given sequence
        
        Parameters
        ----------
        qm : dict
            Query subsequences with the corresponding indexes in the query
            
        seq: str
            Sequence serving as a template to search the indexes

        Returns
        ----------
        hits: list
            Tuple list with the first member of the tuple being the coordinate in the query and the second in the given sequence
        '''        
        if seq is None:
            seq = self.seq
        
        if qm is None:
            qm = self.query_map()
        
        assert type(qm) == dict and all(type(k) == str for k in qm.keys()) and \
               all(type(v) in [list,tuple] for v in qm.values()) and type(seq) == str, "Query map (qm) must be a dictionary and sequence a string"
        
        for q in qm:
            assert all(isinstance(x,int) for x in qm[q]) and all( x >= 0 for x in qm[q]), "Query map (qm) values must only contain integeres and positives"
        seq = seq.upper()
        
        hits = []
        for key, value in qm.items():
            for m in value:
                for i in range(len(seq)):
                    if key == seq[i:i + len(key)]: hits.append((m, i))
        return hits
    
    
    def extend_hit(self, index: list, query = None, seq = None, w = None) -> tuple:
        '''
        Construction of a tuple with the index of the beginning of the extended hit in the query, in the sequence, the respective size, and number of matches

        Parameters
        ----------
        index : list
            Specific hit to be extended in each direction
        
        query: str
            Sequence query
            
        seq: str
            Sequence
        
        w: int
            It determines the size of the subsequences to search in the seq

        Returns
        ----------
        ind_query_left, ind_seq_left, len(seq_final), cont_final: tuple
            Returns respectively index of the extended hit in the query, in the sequence, size of the hit, and the respective number of matches
        '''        
        if query is None: query = self.query
        if seq is None: seq = self.seq
        if w is None: w = self.w
        
        assert type(index) in [list, tuple] and type(seq) == str and type(query) == str and type(w) == int and w > 0 and w < len(query), "Unfeasible Parameters!"
        query = query.upper()
        seq = seq.upper()
        
        ind_query_left = index[0]
        ind_seq_left = index[1]
        ind_query_right = index[0]
        ind_seq_right = index[1]
        cont_left = 0
        cont_left_iguais = 0
        cont_right = 0
        cont_right_iguais = 0
        score = 1
        while score >= 0.5:
            if ind_seq_left == 0 or ind_query_left == 0: break
            cont_left += 1
            if (ind_seq_left - w) < 0 or (ind_query_left - w) < 0:
                ind_seq_left -= 1
                ind_query_left -= 1
            else:
                ind_seq_left -= w
                ind_query_left -= w

            if seq[ind_seq_left] == query[ind_query_left]: cont_left_iguais += 1

            if cont_left > 1: score = cont_left_iguais / cont_left

        while score >= 0.5:
            if ind_seq_right == len(seq) - 1 or ind_query_right == len(query) - 1: break
            cont_right += 1
            if (ind_seq_right + w) > len(seq) - 1 or (ind_query_right + w) > len(query) - 1:
                ind_seq_right += 1
                ind_query_right += 1
            else:
                ind_seq_right += w
                ind_query_right += w

            if seq[ind_seq_right] == query[ind_query_right]: cont_right_iguais += 1

            if cont_right > 1: score = cont_right_iguais / cont_right

        seq_final = seq[ind_seq_left:ind_seq_right + 1]
        query_final = query[ind_query_left:ind_query_right + 1]

        cont_final = 0
        for i in range(len(seq_final)):
            if seq_final[i] == query_final[i]: cont_final += 1
            
        return ind_query_left, ind_seq_left, len(seq_final), cont_final
    
    
    def best_hit(self, hits = None) -> tuple:
        '''
        Returns the best extension of hit

        Parameters
        ----------
        hits: list
            Coordinates of query subsequences in the query and sequence

        Returns
        ----------
        index[0], index[1], size_max, score_max: tuple
            Returns respectively index of the best hit in the query, in the sequence, max size of the hit, and the best hit score
        '''
        if hits is None:
            hits = self.hits()
            
        assert type(hits) == list and len(hits) > 0, "Parameter hits must be a list and cannot be empty"
        
        list_scores = {}
        list_size = {}
        for i in hits:
            ind_query, ind_seq, size, score = self.extend_hit(i)
            list_scores[i] = score
            list_size[i] = size

        score_max = max(list_scores.values())
        keys = []
        for key, value in list_scores.items():
            if value == score_max:
                keys.append(key)

        final_dic = {}
        for i in keys:
            for key, value in list_size.items():
                if i == key: final_dic[key] = value

        size_max = max(list_size.values())
        for key, value in final_dic.items():
            if value == size_max: 
                index = key
                break

        return index[0], index[1], size_max, score_max

    
class SimpleBlast:
    '''
    Class that records all the sequences to analyse and interacts with the classes SimpleBlastHit and SimpleBlastMatch in order to obtain the respective hits and matches for query and each sequence
    '''
    def __init__(self, list_seqs: (list), w: int):
        '''
        Initiation and recording of all the sequences to analyse
        
        Parameters
        ----------
        list_seqs: list, tuple or str
            List of sequences
            
        w: int
            It determines the size of the subsequences to search in the sequece
        '''
        assert type(list_seqs) in [list, tuple, str] and all(isinstance(s, str) for s in list_seqs if type(list_seqs) in [list, tuple]) and type(w) == int, "list_seqs must be a list, tuple or string. In case is a tuple or a list the members should be strings. W should also be an integer."
        
        if type(list_seqs) == str:
            try:
                f = open(list_seqs, 'r')
                self.seqs = f.read().splitlines()
            except:
                raise Exception("File not found! Must be in the same directory as the working space.")
        else:
            self.seqs = list_seqs
        
        self.w = w
            
    def get_hits(self, seq: str, query: str) -> list:
        '''
        Returns the hits of an specific sequence and inserted query
        
        Parameters
        ----------
        seq: str
            Sequence
            
        query: str
            Query sequence
        
        Returns
        ----------
        hits: list
            Tuple list with the first member of the tuple being the coordinate in the query and the second in the given sequence
        '''
        assert type(seq) == str and type(query) == str, "Sequence and query must be strings"
        seq = seq.upper()
        query = query.upper()
        BlastHits = SimpleBlastHit(seq, query, self.w)
        hits = BlastHits.hits
        if len(hits) == 0: raise Exception("No hits found for the given parameters!")
        return hits
    
    def get_matches(self, seq: str, query: str) -> dict:
        '''
        Returns the matches for all the hits of an specific sequence and inserted query
        
        Parameters
        ----------
        seq: str
            Sequence
        
        query: str
            Query sequence
        
        Returns
        ----------
        matches: dict
            All the hits and respective number of matches
        '''
        assert type(seq) == str and type(query) == str, "Sequence and query must be strings"
        seq = seq.upper()
        query = query.upper()
        BlastHits = SimpleBlastHit(seq, query, self.w)
        matches = BlastHits.get_matches()
        if len(matches) == 0: raise Exception("No matches found for the parameters inserted!")
        return matches
    
    def best_alignment(self, query: str) -> str:
        '''
        Returns the sequence that possesses the best match for the query
        
        Parameters
        ----------
        query: str
            Query sequence
        
        Returns
        ----------
        best_seq: seq
            Sequence with the best match for the given query
        '''
        assert type(query) == str, "Sequence and query must be strings"
        query = query.upper()
        seq_match = {}
        for seq in self.seqs:
            BlastHits = SimpleBlastHit(seq.upper(), query, self.w)
            matches = BlastHits.get_matches()
            if len(matches) > 0:
                seq_match[seq] = max(matches.values())
            else:
                seq_match[seq] = None
        
        best_seq = None
        best_value = 0
        for k, v in seq_match.items():
            if v != None and v > best_value:
                best_value = v
                best_seq = k.upper()
        
        if best_seq == None:
            return 'No sequence matches obtained'
        else:
            return best_seq
    
    
class SimpleBlastHit:
    '''
    Class that records all the hits for a given sequence, query and w. Interacts with SimpleBlast and SimpleBlastMatch
    '''
    def __init__(self, seq: str, query: str, w: int):
        '''
        Initiation and recording of all the hits for the given parameters
        
        Parameters
        ----------
        seq: str
            Sequence
            
        query: str
            Query Sequence
            
        w: int
            It determines the size of the subsequences to search in the sequece
        '''
        assert type(seq) == str and type(query) == str and type(w) == int, "Sequence and query must be strings and W an integer"
        self.seq = seq.upper()
        self.query = query.upper()
        self.w = w
        x = Blast(seq, query, w)
        self.hits = x.hits()
    
    def get_matches(self) -> dict:
        '''
        Construction of a dictionary with the matches for all the hits
        
        Returns
        ----------
        matches: dict
            Hits and corresponding match
        '''
        BlastMatch = SimpleBlastMatch(self.seq, self.query, self.hits, self.w)
        matches = BlastMatch.matches
        return matches
    
    def hits(self) -> list:
        '''
        Function to return the constructed list of hits
        
        Returns
        ----------
        self.hits: list
            Hits for the inputted sequence, query and W in the object initialization
        '''
        return self.hits


class SimpleBlastMatch:
    '''
    Class that records all the matches for a given sequence, query, hits and w. Interacts with SimpleBlast and SimpleBlastHit
    '''
    def __init__(self, seq: str, query: str, hits: list, w: int):
        '''
        Initiation and recording of all the matches for the given parameters
        
        Parameters
        ----------
        seq: str
            Sequence
            
        query: str
            Query sequence
        
        hits: list
            Coordinates of query subsequences in the query and sequence
            
        w: int
            It determines the size of the subsequences to search in the sequece
        '''
        assert type(seq) == str and type(query) == str and type(hits) == list and type(w) == int, "Sequence and query must be strings. Hits must be a list and W an integer."
        self.matches = {}
        for hit in hits:
            x = Blast(seq, query, w)
            self.matches[hit] = x.extend_hit(hit)[3]
            
    def matches(self) -> dict:
        '''
        Function to return the constructed dictionary of hits and matches
        
        Returns
        ----------
        self.matches: dict
            Hits and matches for the inputted sequence, query and W in the object initialization
        '''
        return self.matches
    
