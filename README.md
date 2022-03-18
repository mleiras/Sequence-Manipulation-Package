# Sequence Manipulation Package
The Sequence Manipulation package constitutes a set of classes and functions for biological sequences manipulation available for Python.

## What can be found in this package?
The Sequence Manipulation package has some functionalities such as:
- Sequence classes that deals with sequences, namely DNA, RNA and amino acid sequences;
- Functions to perform common operations on sequences such as sequence features (content), transcription, open reading frames and translation;
- Code for dealing with alignments, including a way to deal with local, global and multiple alignment and also search for good local alignment (BLAST);
- Code for dealing with probabilistic profiles;
- Documentation and help with using the classes, including this file.

## What can be done with the Sequence Manipulation package?
This package contains the following classes:
- Class Sequence;
- Class Alignment;
- Class Blast;
- Class Motifs;
- Class Progressive Alignment.
The features of each class will be addressed in the next sections, along with usage examples.

### Class Sequence
The class Sequence provides a set of sequence features, such as the type of sequence (DNA, RNA or amino acid) and its elements. Besides, it allows the transformation of DNA to RNA sequences, determination of the complement inverse of a sequence and complementary DNA, open reading frames, codons and respective amino acids and proteins generated.
**The parameters are:**
- seq: sequence.

#### Usage Example
An example of the get_all_prots function (which returns all the proteins from the sequence) of the class Sequence is:
```
seq_dna = Sequence('ATGTTAGAGTTATTAAAAAGTCTGGTATTCGCCGTAATCATGGTACCTGTCGTGATGGCCATCATCCTGGGTCTGATTTACGGTCTTGGTGAAGTATTCAACATCTTTTCTGGTGTTGGTAAAAAAGACCAGCCCGGACAAAATCATTGAATTTAATTACAAGTCTTCAGAATGCCAGAGATATACAGGATCTAACCA')
print(seq_dna.get_all_prots())
```
An output example of the function is:

```['MLELLKSLVFAVIMVPVVMAIILGLIYGLGEVFNIFSGVGKKDQPGQNH_','MPEIYRI_']```

### Class Alignment
The class Alignment constructs local and global alignments using a blosum62 matrix. Its functions and parameters (into brackets) are:
- seq1, seq2: sequences of which the matrix will be constructed;
- file: file location of the matrix with the characters values;
- g: spacing value;
- l: verifys the use of local alignment;
- a, b: characters intended to verify the value that corresponds to the individual caracteres inputted by the user;
- ties: takes False argument as default and indicates if we should assume the ties in the matrix or not;
- i: position that will correspond to the position to be incremented on the lists
- x, y: coordinates to begin the recursive function in the matrix and correctly implement in the lists the alignments.

### Class Blast
This section is represented by 4 classes that may interact with each other. The classes and respective parameters (in order of appearance) are:

1. The **class Blast** aims to perform a BLAST between a sequence query, the given sequence and the respective size of the query hits that we want to obtain in the sequence. For that, the parameters are:
- seq : sequence;
- query : sequence query;
- w: size of the subsequences;
- q: query subsequence;
- query: query sequence serving as a template to find the indexes;
- qm: query subsequences with the corresponding indexes in the query;
- seq: sequence serving as a template to search the indexes;
- index : specific hit to be extended in each direction;
- hits: coordinates of query subsequences in the query and sequence.

2. The **class SimpleBlast** records all the sequences to analyze and returns the respective hits and matches for the query and each sequence. The parameters are:
● list_seqs: list of sequences;
● w: size of the subsequences;
● seq: sequence;
● query: query sequence.

3. The **class SimpleBlastHit** records all the hits for a given sequence, query and subsequence (w). The parameters are:
● seq: sequence;
● query: query sequence;
● w: size of the subsequences.

4. The **class SimpleBlastMatch** records all the matches for a given sequence, query, hits and subsequences. The parameters are:
● seq: sequence;
● query: query sequence;
● hits: coordinates of query subsequences;
● w: size of the subsequences.

### Class Motifs 
The class Motifs creates and displays probabilistic profiles through Position Weight Matrices (PWM) and/or Position-Specific Scoring Matrices (PSSM). This class has the following parameters:
- list_seq: list of sequences;
- kwargs;
- seq: sequence;
- profile: dictionary of the probabilistic profile (PWM or PSSM).

### Class Progressive Alignments
The class Progressive Alignments performs multiple alignments giving a list of sequences. The parameters are:
- list_seqs: list of the multiple alignments (DNA or Amino Acid);
- kwargs;
- s1: sequence of DNA or Amino Acids;
- s2: sequence of DNA or Amino Acids.
