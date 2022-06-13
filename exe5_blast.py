import numpy as np

from pathlib import Path


"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix (PSSM or substitution
           matrix) parameters. Failure to do so will most likely result in not
           passing the tests.
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}


class BlastDB:

    def __init__(self):
        """
        Initialize the BlastDB class.
        """
        pass

    def add_sequence(self, sequence):
        """
        Add a sequence to the database.

        :param sequence: a protein sequence (string).
        """
        pass

    def get_sequences(self, word):
        """
        Return all sequences in the database containing a given word.

        :param word: a word (string).

        :return: List with sequences.
        """
        return ['NOPE']

    def get_db_stats(self):
        """
        Return some database statistics:
            - Number of sequences in database
            - Number of different words in database
            - Average number of words per sequence (rounded to nearest int)
            - Average number of sequences per word (rounded to nearest int)

        :return: Tuple with four integer numbers corrsponding to the mentioned
                 statistics (in order of listing above).
        """
        return (1, 2, 3, 4)


class Blast:

    def __init__(self, substitution_matrix):
        """
        Initialize the Blast class with the given substitution_matrix.

        :param substitution_matrix: 20x20 amino acid substitution score matrix.
        """
        pass

    def get_words(self, *, sequence=None, pssm=None, T=11):
        """
        Return all words with score >= T for given protein sequence or PSSM.
        Only a sequence or PSSM will be provided, not both at the same time.
        A word may only appear once in the list.

        :param sequence: a protein sequence (string).
        :param pssm: a PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.

        :return: List of unique words.
        """
        return ['AAA', 'YYY']

    def search_one_hit(self, blast_db, *, query=None, pssm=None, T=13, X=5, S=30):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.

        :return: dictionary of target sequences and list of HSP tuples.
        """
        d = dict()
        d['SEQWENCE'] = [(1, 2, 4, 13)]

        return d

    def search_two_hit(self, blast_db, *, query=None, pssm=None, T=11, X=5, S=30, A=40):
        """
        Search a database for target sequences with a given query sequence or
        PSSM. Return a dictionary where the keys are the target sequences for
        which HSPs have been found and the corresponding values are lists of
        tuples. Each tuple is a HSP with the following elements (and order):
            - Start position of HSP in query sequence
            - Start position of HSP in target sequence
            - Length of the HSP
            - Total score of the HSP
        The same HSP may not appear twice in the list (remove duplictes).
        Only a sequence or PSSM will be provided, not both at the same time.

        :param blast_db: BlastDB class object with protein sequences.
        :param query: query protein sequence.
        :param pssm: query PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.
        :param X: drop-off threshold X during extension.
        :param S: score threshold S for the HSP.
        :param A: max distance A between two hits for the two-hit method.

        :return: dictionary of target sequences and list of HSP tuples.
        """
        d = dict()
        d['SEQWENCE'] = [(1, 2, 4, 13)]

        return d
