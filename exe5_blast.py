from collections import defaultdict
from itertools import product

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
WORD_SIZE = 3


class BlastDB:
    def __init__(self):
        """
        Initialize the BlastDB class.
        """
        self.idx_seq = {}
        self.word_num = []
        self.word_seq_idxs = defaultdict(set)

    def add_sequence(self, sequence):
        """
        Add a sequence to the database.

        :param sequence: a protein sequence (string).
        """
        idx = len(self.idx_seq)
        self.idx_seq[idx] = sequence
        words = set(sequence[i:i + WORD_SIZE] for i in range(len(sequence) - WORD_SIZE + 1))
        self.word_num.append(len(words))
        for w in words:
            self.word_seq_idxs[w].add(idx)

    def get_sequences(self, word):
        """
        Return all sequences in the database containing a given word.

        :param word: a word (string).

        :return: List with sequences.
        """
        idxs = self.word_seq_idxs[word]
        return [self.idx_seq[i] for i in idxs]

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
        seq_num = list(map(len, self.word_seq_idxs.values()))
        return (
            len(self.idx_seq),
            len(self.word_seq_idxs),
            round(sum(self.word_num) / len(self.word_num)),
            round(sum(seq_num) / len(seq_num)),
        )


class Blast:
    def __init__(self, substitution_matrix):
        """
        Initialize the Blast class with the given substitution_matrix.

        :param substitution_matrix: 20x20 amino acid substitution score matrix.
        """
        self.sub_mat = substitution_matrix
        self.unique_words = list(''.join(chars) for chars in product(ALPHABET, repeat=3))

    def get_words(self, *, sequence=None, pssm=None, T=11, ret_add=False):
        """
        Return all words with score >= T for given protein sequence or PSSM.
        Only a sequence or PSSM will be provided, not both at the same time.
        A word may only appear once in the list.

        :param sequence: a protein sequence (string).
        :param pssm: a PSSM (Lx20 matrix, where L is length of sequence).
        :param T: score threshold T for the words.

        :return: List of unique words.
        """
        res = set()
        if sequence:
            words = [(i, sequence[i:i + WORD_SIZE]) for i in range(len(sequence) - WORD_SIZE + 1)]
            for i, w1 in words:
                for w2 in self.unique_words:
                    s = sum(self.sub_mat[AA_TO_INT[w1[j]]][AA_TO_INT[w2[j]]] for j in range(WORD_SIZE))
                    if s >= T:
                        res.add((i, s, w2) if ret_add else w2)
        else:
            for i in range(pssm.shape[0] - WORD_SIZE + 1):
                for w in self.unique_words:
                    s = sum(pssm[i + j, AA_TO_INT[w[j]]] for j in range(WORD_SIZE))
                    if s >= T:
                        res.add((i, s, w) if ret_add else w)
        return res

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
        d = defaultdict(set)

        for query_entry, query_score, word in self.get_words(sequence=query, pssm=pssm, T=T, ret_add=True):
            for target in blast_db.get_sequences(word):
                target_entries = [i for i in range(len(target) - WORD_SIZE + 1) if target[i:i + WORD_SIZE] == word]
                for target_entry in target_entries:
                    score = query_score
                    max_score = query_score

                    query_left = query_entry  # query left end
                    target_left = target_entry
                    query_right = query_entry + 2  # query right end

                    query_pos = query_entry + 2 + 1  # query position
                    target_pos = target_entry + 2 + 1

                    while query_pos < len(query) and target_pos < len(target):
                        score += self.sub_mat[AA_TO_INT[query[query_pos]]][AA_TO_INT[target[target_pos]]]

                        if max_score - score >= X:
                            break

                        if score > max_score:
                            max_score = score
                            query_right = query_pos

                        query_pos += 1
                        target_pos += 1

                    score = max_score
                    query_pos = query_entry - 1
                    target_pos = target_entry - 1

                    while query_pos >= 0 and target_pos >= 0:
                        score += self.sub_mat[AA_TO_INT[query[query_pos]]][AA_TO_INT[target[target_pos]]]

                        if max_score - score >= X:
                            break

                        if score > max_score:
                            max_score = score
                            query_left = query_pos
                            target_left = target_pos

                        query_pos -= 1
                        target_pos -= 1

                    if max_score >= S:
                        d[target].add((
                            query_left,
                            target_left,
                            query_right - query_left + 1,
                            max_score
                        ))

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
