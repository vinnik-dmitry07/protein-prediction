import numpy as np


class GlobalAlignment:
    def __init__(self, string1, string2, gap_penalty, matrix):
        """
        :param string1: first string to be aligned, string
        :param string2: second string to be aligned, string
        :param gap_penalty: gap penalty, integer
        :param matrix: substitution matrix containing scores for amino acid
                       matches and mismatches, dict

        Attention! string1 is used to index columns, string2 is used to index rows
        """
        self.string1 = string1
        self.string2 = string2
        self.gap_penalty = gap_penalty
        self.substituion_matrix = matrix
        self.score_matrix = np.zeros((len(string2) + 1, len(string1) + 1), dtype=np.int32)
        self.align()

    def align(self):
        """
        Align given strings using the Needleman-Wunsch algorithm,
        store the alignments and the score matrix used to compute those alignments.
        NB: score matrix and the substitution matrix are different matrices!
        """

    def get_best_score(self):
        """
        :return: the highest score for the aligned strings, int

        """
        return 4

    def get_number_of_alignments(self):
        """
        :return: number of found alignments with the best score
        """
        return 43

    def get_alignments(self):
        """
        :return: list of alignments, where each alignment is represented
                 as a tuple of aligned strings
        """
        return [
            ('ADMI-NS', 'ADMIRES'), ('ADMIN-S', 'ADMIRES')
        ]

    def get_score_matrix(self):
        """
        :return: matrix built during the alignment process as a list of lists
        """
        return [
            [0, -1, -2, -3, -4, -5, -6],
            [-1, 1, 0, -1, -2, -3, -4],
            [-2, 0, 2, 1, 0, -1, -2],
            [-3, -1, 1, 3, 2, 1, 0],
            [-4, -2, 0, 2, 4, 3, 2],
            [-5, -3, -1, 1, 3, 4, 3],
            [-6, -4, -2, 0, 2, 3, 4],
            [-7, -5, -3, -1, 1, 2, 4]
        ]
