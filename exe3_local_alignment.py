from collections import defaultdict
from copy import deepcopy

import numpy as np

DIAG = 2
LEFT = 3
UP = 4

CONTAINS_DIAG = [DIAG, DIAG + LEFT, DIAG + UP, DIAG + LEFT + UP]
CONTAINS_LEFT = [LEFT, LEFT + DIAG, LEFT + UP, LEFT + DIAG + UP]
CONTAINS_UP = [UP, UP + DIAG, UP + LEFT, UP + DIAG + LEFT]


class LocalAlignment:
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
        self.subst_mat = matrix
        self.score_mat = None
        self.alignments = []
        self.align()

    def align(self):
        """
        Align given strings using the Smith-Waterman algorithm.
        NB: score matrix and the substitution matrix are different matrices!
        """
        x = self.string1
        y = self.string2

        nx = len(x)
        ny = len(y)
        # Optimal score at each possible pair of characters.
        score_mat = np.zeros((nx + 1, ny + 1), dtype=int)
        # Pointers to trace through an optimal aligment.
        point_mat = np.zeros((nx + 1, ny + 1), dtype=int)
        point_mat[:, 0] = UP
        point_mat[0, :] = LEFT
        # Temporary scores.
        t = np.zeros(4)
        for i in range(nx):
            for j in range(ny):
                t[0] = score_mat[i, j] + self.subst_mat[x[i]][y[j]]
                t[1] = score_mat[i, j + 1] + self.gap_penalty
                t[2] = score_mat[i + 1, j] + self.gap_penalty
                t[3] = 0
                max_t = np.max(t)
                score_mat[i + 1, j + 1] = max_t
                if t[0] == max_t:
                    point_mat[i + 1, j + 1] += DIAG
                if t[1] == max_t:
                    point_mat[i + 1, j + 1] += UP
                if t[2] == max_t:
                    point_mat[i + 1, j + 1] += LEFT

        def rec(i, j, rx, ry, idx, ridx):
            if i <= 0 and j <= 0:
                return rx, ry

            if score_mat[i, j] == 0:
                ridx.append(idx)
                return rx, ry

            idx1 = idx + (point_mat[i, j],)

            rx.setdefault(idx, {''})
            ry.setdefault(idx, {''})
            rx.setdefault(idx1, set())
            ry.setdefault(idx1, set())

            rx1_ry1_i = []
            if point_mat[i, j] in CONTAINS_DIAG:
                rx_diag = deepcopy(rx)
                ry_diag = deepcopy(ry)
                rx_diag[idx1] |= {x[i - 1] + rx_idx_i for rx_idx_i in rx_diag[idx]}
                ry_diag[idx1] |= {y[j - 1] + ry_idx_i for ry_idx_i in ry_diag[idx]}
                rx1_ry1_i.append(rec(i - 1, j - 1, rx_diag, ry_diag, idx1, ridx))
            if point_mat[i, j] in CONTAINS_LEFT:
                rx_left = deepcopy(rx)
                ry_left = deepcopy(ry)
                rx_left[idx1] |= {'-' + rx_idx_i for rx_idx_i in rx_left[idx]}
                ry_left[idx1] |= {y[j - 1] + ry_idx_i for ry_idx_i in ry_left[idx]}
                rx1_ry1_i.append(rec(i, j - 1, rx_left, ry_left, idx1, ridx))
            if point_mat[i, j] in CONTAINS_UP:
                rx_up = deepcopy(rx)
                ry_up = deepcopy(ry)
                rx_up[idx1] |= {x[i - 1] + rx_idx_i for rx_idx_i in rx_up[idx]}
                ry_up[idx1] |= {'-' + ry_idx_i for ry_idx_i in ry_up[idx]}
                rx1_ry1_i.append(rec(i - 1, j, rx_up, ry_up, idx1, ridx))

            rx1 = defaultdict(set)
            ry1 = defaultdict(set)
            for rx1_i, ry1_i in rx1_ry1_i:
                for key, value in rx1_i.items():
                    rx1[key] |= value
                for key, value in ry1_i.items():
                    ry1[key] |= value

            return rx1, ry1

        ridx = []
        start_i, start_j = np.unravel_index(score_mat.argmax(), score_mat.shape)
        rx, ry = rec(i=start_i, j=start_j, rx={}, ry={}, idx=(), ridx=ridx)

        for idx in ridx:
            self.alignments.append((''.join(rx[idx]), ''.join(ry[idx])))

        self.score_mat = score_mat

    def has_alignment(self):
        """
        :return: True if a local alignment has been found, False otherwise
        """
        return bool(len(self.alignments))

    def get_alignment(self):
        """
        :return: alignment represented as a tuple of aligned strings
        """
        return self.alignments[0] if self.alignments else ('', '')

    def is_residue_aligned(self, string_number, residue_index):
        """
        :param string_number: number of the string (1 for string1, 2 for string2) to check
        :param residue_index: index of the residue to check
        :return: True if the residue with a given index in a given string has been aligned
                 False otherwise
        """
        start_i, start_j = np.unravel_index(self.score_mat.argmax(), self.score_mat.shape)
        if string_number == 1:
            return start_i + 1 - len(self.alignments[0][0]) <= residue_index <= start_i - 1
        elif string_number == 2:
            return start_j + 1 - len(self.alignments[0][0]) <= residue_index <= start_j - 1
        else:
            raise ValueError(string_number)
