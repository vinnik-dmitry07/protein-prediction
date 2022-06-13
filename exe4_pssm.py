import numpy as np

"""
ATTENTION: Use the following dictionaries to get the correct index for each
           amino acid when accessing any type of matrix or array provided as
           parameters. Further, use those indices when generating or returning
           any matrices or arrays. Failure to do so will most likely result in
           not passing the tests.
EXAMPLE: To access the substitution frequency from alanine 'A' to proline 'P'
         in the bg_matrix use bg_matrix[AA_TO_INT['A'], AA_TO_INT['P']].
"""
ALPHABET = 'ACDEFGHIKLMNPQRSTVWY-'
AA_TO_INT = {aa: index for index, aa in enumerate(ALPHABET)}
INT_TO_AA = {index: aa for index, aa in enumerate(ALPHABET)}
GAP_INDEX = AA_TO_INT['-']


def all_equal(lst):
    return lst.count(lst[0]) == len(lst)


blosum62 = {
    'A': {'A': 4, 'C': 0, 'B': -2, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1,
          'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2, 'X': 0, 'Z': -1},
    'C': {'A': 0, 'C': 9, 'B': -3, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1,
          'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2, 'X': -2, 'Z': -3},
    'B': {'A': -2, 'C': -3, 'B': 4, 'E': 1, 'D': 4, 'G': -1, 'F': -3, 'I': -3, 'H': 0, 'K': 0, 'M': -3, 'L': -4, 'N': 3,
          'Q': 0, 'P': -2, 'S': 0, 'R': -1, 'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': 1},
    'E': {'A': -1, 'C': -4, 'B': 1, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0,
          'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 4},
    'D': {'A': -2, 'C': -3, 'B': 4, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4,
          'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3, 'X': -1, 'Z': 1},
    'G': {'A': 0, 'C': -3, 'B': -1, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4,
          'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3, 'X': -1, 'Z': -2},
    'F': {'A': -2, 'C': -2, 'B': -3, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0,
          'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3, 'X': -1, 'Z': -3},
    'I': {'A': -1, 'C': -1, 'B': -3, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2,
          'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1, 'X': -1, 'Z': -3},
    'H': {'A': -2, 'C': -3, 'B': 0, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3,
          'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2, 'X': -1, 'Z': 0},
    'K': {'A': -1, 'C': -3, 'B': 0, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2,
          'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 1},
    'M': {'A': -1, 'C': -1, 'B': -3, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2,
          'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1, 'X': -1, 'Z': -1},
    'L': {'A': -1, 'C': -1, 'B': -4, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4,
          'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1, 'X': -1, 'Z': -3},
    'N': {'A': -2, 'C': -3, 'B': 3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6,
          'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2, 'X': -1, 'Z': 0},
    'Q': {'A': -1, 'C': -3, 'B': 0, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0,
          'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1, 'X': -1, 'Z': 3},
    'P': {'A': -1, 'C': -3, 'B': -2, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3,
          'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3, 'X': -2, 'Z': -1},
    'S': {'A': 1, 'C': -1, 'B': 0, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1,
          'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2, 'X': 0, 'Z': 0},
    'R': {'A': -1, 'C': -3, 'B': -1, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2,
          'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2, 'X': -1, 'Z': 0},
    'T': {'A': 0, 'C': -1, 'B': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1,
          'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2, 'X': 0, 'Z': -1},
    'W': {'A': -3, 'C': -2, 'B': -4, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2,
          'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2, 'X': -2, 'Z': -3},
    'V': {'A': 0, 'C': -1, 'B': -3, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1,
          'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1, 'X': -1, 'Z': -2},
    'Y': {'A': -2, 'C': -2, 'B': -3, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1,
          'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7, 'X': -1, 'Z': -2},
    'X': {'A': 0, 'C': -2, 'B': -1, 'E': -1, 'D': -1, 'G': -1, 'F': -1, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -1,
          'N': -1, 'Q': -1, 'P': -2, 'S': 0, 'R': -1, 'T': 0, 'W': -2, 'V': -1, 'Y': -1, 'X': -1, 'Z': -1},
    'Z': {'A': -1, 'C': -3, 'B': 1, 'E': 4, 'D': 1, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -1, 'L': -3, 'N': 0,
          'Q': 3, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2, 'X': -1, 'Z': 4}
}

blosum_freq = np.array([
    [0.0215, 0.0023, 0.0019, 0.0022, 0.0016, 0.0019, 0.003, 0.0058, 0.0011, 0.0032, 0.0044, 0.0033, 0.0013, 0.0016,
     0.0022, 0.0063, 0.0037, 0.0004, 0.0013, 0.0051],
    [0.0023, 0.0178, 0.002, 0.0016, 0.0004, 0.0025, 0.0027, 0.0017, 0.0012, 0.0012, 0.0024, 0.0062, 0.0008, 0.0009,
     0.001, 0.0023, 0.0018, 0.0003, 0.0009, 0.0016],
    [0.0019, 0.002, 0.0141, 0.0037, 0.0004, 0.0015, 0.0022, 0.0029, 0.0014, 0.001, 0.0014, 0.0024, 0.0005, 0.0008,
     0.0009, 0.0031, 0.0022, 0.0002, 0.0007, 0.0012],
    [0.0022, 0.0016, 0.0037, 0.0213, 0.0004, 0.0016, 0.0049, 0.0025, 0.001, 0.0012, 0.0015, 0.0024, 0.0005, 0.0008,
     0.0012, 0.0028, 0.0019, 0.0002, 0.0006, 0.0013],
    [0.0016, 0.0004, 0.0004, 0.0004, 0.0119, 0.0003, 0.0004, 0.0008, 0.0002, 0.0011, 0.0016, 0.0005, 0.0004, 0.0005,
     0.0004, 0.001, 0.0009, 0.0001, 0.0003, 0.0014],
    [0.0019, 0.0025, 0.0015, 0.0016, 0.0003, 0.0073, 0.0035, 0.0014, 0.001, 0.0009, 0.0016, 0.0031, 0.0007, 0.0005,
     0.0008, 0.0019, 0.0014, 0.0002, 0.0007, 0.0012],
    [0.003, 0.0027, 0.0022, 0.0049, 0.0004, 0.0035, 0.0161, 0.0019, 0.0014, 0.0012, 0.002, 0.0041, 0.0007, 0.0009,
     0.0014, 0.003, 0.002, 0.0003, 0.0009, 0.0017],
    [0.0058, 0.0017, 0.0029, 0.0025, 0.0008, 0.0014, 0.0019, 0.0378, 0.001, 0.0014, 0.0021, 0.0025, 0.0007, 0.0012,
     0.0014, 0.0038, 0.0022, 0.0004, 0.0008, 0.0018],
    [0.0011, 0.0012, 0.0014, 0.001, 0.0002, 0.001, 0.0014, 0.001, 0.0093, 0.0006, 0.001, 0.0012, 0.0004, 0.0008, 0.0005,
     0.0011, 0.0007, 0.0002, 0.0015, 0.0006],
    [0.0032, 0.0012, 0.001, 0.0012, 0.0011, 0.0009, 0.0012, 0.0014, 0.0006, 0.0184, 0.0114, 0.0016, 0.0025, 0.003,
     0.001, 0.0017, 0.0027, 0.0004, 0.0014, 0.012],
    [0.0044, 0.0024, 0.0014, 0.0015, 0.0016, 0.0016, 0.002, 0.0021, 0.001, 0.0114, 0.0371, 0.0025, 0.0049, 0.0054,
     0.0014, 0.0024, 0.0033, 0.0007, 0.0022, 0.0095],
    [0.0033, 0.0062, 0.0024, 0.0024, 0.0005, 0.0031, 0.0041, 0.0025, 0.0012, 0.0016, 0.0025, 0.0161, 0.0009, 0.0009,
     0.0016, 0.0031, 0.0023, 0.0003, 0.001, 0.0019],
    [0.0013, 0.0008, 0.0005, 0.0005, 0.0004, 0.0007, 0.0007, 0.0007, 0.0004, 0.0025, 0.0049, 0.0009, 0.004, 0.0012,
     0.0004, 0.0009, 0.001, 0.0002, 0.0006, 0.0023],
    [0.0016, 0.0009, 0.0008, 0.0008, 0.0005, 0.0005, 0.0009, 0.0012, 0.0008, 0.003, 0.0054, 0.0009, 0.0012, 0.0183,
     0.0005, 0.0012, 0.0012, 0.0008, 0.0042, 0.0026],
    [0.0022, 0.001, 0.0009, 0.0012, 0.0004, 0.0008, 0.0014, 0.0014, 0.0005, 0.001, 0.0014, 0.0016, 0.0004, 0.0005,
     0.0191, 0.0017, 0.0014, 0.0001, 0.0005, 0.0012],
    [0.0063, 0.0023, 0.0031, 0.0028, 0.001, 0.0019, 0.003, 0.0038, 0.0011, 0.0017, 0.0024, 0.0031, 0.0009, 0.0012,
     0.0017, 0.0126, 0.0047, 0.0003, 0.001, 0.0024],
    [0.0037, 0.0018, 0.0022, 0.0019, 0.0009, 0.0014, 0.002, 0.0022, 0.0007, 0.0027, 0.0033, 0.0023, 0.001, 0.0012,
     0.0014, 0.0047, 0.0125, 0.0003, 0.0009, 0.0036],
    [0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0002, 0.0003, 0.0004, 0.0002, 0.0004, 0.0007, 0.0003, 0.0002, 0.0008,
     0.0001, 0.0003, 0.0003, 0.0065, 0.0009, 0.0004],
    [0.0013, 0.0009, 0.0007, 0.0006, 0.0003, 0.0007, 0.0009, 0.0008, 0.0015, 0.0014, 0.0022, 0.001, 0.0006, 0.0042,
     0.0005, 0.001, 0.0009, 0.0009, 0.0102, 0.0015],
    [0.0051, 0.0016, 0.0012, 0.0013, 0.0014, 0.0012, 0.0017, 0.0018, 0.0006, 0.012, 0.0095, 0.0019, 0.0023, 0.0026,
     0.0012, 0.0024, 0.0036, 0.0004, 0.0015, 0.0196],
])

mapping = [0, 4, 3, 6, 13, 7, 8, 9, 11, 10, 12, 2, 14, 5, 1, 15, 16, 19, 17, 18]
blosum_freq = blosum_freq[mapping][:, mapping]


class MSA:
    def __init__(self, sequences):
        if not (
                len(sequences) > 0
                and all_equal([len(s) for s in sequences])
                and all(all(si in ALPHABET for si in s) for s in sequences)
        ):
            raise TypeError

        self.sequences = sequences
        self.sequences_matrix = np.array(sequences, dtype=bytes).view('S1').reshape((len(sequences), -1))

        self.r = np.apply_along_axis(lambda col: np.unique(col).size, axis=0, arr=self.sequences_matrix)

        self.seqs_elems_entries = self.get_seqs_elems_entries()

    def get_seqs_elems_entries(self, weights=None):
        a = np.array(list(ALPHABET), dtype='S1')
        seqs_elems_entries = np.empty((self.sequences_matrix.shape[1], a.size))

        for i in range(seqs_elems_entries.shape[0]):
            b = self.sequences_matrix[:, i]
            idx_sorted = a.argsort()
            idx = idx_sorted[np.searchsorted(a, b, sorter=idx_sorted)]
            mask = a[idx] == b
            seqs_elems_entries[i] = np.bincount(idx[mask], minlength=a.size, weights=weights)

        return seqs_elems_entries

    def get_pssm(self, *, bg_matrix=None, beta=10, use_sequence_weights=False,
                 redistribute_gaps=False, add_pseudocounts=False):
        """
        Return a PSSM for the underlying MSA. Use the appropriate refinements 
        according to the parameters. If no bg_matrix is specified, use uniform 
        background and pair frequencies.
        Every row in the resulting PSSM corresponds to a non-gap position in 
        the primary sequence of the MSA (i.e. the first one).
        Every column in the PSSM corresponds to one of the 20 amino acids.
        Values that would be -inf must be replaced by -20 in the final PSSM.
        Before casting to dtype=numpy.int64, round all values to the nearest
        integer (do not just FLOOR all values).

        :param bg_matrix: Amino acid pair frequencies as numpy array (20, 20).
                          Access the matrix using the indices from AA_TO_INT.
        :param beta: Beta value (float) used to weight the pseudocounts 
                     against the observed amino acids in the MSA.
        :param use_sequence_weights: Calculate and apply sequence weights.
        :param redistribute_gaps: Redistribute the gaps according to the 
                                  background frequencies.
        :param add_pseudocounts: Calculate and add pseudocounts according 
                                 to the background frequencies.

        :return: PSSM as numpy array of shape (L x 20, dtype=numpy.int64).
                 L = ungapped length of the primary sequence.
        """
        if bg_matrix is None:
            bg_freq = 0.05 * np.ones(20)
        else:
            bg_freq = np.sum(bg_matrix, axis=0)

        non_gap_mask = self.sequences_matrix[0] != b'-'

        if use_sequence_weights:
            f_matrix = self.get_seqs_elems_entries(self.get_sequence_weights())
        else:
            f_matrix = self.seqs_elems_entries

        pssm = f_matrix[non_gap_mask, :-1]

        if redistribute_gaps:
            pssm += f_matrix[non_gap_mask, -1:] * bg_freq

        if add_pseudocounts:
            ps = np.empty_like(pssm)
            for i in range(ps.shape[0]):
                for a in range(ps.shape[1]):
                    ps[i, a] = np.sum(blosum_freq[a] * pssm[i] / bg_freq)
                    alpha = self.get_number_of_observations() - 1
            pssm = (alpha * pssm + beta * ps) / (alpha + beta)

        pssm /= pssm.sum(axis=1, keepdims=True)
        pssm = 2 * np.log2(pssm / bg_freq)
        pssm[np.isinf(pssm)] = -20
        res = np.rint(pssm).astype(np.int64)
        return res

    def get_size(self):
        """
        Return the number of sequences in the MSA and the MSA length, i.e.
        the number of columns in the MSA. This includes gaps.

        :return: Tuple of two integers. First element is the number of
                 sequences in the MSA, second element is the MSA length.
        """
        return self.sequences_matrix.shape

    def get_primary_sequence(self):
        """
        Return the primary sequence of the MSA. In this exercise, the primary
        sequence is always the first sequence of the MSA. The returned 
        sequence must NOT include gap characters.

        :return: String containing the ungapped primary sequence.
        """
        return self.sequences[0].replace('-', '')

    def get_sequence_weights(self):
        """
        Return the calculated sequence weights for all sequences in the MSA.
        The order of weights in the array must be equal to the order of the
        sequences in the MSA.

        :return: Numpy array (dtype=numpy.float64) containing the weights for
                 all sequences in the MSA.
        """

        a = np.array(list(ALPHABET), dtype='S1')
        matrix_alphabet_indexes_sorted = np.searchsorted(np.sort(a), self.sequences_matrix)
        matrix_alphabet_indexes = np.argsort(a)[matrix_alphabet_indexes_sorted]

        s = np.take_along_axis(self.seqs_elems_entries, matrix_alphabet_indexes.T, axis=1)
        w = 1 / (self.r[:, np.newaxis] * s)
        return w[self.r != 1].sum(axis=0)

    def get_number_of_observations(self):
        """
        Return the estimated number of independent observations in the MSA.

        :return: Estimate of independent observation (dtype=numpy.float64).
        """
        return np.mean(self.r)
