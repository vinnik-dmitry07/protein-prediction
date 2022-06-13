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
            bg_freq = 0.05
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

        pssm /= pssm.sum(axis=1, keepdims=True)
        pssm = 2 * np.log2(pssm / bg_freq)
        pssm[np.isinf(pssm)] = -20
        return np.rint(pssm).astype(np.int64)

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


# MSA(['SE-AN', 'SE-ES', 'SEVEN', 'SE-AS']).get_pssm(use_sequence_weights=True)
