class Genome:

    def __init__(self, genome):
        """
        Initialize the Genome class with the provided genome sequence.

        :param genome: String with the genome sequence.
        """
        pass

    def get_at_content(self):
        """
        Return the AT content of the genome sequence, i.e. the combined
        fraction of 'A' and 'T' in the entire genome sequence.

        :return: AT content (float, rounded to 6 digits)
        """
        return 0.123456

    def get_codon_dist(self):
        """
        Return the expected codon distribution (fractions) based on the
        distribution (fractions) of the four different nucleotides (ATGC).

        :return: Tree-like structure made out of nested dictionaries. The nodes
                 represent the different nucleotides and the path down the tree
                 forms the corresponding codons. The leafs contain the expected
                 codon frequencies (rounded to 6 digits).
        """
        return {'A': {'A': {'A': 0.123456}}}

    def get_amino_acid_dist(self):
        """
        Return the expected amino acid distribution (fractions) based on the
        expected distribution (fractions) of the different codons.

        :return: Dictionary that contains the expected amino acid distribution.
                 The keys are the 20 different amino acids, the values are the
                 corresponding frequencies (rounded to 6 digits).
        """
        return {'A': 0.123456}
