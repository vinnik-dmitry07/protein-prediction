from Bio.Data.CodonTable import standard_dna_table


class Genome:

    def __init__(self, genome):
        """
        Initialize the Genome class with the provided genome sequence.

        :param genome: String with the genome sequence.
        """
        self.genome = genome

    def get_at_content(self):
        """
        Return the AT content of the genome sequence, i.e. the combined
        fraction of 'A' and 'T' in the entire genome sequence.

        :return: AT content (float, rounded to 6 digits)
        """
        return (self.genome.count('A') + self.genome.count('T')) / len(self.genome)

    def get_codon_dist(self):
        """
        Return the expected codon distribution (fractions) based on the
        distribution (fractions) of the four different nucleotides (ATGC).

        :return: Tree-like structure made out of nested dictionaries. The nodes
                 represent the different nucleotides and the path down the tree
                 forms the corresponding codons. The leafs contain the expected
                 codon frequencies (rounded to 6 digits).
        """
        fractions = {k: self.genome.count(k) / len(self.genome) for k in 'ATGC'}
        fraction_tree = {
            k1: {
                k2: {
                    k3: round(v1 * v2 * v3, 6) for k3, v3 in fractions.items()
                } for k2, v2 in fractions.items()
            } for k1, v1 in fractions.items()
        }
        return fraction_tree

    def get_amino_acid_dist(self):
        """
        Return the expected amino acid distribution (fractions) based on the
        expected distribution (fractions) of the different codons.

        :return: Dictionary that contains the expected amino acid distribution.
                 The keys are the 20 different amino acids, the values are the
                 corresponding frequencies (rounded to 6 digits).
        """
        codon_dist = self.get_codon_dist()
        acid_dist = {acid: 0 for acid in standard_dna_table.back_table.keys() - {None}}
        stop_codons_sum = sum(codon_dist[c[0]][c[1]][c[2]] for c in standard_dna_table.stop_codons)
        for codon, acid in standard_dna_table.forward_table.items():
            if codon is None:
                continue
            acid_dist[acid] += codon_dist[codon[0]][codon[1]][codon[2]]
        acid_dist = {acid: round(p / (1 - stop_codons_sum), 6) for acid, p in acid_dist.items()}
        return acid_dist
