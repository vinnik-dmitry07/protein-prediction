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
        acids = ['K', 'N', 'T', 'R', 'S', 'I', 'M', 'Q', 'H', 'P', 'L', 'E', 'D', 'A', 'G', 'V', 'Y', 'C', 'W', 'F']
        stop_codons = ['TAA', 'TAG', 'TGA']
        forward_table = {
            'TTT': 'F',
             'TTC': 'F',
             'TTA': 'L',
             'TTG': 'L',
             'TCT': 'S',
             'TCC': 'S',
             'TCA': 'S',
             'TCG': 'S',
             'TAT': 'Y',
             'TAC': 'Y',
             'TGT': 'C',
             'TGC': 'C',
             'TGG': 'W',
             'CTT': 'L',
             'CTC': 'L',
             'CTA': 'L',
             'CTG': 'L',
             'CCT': 'P',
             'CCC': 'P',
             'CCA': 'P',
             'CCG': 'P',
             'CAT': 'H',
             'CAC': 'H',
             'CAA': 'Q',
             'CAG': 'Q',
             'CGT': 'R',
             'CGC': 'R',
             'CGA': 'R',
             'CGG': 'R',
             'ATT': 'I',
             'ATC': 'I',
             'ATA': 'I',
             'ATG': 'M',
             'ACT': 'T',
             'ACC': 'T',
             'ACA': 'T',
             'ACG': 'T',
             'AAT': 'N',
             'AAC': 'N',
             'AAA': 'K',
             'AAG': 'K',
             'AGT': 'S',
             'AGC': 'S',
             'AGA': 'R',
             'AGG': 'R',
             'GTT': 'V',
             'GTC': 'V',
             'GTA': 'V',
             'GTG': 'V',
             'GCT': 'A',
             'GCC': 'A',
             'GCA': 'A',
             'GCG': 'A',
             'GAT': 'D',
             'GAC': 'D',
             'GAA': 'E',
             'GAG': 'E',
             'GGT': 'G',
             'GGC': 'G',
             'GGA': 'G',
             'GGG': 'G'
        }

        codon_dist = self.get_codon_dist()
        acid_dist = {a: 0 for a in acids}
        stop_codons_sum = sum(codon_dist[c[0]][c[1]][c[2]] for c in stop_codons)
        for codon, acid in forward_table.items():
            acid_dist[acid] += codon_dist[codon[0]][codon[1]][codon[2]]
        acid_dist = {acid: round(p / (1 - stop_codons_sum), 6) for acid, p in acid_dist.items()}
        return acid_dist
