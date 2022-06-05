def complementary(strand):
    return strand.translate(
        str.maketrans({
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
        })
    )

gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }


def get_orfs(sequence, min_num_aa):
    """
    Find and return all ORFs within the provided genome sequence.

    :param genome: String with the genome sequence.
    :param min_num_aa: Minimum number of amino acids in the translated protein
                       sequence.

    :return: List of tuples representing the ORFs (details in worksheet).
    """
    if set(sequence) != {'A', 'C', 'T', 'G'}:
        raise TypeError

    res = []
    for s, reverse in [(sequence, False), (complementary(sequence)[::-1], True)]:
        for i in range(3):
            s1 = s[i:] + s[:i] + s[i:]
            aa = ''.join(gencode[s1[j:j+3]] for j in range(0, 3 * (len(s1) // 3), 3))
            for t in aa.split('*')[:-1]:
                if 'M' not in t:
                    continue
                a = t[t.index('M'):]
                if len(a) >= min_num_aa:
                    index = aa.index(a) * 3 + i
                    index = index % len(s)
                    if reverse:
                        index = len(s) - index - 1
                    res.append((index, (len(a) + 1) * 3, str(a), reverse))
    res = list(set(res))
    return res
