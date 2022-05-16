import Bio.Seq


def get_orfs(genome, min_num_aa):
    """
    Find and return all ORFs within the provided genome sequence.

    :param genome: String with the genome sequence.
    :param min_num_aa: Minimum number of amino acids in the translated protein
                       sequence.

    :return: List of tuples representing the ORFs (details in worksheet).
    """
    if set(genome) != {'A', 'C', 'T', 'G'}:
        raise TypeError

    sequence = Bio.Seq.Seq(genome)
    res = []
    for s, reverse in [(sequence, False), (sequence.reverse_complement(), True)]:
        for i in range(3):
            aa = s[i:].translate(to_stop=False)
            for t in aa.split('*'):
                if 'M' not in t:
                    continue
                a = t[t.index('M'):]
                if len(a) >= min_num_aa:
                    index = aa.index(a) * 3 + i
                    if reverse:
                        index = len(s) - index - 1
                    res.append((index, (len(a) + 1) * 3, str(a), reverse))
    return res
