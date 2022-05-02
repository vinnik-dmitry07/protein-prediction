# This is a standard implementation for the levenshtein distance. 
# You can use it to solve the exercise but you do not have to. 

def levenshteinDistance(s, t):
    """
    Levenshtein implementation by Christopher P. Matthews, taken from  
    https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python

    Args:
        s (str): string1
        t (str): string2

    Returns:
        int: computed levenshtein distance
    """
    if s == t: return 0
    elif len(s) == 0: return len(t)
    elif len(t) == 0: return len(s)
    v0 = [None] * (len(t) + 1)
    v1 = [None] * (len(t) + 1)
    for i in range(len(v0)):
        v0[i] = i
    for i in range(len(s)):
        v1[0] = i + 1
        for j in range(len(t)):
            cost = 0 if s[i] == t[j] else 1
            v1[j + 1] = min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
        for j in range(len(v0)):
            v0[j] = v1[j]
    return v1[len(t)]



# The following function should take max edit distance, a motif and a dna-string and return
# a list of starting positions and length of dna-substrings that match the motif with at
# most max edit distance

def find_similar_motifs(max_edit_distance, motif, dna):
    return []