# Function to return all LCS of substrings `X[0…m-1]`, `Y[0…n-1]`
def LCS(X, Y, m, n, lookup):
    # if the end of either sequence is reached
    if m == 0 or n == 0:
        # create a list with one empty string and return
        return ['']

    # if the last character of `X` and `Y` matches
    if X[m - 1] == Y[n - 1]:

        # ignore the last characters of `X` and `Y` and find all LCS of substring
        # `X[0…m-2]`, `Y[0…n-2]` and store it in a list
        lcs = LCS(X, Y, m - 1, n - 1, lookup)

        # append current character `X[m-1]` or `Y[n-1]`
        # to all LCS of substring `X[0…m-2]` and `Y[0…n-2]`
        for i in range(len(lcs)):
            lcs[i] = lcs[i] + (X[m - 1])

        return lcs

    # we reach here when the last character of `X` and `Y` don't match

    # if a top cell of the current cell has more value than the left cell,
    # then ignore the current character of string `X` and find all LCS of
    # substring `X[0…m-2]`, `Y[0…n-1]`
    if lookup[m - 1][n] > lookup[m][n - 1]:
        return LCS(X, Y, m - 1, n, lookup)

    # if a left cell of the current cell has more value than the top cell,
    # then ignore the current character of string `Y` and find all LCS of
    # substring `X[0…m-1]`, `Y[0…n-2]`
    if lookup[m][n - 1] > lookup[m - 1][n]:
        return LCS(X, Y, m, n - 1, lookup)

    # if the top cell has equal value to the left cell, then consider both characters

    top = LCS(X, Y, m - 1, n, lookup)
    left = LCS(X, Y, m, n - 1, lookup)

    # merge two lists and return
    return top + left


# Function to fill the lookup table by finding the length of LCS
# of substring `X` and `Y`
def LCSLength(X, Y, lookup):
    # fill the lookup table in a bottom-up manner
    for i in range(1, len(X) + 1):
        for j in range(1, len(Y) + 1):
            # if current character of `X` and `Y` matches
            if X[i - 1] == Y[j - 1]:
                lookup[i][j] = lookup[i - 1][j - 1] + 1

            # otherwise, if the current character of `X` and `Y` don't match
            else:
                lookup[i][j] = max(lookup[i - 1][j], lookup[i][j - 1])


# Function to find all LCS of string `X[0…m-1]` and `Y[0…n-1]`
def findLCS(X, Y):
    # lookup[i][j] stores the length of LCS of substring `X[0…i-1]` and `Y[0…j-1]`
    lookup = [[0 for x in range(len(Y) + 1)] for y in range(len(X) + 1)]

    # fill lookup table
    LCSLength(X, Y, lookup)

    # find all the longest common subsequences
    lcs = LCS(X, Y, len(X), len(Y), lookup)

    # since a list can contain duplicates, "convert" it to a set and return
    return set(lcs)


def longest_common_subseq(seq_1, seq_2):
    """
    Return all longest common subsequences between two sequences.
    Problem description: https://en.wikipedia.org/wiki/Longest_common_subsequence_problem

    :param seq_1: String representing the first sequence
    :param seq_2: String representing the second sequence

    :return: List or set of longest common subsequences
    """
    return findLCS(seq_1, seq_2)


def longest_common_subseq_length(seq_1, seq_2):
    """
    Return the length of the longest common subsequences between two sequences.
    Problem description: https://en.wikipedia.org/wiki/Longest_common_subsequence_problem

    :param seq_1: String representing the first sequence
    :param seq_2: String representing the second sequence

    :return: int length of longest common subsequences
    """
    lookup = [[0 for x in range(len(seq_2) + 1)] for y in range(len(seq_1) + 1)]
    LCSLength(seq_1, seq_2, lookup)
    return lookup[-1][-1]


def longest_common_subseq_count(seq_1, seq_2):
    """
    Return the total number of all longest common subsequences between two sequences.
    Problem description: https://en.wikipedia.org/wiki/Longest_common_subsequence_problem

    :param seq_1: String representing the first sequence
    :param seq_2: String representing the second sequence

    :return: int number of longest common subsequences
    """
    return len(findLCS(seq_1, seq_2))


def SCS(X, Y, m, n, lookup):
    # if the end of the first string is reached, create a list
    # containing the second substring and return
    if m == 0:
        return {Y[:n]}

    # if the end of the second string is reached, create a list
    # containing the first substring and return
    elif n == 0:
        return {X[:m]}

    # if the last character of `X` and `Y` is the same, it should occur
    # only one time in SSC
    if X[m - 1] == Y[n - 1]:
        # find all SCS of substring `X[0…m-2]` and `Y[0…n-2]`
        # and append the current character `X[m-1]` or `Y[n-1]` to all SCS of
        # substring `X[0…m-2]` and `Y[0…n-2]`

        scs = SCS(X, Y, m - 1, n - 1, lookup)
        return {s + X[m - 1] for s in scs}

    # we reach here when the last character of `X` and `Y` don't match

    # if a top cell of the current cell has less value than the left cell,
    # then append the current character of string `X` to all SCS of
    # substring `X[0…m-2]`, `Y[0…n-1]`

    if lookup[m - 1][n] < lookup[m][n - 1]:
        scs = SCS(X, Y, m - 1, n, lookup)
        return {s + X[m - 1] for s in scs}

    # if a left cell of the current cell has less value than the top cell,
    # then append the current character of string `Y` to all SCS of
    # substring `X[0…m-1]`, `Y[0…n-2]`

    if lookup[m][n - 1] < lookup[m - 1][n]:
        scs = SCS(X, Y, m, n - 1, lookup)
        return {s + Y[n - 1] for s in scs}

    # if the top cell value is the same as the left cell, then go in both
    # top and left directions

    # append the current character of string `X` to all SCS of
    # substring `X[0…m-2]`, `Y[0…n-1]`
    top = SCS(X, Y, m - 1, n, lookup)

    # append the current character of string `Y` to all SCS of
    # substring `X[0…m-1]`, `Y[0…n-2]`
    left = SCS(X, Y, m, n - 1, lookup)

    return set([s + X[m - 1] for s in top] + [s + Y[n - 1] for s in left])


# Function to fill the lookup table by finding the length of SCS of
# sequences `X[0…m-1]` and `Y[0…n-1]`
def SCSLength(X, Y, m, n, lookup):
    # initialize the first column of the lookup table
    for i in range(m + 1):
        lookup[i][0] = i

    # initialize the first row of the lookup table
    for j in range(n + 1):
        lookup[0][j] = j

    # fill the lookup table in a bottom-up manner
    for i in range(1, m + 1):

        for j in range(1, n + 1):
            # if the current character of `X` and `Y` matches
            if X[i - 1] == Y[j - 1]:
                lookup[i][j] = lookup[i - 1][j - 1] + 1
            # otherwise, if the current character of `X` and `Y` don't match
            else:
                lookup[i][j] = min(lookup[i - 1][j] + 1, lookup[i][j - 1] + 1)


# Function to find all shortest common supersequence of string `X` and `Y`
def findAllSCS(X, Y):
    m = len(X)
    n = len(Y)

    # lookup[i][j] stores the length of SCS of substring `X[0…i-1]` and `Y[0…j-1]`
    lookup = [[0 for x in range(n + 1)] for y in range(m + 1)]

    # fill lookup table
    SCSLength(X, Y, m, n, lookup)

    # find and return all shortest common supersequence
    return SCS(X, Y, m, n, lookup)


def shortest_common_supersequence(seq_1, seq_2):
    """
    Return all shortest common supersequences between two sequences.
    Problem description: https://en.wikipedia.org/wiki/Shortest_common_supersequence_problem

    :param seq_1: String representing the first sequence
    :param seq_2: String representing the second sequence

    :return: List or set of shortest common supersequences
    """
    return findAllSCS(seq_1, seq_2)


def shortest_common_supersequence_length(seq_1, seq_2):
    m = len(seq_1)
    n = len(seq_2)

    # lookup[i][j] stores the length of SCS of substring `X[0…i-1]` and `Y[0…j-1]`
    lookup = [[0 for x in range(n + 1)] for y in range(m + 1)]

    # fill lookup table
    SCSLength(seq_1, seq_2, m, n, lookup)

    return lookup[-1][-1]


def shortest_common_supersequence_count(seq_1, seq_2):
    """
        Return the total number of all shortest common supersequences between two sequences.
        Problem description: https://en.wikipedia.org/wiki/Shortest_common_supersequence_problem

        :param seq_1: String representing the first sequence
        :param seq_2: String representing the second sequence

        :return: int number of shortest common supersequences
    """
    return len(findAllSCS(seq_1, seq_2))
