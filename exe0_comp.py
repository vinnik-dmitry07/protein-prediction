# The following method should take a strand as an argument and
# return a complementary strand

def complementary(strand):
    return ""# The following method should take a strand as an argument and
# return a complementary strand

def complementary(strand):
    return strand.translate(
        str.maketrans({
            'A': 'T',
            'T': 'A',
            'G': 'C',
            'C': 'G',
        })
    )
