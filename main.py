from MSA import MSA

GLOBAL = True
MATCH = 5
MISMATCH = -2
INDEL = -4
GAPGAP = 0

# could modify this to use a substitution matrix
def replacementScore(A, B):
    return MATCH if A == B else MISMATCH


def globalAlignTest():
    sequences = [
        "ACTGGTCA",
        "CAGGGTCA",
        "CCAGGGACCA",
    ]
    alignments = MSA(sequences, gap_penalty=INDEL, gap_gap_penalty=GAPGAP,
                     substitution_func=replacementScore)

    for i, a in enumerate(alignments):
        print(f's{i + 1}: {a}')
    print()


def localAlignTest():
    sequences = [
        "ACTGGTCA",
        "CAGGGTCA",
        "CCAGGGACCA",
    ]
    alignments, lpad, rpad = MSA(sequences, gap_penalty=INDEL, gap_gap_penalty=GAPGAP,
                                 substitution_func=replacementScore, global_align=False)

    indicator = " "*(4 + lpad) + "[" + " "*(len(alignments[0])-rpad-lpad-2) + "]"
    print(indicator)
    maxj = len(alignments[0]) - rpad
    for i, a in enumerate(alignments):
        repr = f's{i + 1}: {a}'
        if i != len(alignments)-1:
            next_a = alignments[i+1]
            repr+= "\n    " + (" " * lpad)
            for j in range(lpad, maxj):
                repr += "|" if a[j] == next_a[j] else " "

        print(repr)
    print()


def main():
    globalAlignTest()
    localAlignTest()

if __name__ == '__main__':
    main()
