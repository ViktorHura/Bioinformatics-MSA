from math import inf
from itertools import product, combinations

MATCH = 5
MISMATCH = -2
INDEL = -4
GAPGAP = 0

# could modify this to use a substitution matrix
def replacementScore(A, B):
    return MATCH if A == B else MISMATCH


class NMatrix:
    # dimension_length in the form of [..., #pages, #rows, #cols]
    def __init__(self, dimension_lengths):
        self.dims = dimension_lengths

        self.total_size = 1
        for l in dimension_lengths:
            self.total_size *= l

        # each cell in the form of (score, parent coordinate mask)
        self.data = [(None, [])] * self.total_size

    # convert n-tuple to index in our 1D array
    def __tupleToIndex__(self, key):
        index = 0
        # sum of products, in 3D: index(p,r,c) = (p * #rows * #cols) + (r * #cols) + c
        for i, coord in enumerate(key):
            # last term is not a product
            if i + 1 == len(self.dims):
                index += coord
                continue
            product = coord
            # product of all subsequent dimension lengths
            for j in range(i + 1, len(self.dims)):
                product *= self.dims[j]
            index += product
        return index

    def __getitem__(self, key):
        for k in key:
            if k < 0:
                return None
        index = self.__tupleToIndex__(key)
        if index < 0 or index >= len(self.data):
            return None

        return self.data[index]

    def __setitem__(self, key, value):
        index = self.__tupleToIndex__(key)
        self.data[index] = value


# relative coordinates of each to consider neighbour in n dimensions
# will yield 3 neighbours in 2D, 7 in 3D, 2^n - 1 neighbours in N-D
def neighbourMask(dim_count):
    masks = list(product([0, -1], repeat=dim_count))
    masks.remove(tuple([0]*dim_count))
    return [list(m) for m in masks]


# increments the current matrix position in order
def nextMatrixPosition(current_position, dimensions):
    new_position = current_position
    i = len(current_position) - 1   # start with last value in tuple
    while(i != -1):                 # increment from back to front and carry over if needed
        new_position[i] += 1
        if new_position[i] == dimensions[i]:
            new_position[i] = 0
            i -= 1
        else:
            break
    return new_position


# element wise sum of two positions
def vectorSum(v1, v2):
    return [v + v2[i] for i, v in enumerate(v1)]


# returns score derived from a neighbour cell
# given the current position, sequences, neighbour score and mask for neighbour position
def parseNeighbour(position, sequences, n_score, n_mask):
    score = n_score

    pairwise_combos = list(combinations(range(len(n_mask)), 2))
    for combo in pairwise_combos:
        aI, bI = combo
        im = n_mask[aI]
        jm = n_mask[bI]

        if (im, jm) == (0,0):
            score += GAPGAP
        elif (im, jm) == (0, -1) or (im, jm) == (-1, 0):
            score += INDEL
        elif (im, jm) == (-1, -1):
            seqA = sequences[aI]
            seqB = sequences[bI]
            A = seqA[position[aI]]
            B = seqB[position[bI]]
            score += replacementScore(A,B)
    return score


# get global alignment path in a filled NMatrix
# path in the form of [ (parent mask, cell position) ,...]
def globalAlignmentPath(matrix):
    path = []
    cur_pos = [d-1 for d in matrix.dims]        # start from the end
    origin = [0] * len(matrix.dims)             # stop when back at origin
    while cur_pos != origin:
        cell = matrix[cur_pos]
        path.append((cell[1], cur_pos))         # (parent mask, cell position)
        cur_pos = vectorSum(cur_pos, cell[1])   # next position = parent position (via mask)
    path.reverse()
    return path


def parseGlobalAlignments(sequences, path):
    alignments = ["" for s in sequences]    # store alignments as list of strings
    maxlen = 0                              # longest alignment length, to pad out the other ones

    for node in path:
        delta = node[0]
        for pos, coord in enumerate(delta):
            sequence = sequences[pos]
            letter_index = node[1][pos]
            letter_to_consume = sequence[letter_index]
            alignments[pos] += "." if coord == 0 else letter_to_consume

    for i, a in enumerate(alignments):
        # take sequence and replace each letter that is consumed already in alignment
        # this gives us the leftover letters
        residue = sequences[i][1:]
        for c in a:
            residue = residue.replace(c, '', 1)

        # append leftover to the end of the alignment
        alignments[i] += residue

        # keep track of longest alignment length
        length = len(alignments[i])
        if length > maxlen:
            maxlen = length

    # pad ends of alignments if they are not maxlen length
    for i, a in enumerate(alignments):
        dif = maxlen - len(a)
        for j in range(dif):
            alignments[i] += "."

    return alignments


def MSA(sequences, global_align=True):
    seqs = ["." + s for s in sequences]
    dimensions = [len(s) for s in seqs]

    matrix = NMatrix(dimensions)
    best_cell = (-inf, None)

    # origin = 0
    cur_pos = [0] * len(dimensions)
    matrix[cur_pos] = (0, [])
    # next pos
    cur_pos[-1] = 1

    final_pos = [d-1 for d in dimensions]
    nMasks = neighbourMask(len(dimensions))

    print(f'Filling matrix of dimensions: {dimensions}\n')

    # Fill each cell in matrix in order
    while(True):
        # neighbouring cells to consider and a place to store their scores
        neighbours = [vectorSum(cur_pos, mask) for mask in nMasks]
        scores = [-inf]*len(neighbours)

        # parse scores from each neighbour (if exists)
        for i, nPos in enumerate(neighbours):
            nVal = matrix[nPos]
            if nVal is None:
                continue
            scores[i] = parseNeighbour(cur_pos, seqs, nVal[0], nMasks[i])

        # get max score and the neighbour mask that got us there
        maxScore = max(scores)
        maxI = scores.index(maxScore)
        delta = nMasks[maxI]

        # current cel = (score, parent coordinate mask)
        matrix[cur_pos] = (maxScore, delta)

        # keep track of best cell
        if maxScore > best_cell[0]:
            best_cell = (maxScore, cur_pos)

        # stop if matrix filled, else next cel position
        if cur_pos == final_pos:
            break
        cur_pos = nextMatrixPosition(cur_pos, dimensions)

    if global_align:
        print(f'Final alignment score: {matrix[final_pos][0]}\n')
        path = globalAlignmentPath(matrix)

        return parseGlobalAlignments(seqs, path)


def main():
    sequences = [
        "CAG",
        "ACT",
    ]
    alignments = MSA(sequences)

    for i, a in enumerate(alignments):
        print(f's{i+1}: {a}')

if __name__ == '__main__':
    main()
