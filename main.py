from math import inf
from itertools import product, combinations

MATCH = 5
MISMATCH = -2
INDEL = -4
MULTI_GAP_THRESHOLD = 2
MULTI_GAP = 0

class NMatrix:
    def __init__(self, dimension_lengths):
        self.dims = dimension_lengths

        self.total_size = 1
        for l in dimension_lengths:
            self.total_size *= l

        self.data = [(None, [])] * self.total_size

    def __keyToIndex__(self, key):
        index = 0
        for i, coord in enumerate(key):
            if i + 1 == len(self.dims):
                index += coord
                continue
            product = coord
            for j in range(i + 1, len(self.dims)):
                product *= self.dims[j]
            index += product
        return index

    def __getitem__(self, key):
        index = self.__keyToIndex__(key)
        if index < 0 or index >= len(self.data):
            return None
        for k in key:
            if k < 0:
                return None

        return self.data[index]

    def __setitem__(self, key, value):
        index = self.__keyToIndex__(key)
        self.data[index] = value


    def __repr__(self):
        if len(self.dims) == 3:
            result = ""
            for p in range(self.dims[0]):

                for r in range(self.dims[1]):
                    rowStr = ""
                    for c in range(self.dims[2]):
                        rowStr += str(self.__getitem__((p, r, c))[0]) + ", "
                    result += rowStr + '\n'

                result += '\n'
            return result
        else:
            return str(self.data)


def neighbourMask(dim_count):
    masks = list(product([0, -1], repeat=dim_count))
    masks.remove(tuple([0]*dim_count))
    return [list(m) for m in masks]


def nextMatrixPosition(current_position, dimensions):
    new_position = current_position
    i = len(current_position) - 1
    while(i != -1):
        new_position[i] += 1
        if new_position[i] == dimensions[i]:
            new_position[i] = 0
            i -= 1
        else:
            break
    return new_position


def vectorSum(v1, v2):
    return [v + v2[i] for i, v in enumerate(v1)]


def MSA(sequences):
    seqs = ["." + s for s in sequences]
    dimensions = [len(s) for s in seqs]

    matrix = NMatrix(dimensions)
    best_cell = (-inf, None)

    cur_pos = [0] * len(dimensions)
    matrix[cur_pos] = (0, [])
    #next pos
    cur_pos[-1] = 1

    final_pos = [d-1 for d in dimensions]
    #final_pos = [0,0,1]

    nMasks = neighbourMask(len(dimensions))

    print(f'Dimensions: {dimensions}\n')

    #TODO CLEAN CODE
    #TODO BETTER SUM OF PAIRS

    while(True):
        neighbours = [vectorSum(cur_pos, mask) for mask in nMasks]
        scores = [-inf]*len(neighbours)

        for i, nPos in enumerate(neighbours):
            nVal = matrix[nPos]
            if nVal is None:
                continue

            scores[i] = nVal[0]

            mask = nMasks[i]
            gap_positions = []
            match_positions = []
            for j, v in enumerate(mask):
                gap_positions.append(j) if v == 0 else match_positions.append(j)

            scores[i] += len(gap_positions)*INDEL if len(gap_positions) < MULTI_GAP_THRESHOLD else MULTI_GAP

            if len(match_positions) > 1:
                combos = list(combinations(match_positions, 2))
                for combo in combos:
                    seqA = seqs[combo[0]]
                    seqB = seqs[combo[1]]
                    a = seqA[cur_pos[combo[0]]]
                    b = seqB[cur_pos[combo[1]]]
                    scores[i] += MATCH if a == b else MISMATCH

        maxScore = max(scores)
        maxI = scores.index(maxScore)
        delta = nMasks[maxI]

        matrix[cur_pos] = (maxScore, delta)

        if maxScore > best_cell[0]:
            best_cell = (maxScore, cur_pos)

        if cur_pos == final_pos:
            break
        cur_pos = nextMatrixPosition(cur_pos, dimensions)

    print(f'Alignment score: {matrix[final_pos][0]}')

    path = []
    cur_pos = final_pos
    while(True):
        if cur_pos == [0] * len(dimensions):
            break
        cell = matrix[cur_pos]
        path.append((cell[1], cur_pos))
        cur_pos = vectorSum(cur_pos, cell[1])
    path.reverse()

    #print(path)
    alignments = ["" for s in seqs]
    for node in path:
        delta = node[0]
        for pos, coord in enumerate(delta):
            alignments[pos] += "." if coord == 0 else seqs[coord][node[1][coord]]

    maxlen = 0
    for i, a in enumerate(alignments):
        residue = seqs[i]
        for c in a:
            residue = residue.replace(c, '', 1)
        alignments[i] += residue
        length = len(alignments[i])
        if length > maxlen:
            maxlen = length

    for i, a in enumerate(alignments):
        dif = maxlen - len(a)
        for j in range(dif):
            alignments[i] += "."


    for i, a in enumerate(alignments):
        print(f's{i+1}: {a}')

def main():
    sequences = [
        "A",
        "C",
        "C"
    ]
    MSA(sequences)

if __name__ == '__main__':
    main()
