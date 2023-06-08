import sys
import configparser
from functools import partial
import fastaparser
from MSA import MSA


def parseConfig(path):
    # Create configparser object
    config_object = configparser.ConfigParser()
    with open(path, "r") as file_object:
        config_object.read_file(file_object)
        return config_object['config']


def parseInput(path):
    with open(path) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        return [seq for seq in parser]


def parseGlobalOutput(input, output):
    score, alignments = output
    outstr = f'Score: {score}\n'

    maxIDlen = max([len(i.id) for i in input])
    for i, seq in enumerate(input):
        input[i] = seq.id.ljust(maxIDlen, ' ')

    for i, a in enumerate(alignments):
        outstr += f'{input[i]}: {a}\n'
    return outstr


def parseLocalOutput(input, output):
    score, alignments, lpad, rpad = output

    maxIDlen = max([len(i.id) for i in input])
    for i, seq in enumerate(input):
        input[i] = seq.id.ljust(maxIDlen, ' ')

    # indicator line to show which region was aligned
    outstr = f'Score: {score}\n'
    outstr += " "*(maxIDlen+2+lpad) + "[" + " "*(len(alignments[0])-rpad-lpad-2) + "]\n"
    maxj = len(alignments[0]) - rpad            # aligned region = [lpad, maxj[
    for i, a in enumerate(alignments):
        repr = f'{input[i]}: {a}'                 # print aligned string
        if i != len(alignments)-1:              # print | matches for the next alignment
            next_a = alignments[i+1]
            repr += "\n" + " "*(maxIDlen+2) + (" " * lpad)
            for j in range(lpad, maxj):
                repr += "|" if a[j] == next_a[j] and a[j] != '.' else " "
        outstr += repr + "\n"
    return outstr


def main():
    configPath = sys.argv[1]
    inputPath = sys.argv[2]
    outputPath = sys.argv[3]

    config = parseConfig(configPath)
    input = parseInput(inputPath)

    # could modify this to use a substitution matrix
    def replacementScore(A, B):
        return config.getfloat('match') if A == B else config.getfloat('mismatch')

    output = MSA([seq.sequence_as_string() for seq in input],
                 gap_penalty=config.getfloat('indel'),
                 gap_gap_penalty=config.getfloat('gapgap'),
                 global_align=config.getboolean('global'),
                 substitution_func=replacementScore)

    outputString = parseGlobalOutput(input, output) if config.getboolean('global')\
        else parseLocalOutput(input, output)

    print("=== Output ===")
    print(outputString)

    f = open(outputPath, "w")
    f.write(outputString)
    f.close()


if __name__ == '__main__':
    main()
