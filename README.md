
# Multiple Sequence Alignment algorithms

Multiple Sequence Alignment algorithms for the Bioinformatics course.

Needleman-Wunsch global alignment and Smith-Waterman local alignment for an arbitrary number of DNA/protein sequences.

## Installation

Make sure to clone the project and initialize the submodule [FastaParser](https://github.com/Kronopt/FastaParser)

```bash
  git clone https://github.com/ViktorHura/Bioinformatics-MSA
  cd Bioinformatics-MSA
  git submodule update --init --recursive
```
## Usage

```bash
python3 main.py config_file.ini input_file.fasta output_file.txt
```

`config_file` must be path to an ini file, containing the parameters for the algorithm in the following format

```ini
[config]
global = yes
match = 5
mismatch = -2
indel = -4
gapgap = 0
```

`input_file` must be path to [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file, containing the sequences you wish to align

`output_file` must be the path to where the output text file will be written (in addition to the console output)

###Examples

Example config files can be found in `configs/` with example inputs in `inputs/` and their generated outputs in `outputs/`

##Optional Improvements

While the algorithms implemented are not intended for practical use due to their complexity, the following improvements can be nice to have when evaluating them.
#### Substitutions matrix

The following function in `main.py` can be easily modified to make use of a substitution matrix:

```python
def replacementScore(A, B):
    return config.getfloat('match') if A == B else config.getfloat('mismatch')
```

#### Just-in-time compilation

The [PyPy](https://www.pypy.org/) JIT-compiler can be used as a drop-in replacement for your default python interpreter. In testing, this resulted in about 5x faster execution time for these algorithms.

```bash
pypy main.py config_file.ini input_file.fasta output_file.txt
```