# cbb752_WatermanSmith

Implementation of the Smith-Waterman Sequence Alignment Algorithm

## Cloning

Run the following code in order to set up the code on a local computer:


```git clone https://github.com/RafiBrent/cbb752_WatermanSmith.git```

```cd SmithWaterman```

## Input/Output files

The code requires two input files, one containing the sequences to be aligned (see sample_input.txt for an example) and another containing the pairwise residue similarity matrix (see blosum62.txt for an example).

The code will produce an output file named sw_output_rbrent.txt (see sample_sw_output_rbrent.txt for an example).

## Optional parameters

The code takes optional parameters specifying the gap extension penalty and the gap opening penalty. These have default values of -1 and -2, respectively.

## Code examples

General usage:

```python3 hw1.py -i <input file> -s <score file> -o <gap opening penalty> -e <gap extension penalty>```

Specific example:

```python3 hw1.py -i sample-input1.txt -s blosum62.txt -o -3 -e -2```

