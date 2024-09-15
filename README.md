# matchaSeq 

This is a tool that finds methyltransferase recognition sites within a given sequence. The tool takes a FASTA file as input and a search sequence as well as a context size as command-line arguments.
For each sequence in the FASTA file, the tool will:
- Identify the number of occurrences of the search sequence.
- Report the positions where the search sequence was found.
- Output the sequence context around each occurrence of the search sequence, where the context size is specified by the command-line argument.

The reporting of sequence contexts around each occurrence of the search sequence further informs the researchers of sequence variations between recognition sites. This additional information can be valuable for researchers working with methyltransferase enzymes, as it can provide insights into the specificity and flexibility of the recognition sites.

The output will be txt file containing the following information for each sequence:
- The name of the sequence
- The search sequence used
- The number of occurrences of the search sequence
- The positions where the search sequence was found
- The sequence context around each occurrence of the search sequence

If an output file is not specified in the command line argument, then results are printed out to screen.
```
name: ecorI_test
search string: GACT
number of times found: 2
positions of gact: [5, 88]
sequence contexts found: ['CGGACTGT', 'CGGACTAG']
```

This tool can be useful for researchers working with methyltransferase enzymes, as it allows them to quickly identify potential recognition sites within genomic or transcriptomic sequences.  

## Set up Environment
matchaSeq requires BioPython to run.  To set up a conda environment that contains Biopython, follow the instructions below:

```
# Create new conda env:
conda create -n matchaseq python=3.10
conda activate matchaseq

# install Biopython
conda install -n conda-forge biopython

# verify installation
python -c "import Bio; print(Bio.__version__)"    # should return the version number
1.84

## Usage
Command line usage:
```
% python3 matchaSeq.py -h
usage: matchaSeq.py [-h] [-n NUMBASES] -i INPUT -f FIND [-o OUTPUT]

Given a methylase (methyltransferase) recognition site and the number of upstream and downstream sequence from it, find the
number and locations of recognitions sites in a given DNA template sequence.

optional arguments:
  -h, --help            show this help message and exit
  -n NUMBASES, --numbases NUMBASES
                        number of bases before or after the search string
  -i INPUT, --input INPUT
                        The full path to your fasta file to search.
  -f FIND, --find FIND  The search string or the methylase recognition sequence
  -o OUTPUT, --output OUTPUT
                        path to the output file (optional)
```

### Example output
```
% python3 matchaSeq.py -i fasta_files/matchaseq_test.fasta -f gact -n 2
Processing sequence: dam_test1
Processing sequence: dam_test2
Processing sequence: ecorI_test
Processing sequence: negative_control_test
name: dam_test1
search string: GACT
number of times found: 0
positions of gact: []
sequence contexts found: []

name: dam_test2
search string: GACT
number of times found: 0
positions of gact: []
sequence contexts found: []

name: ecorI_test
search string: GACT
number of times found: 2
positions of gact: [5, 88]
sequence contexts found: ['CGGACTGT', 'CGGACTAG']

name: negative_control_test
search string: GACT
number of times found: 0
positions of gact: []
sequence contexts found: []
```


