import argparse
from Bio import SeqIO

""" Given a methylase (methyltransferase) recognition site and the number of upstream and downstream sequence from it,
find the number and locations of recognitions sites in a given DNA template sequence. """

# pseudo code
# Add arguments for commandline usage.
parser = argparse.ArgumentParser(description = 'Given a methylase (methyltransferase) recognition site and the number of upstream and downstream sequence from it, find the number and locations of recognitions sites in a given DNA template sequence.')
parser.add_argument('-n', '--numbases', type=int, help='number of bases before or after the search string', default=2)
parser.add_argument('-i', '--input', help='The full path to your fasta file to search.', required=True)
parser.add_argument('-f', '--find', help='The search string or the methylase recognition sequence', required=True)
parser.add_argument('-o', '--output', default=False, help="path to the output file (optional)")
args = parser.parse_args()


# unpack the arguments:
fastas = args.input
recseq = args.find.upper()
numbases = args.numbases
outfile = args.output

if args.output == False:
    output_file = None
else:
    try:
        output_file = open(args.output, 'w')
    except OSError as e:
        print(f"Error:  Unable to create or write to the output file '{args.output}'. Error: {e}")
        exit(1)


# Read in sequence template and find matches
for fasta in SeqIO.parse(fastas, "fasta"):
    name, sequence = fasta.id, str(fasta.seq.upper())
    print(f"Processing sequence: {name}")

    # count the number of occurrences found of the search string.
    num_occurrance = sequence.count(recseq)
    pos_list = []
    seqstr_list = []

    # find all positions of the search string
    start = 0
    while start < len(sequence):
        start = sequence.find(recseq, start)
        if start == -1:
            break
        pos_list.append(start)
        
        # Handle the case where the pattern is at the beginning of the sequence:
        if start-numbases < 0:
            seqstr_list.append(sequence[0:start+len(recseq)+numbases])
        else:
            seqstr_list.append(sequence[start-numbases:start+len(recseq)+numbases])
        start += len(recseq)   # move start index to the end of the found substring.


    if output_file is not None:
        try:
            output_file.write(f"name: {name}\n")
            output_file.write(f"search string: {recseq}\n")
            output_file.write(f"positions of {recseq}: {pos_list}\n")
            output_file.write(f"sequence contexts found: {seqstr_list}\n\n")
        except IndexError:
            print(f"Warning: The pattern '{recseq}' was not found in sequence '{name}'.")
    else:
        try:
            print(f"name: {name}")
            print(f"search string: {recseq}")
            print(f"positions of {recseq}: {pos_list}")
            print(f"sequence contexts found: {seqstr_list}")
            print()
        except IndexError:
            print(f"Warning: The pattern '{recseq}' was not found in sequence '{name}'.")


if output_file is not None:
    try:
        output_file.close()
    except OSError as e:
        print(f"Error: Failed to close the output file '{args.output}'. Error: {e}")





    