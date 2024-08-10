import argparse
from Bio import SeqIO


def find_methylase_recognition_sites(fasta_file, search_string, num_bases):
    """ Given a methylase (methyltransferase) recognition site and the number of upstream and downstream sequence from it,
find the number and locations of recognitions sites in a given DNA template sequence. 
    Args:
        fasta_file (str): The full path to the FASTA file to search.
        search_string (str): The search string or the methylase recognition sequence.
        num_bases (int): The number of bases before or after the search string to include.

    Returns:
        tuple: A tuple containing the following elements:
            - name (str): The name of the sequence.
            - num_occurrences (int): The number of occurrences of the search string.
            - positions (list): A list of positions where the search string was found.
            - sequence_contexts (list): A list of sequence contexts surrounding the search string.
            
    """
    results = []

    try:
        for fasta in SeqIO.parse(fasta_file, "fasta"):
            name, sequence = fasta.id, str(fasta.seq.upper())
            print(f"Processing sequence: {name}")

            # count the number of occurrences found of the search string.
            num_occurrences = sequence.count(search_string)
            positions = []
            sequence_contexts = []

            # find all positions of the search string
            start = 0
            while start < len(sequence):
                start = sequence.find(search_string, start)
                if start == -1:
                    break
                positions.append(start)
                
                # Handle the case where the pattern is at the beginning of the sequence:
                if start-num_bases < 0:
                    sequence_contexts.append(sequence[0:start+len(search_string)+num_bases])
                else:
                    sequence_contexts.append(sequence[start-num_bases:start+len(search_string)+num_bases])
                start += len(search_string)   # move start index to the end of the found substring.
        
            results.append((name, num_occurrences, positions, sequence_contexts))
        return results
    except FileNotFoundError:
        print("Error: The input file was not found.")
        exit(1)


# Run main function:
if __name__ == "__main__":
    # Parse the commandline arguments. 
    parser = argparse.ArgumentParser(description = 'Given a methylase (methyltransferase) recognition site and the number of upstream and downstream sequence from it, find the number and locations of recognitions sites in a given DNA template sequence.')
    parser.add_argument('-n', '--numbases', type=int, help='number of bases before or after the search string', default=2)
    parser.add_argument('-i', '--input', help='The full path to your fasta file to search.', required=True)
    parser.add_argument('-f', '--find', help='The search string or the methylase recognition sequence', required=True)
    parser.add_argument('-o', '--output', default=False, help="path to the output file (optional)")
    args = parser.parse_args()

    # Call the function to find the sequence contexts:
    results = find_methylase_recognition_sites(args.input, args.find.upper(), args.numbases)

    # Print or write the results to a file:
    if args.output:
        try:
            with open(args.output, "w") as output_file:
                for name, num_occurrences, positions, sequence_contexts in results:
                    output_file.write(f"name: {name}\n")
                    output_file.write(f"search string: {args.find}\n")
                    output_file.write(f"number of times found: {num_occurrences}\n")
                    output_file.write(f"positions of {args.find}: {positions}\n")
                    output_file.write(f"sequence contexts found: {sequence_contexts}\n\n")
        except OSError as e:
            print(f"Error: Unable to create or write to the output file '{args.output}'. Error: {e}")
            exit(1)
    else:
        for name, num_occurrences, positions, sequence_contexts in results:
            print(f"name: {name}")
            print(f"search string: {args.find}")
            print(f"number of times found: {num_occurrences}")
            print(f"positions of {args.find}: {positions}")
            print(f"sequence contexts found: {sequence_contexts}")
            print()




    