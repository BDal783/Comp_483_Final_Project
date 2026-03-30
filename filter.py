


import argparse
import sys

def filter(read: str):
    """Filters out any reads with ambigious amino acids
    5% or greater. 

    :param read: _description_
    :return: _description_
    """

    # Intialize read length and number of ambigious amino acids
    read_length = len(read)
    ambigious_aa = 0

    # Count number of X's (ambigious amino acids)
    for letter in read:
        if letter == 'X':
            ambigious_aa += 1
    
    # If there are too many ambigious amino acids return True
    if (ambigious_aa / read_length) >= 0.05:
        return True
    else:
        return False


def main():
    # Set up parser
    parser = argparse.ArgumentParser(description='Runs the AEGIS pipeline')
    parser.add_argument('-i', '--input', help='input file', required=True)
    parser.add_argument('-o', '--output', help='output file', required=True)

    # Define infile
    args = parser.parse_args(sys.argv[1:])
    infile = args.input
    outfile = args.output

    # Collect fasta fil4es in a dictionary
    # With header as key and sequence as value
    with open(infile, 'r') as f:
        header = f.readline().strip()
        fastas = {}
        fasta = ''
        for line in f:
            if not line.startswith('>'):
                fasta += line.strip()
            else:
                fastas[header] = fasta
                header = line.strip()
                fasta = ''
        fastas[header] = fasta

    print(f'Number of sequences before filter: {len(fastas)}')

    # Filter dictionary only including fastas that pass the threshold
    filtered_fastas = {}
    for each_read in fastas:
        if not filter(fastas[each_read]):

            filtered_fastas[each_read] = fastas[each_read]
    
    print(f'Number of sequences after filter: {len(filtered_fastas)}')

    # Write to new fasta file
    with open(outfile, 'w') as f:
        for key in filtered_fastas:
            # Write the FASTA header
            f.write(f"{key}\n")
            
            # Write the sequence in 60-character lines
            sequence = filtered_fastas[key]
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')

if __name__ == '__main__':
    main()
