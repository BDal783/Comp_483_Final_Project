mport argparse

def filter(read: str):
    """

    :param read: _description_
    :return: _description_
    """

    read_length = len(read)
    ambigious_aa = 0
    for letter in read:
        if letter == 'X':
            ambigious_aa += 1
    
    if (ambigious_aa / read_length) > 0.01:
        return True
    else:
        return False

with open('output.txt', 'r') as f:
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

print(len(fastas))



filtered_fastas = {}

for each_read in fastas:
    if not filter(fastas[each_read]):

        filtered_fastas[each_read] = fastas[each_read]
         
with open('filter_reads.txt', 'w') as f: 
    for key in filtered_fastas: