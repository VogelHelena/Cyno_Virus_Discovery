#!/usr/bin/env python3

from Bio import SeqIO
import argparse

# Initialize the count dictionary
count = {}

# Get the input file name from command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('input', metavar='N', type=str, help='Input fasta')
args = parser.parse_args()

# Define the sequence length
n = 10

# Create a SeqIO object to read the FASTA file
seqio = SeqIO.parse(args.input, "fasta")

for seq in seqio:
    print(f"Analyzing {seq.id}... ", end="")
    s = str(seq.seq)

    m1 = m2 = None
    for c in s:
        c = c.upper()
        if m2 is None:
            m2 = c
            continue
        if m1 is None:
            m1 = c
            continue
        trint = m2 + m1 + c
        dint = m1 + c
        if 'N' in trint:
            continue
        if trint in count:
            count[trint] += 1
            count[dint] += 1
        else:
            count[trint] = 1
            count[dint] = 1
        m2 = m1
        m1 = c

    print("done")

# Print the count dictionary
for key, value in count.items():
    print(f"{key}\t{value}")
