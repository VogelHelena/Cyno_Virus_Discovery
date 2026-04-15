#!/usr/bin/env python3
import argparse
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
from collections import defaultdict

def get_reverse_complement(seq):
    """Get reverse complement of a sequence."""
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp_dict[base] for base in seq[::-1])

def count_pattern(seq, pattern):
    """Count occurrences of a pattern in a sequence."""
    fwd_cnt = len(re.findall(pattern, seq))
    rev_cnt = len(re.findall(pattern, get_reverse_complement(seq)))
    return fwd_cnt+rev_cnt
def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Extract and analyze variant flanking regions from VCF and reference genome'
    )
    parser.add_argument(
        'reference',
        type=argparse.FileType('r'),
        help='Reference genome in FASTA format'
    )
    parser.add_argument(
        'variants',
        type=argparse.FileType('r'),
        help='Variants in VCF format'
    )
    parser.add_argument(
        '-o', '--output',
        default=sys.stdout,
        type=argparse.FileType('w'),
        help='Output file (default: stdout)'
    )
    
    args = parser.parse_args()

    # Create FASTA database
    try:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(args.reference, 'fasta'))
    except Exception as e:
        print(f"Error reading reference file: {e}", file=sys.stderr)
        sys.exit(1)

    # Initialize output string
    header = "#chr\tpos\t5'tetranuc\t3'tetranuc\ttrinuc\tmut\ttrinuc_mut\tstrand\tflank41bp\tCcount\tTCcount\tTCAcount\tTCTcount\tYTCAcount\tRTCAcount\n"
    args.output.write(header)

    # Process variants file
    try:
        for line in args.variants:
            if line.startswith('#'):
                continue
                
            line = line.strip()
            if not line:
                continue

            # Parse VCF format
            fields = line.split('\t')
            if len(fields) < 8:
                continue

            chr_pos = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            info = fields[7]

            # Skip if not single nucleotide variant
            if len(alt) != 1 or len(ref) != 1:
                continue
            if len(fasta_dict[chr_pos]) < pos+20 or pos-21 < 0:
                continue
                # Process based on reference nucleotide
            if ref in ['T', 'C']:
                # Forward strand processing
                flank40 = str(fasta_dict[chr_pos][pos-21:pos+20].seq).upper()
                counts = {
                    'C': count_pattern(flank40, 'C'),
                    'TC': count_pattern(flank40, 'TC'),
                    'TCA': count_pattern(flank40, 'TCA'),
                    'TCT': count_pattern(flank40, 'TCT'),
                    'YTCA': count_pattern(flank40, '[CT]TCA'),
                    'RTCA': count_pattern(flank40, '[GA]TCA')
                }
    
                seq = str(fasta_dict[chr_pos][pos-3:pos+2].seq).upper()
                bases = list(seq)
                seq = bases[1] + "x" + bases[3]
                tetraseq5 = bases[0] + bases[1] + "x" + bases[3]
                tetraseq3 = bases[1] + "x" + bases[3] + bases[4]
                fullstring = f"{bases[1]}[{ref}>{alt}]{bases[3]}"
                strand = "1"

            elif ref in ['A', 'G']:
                # Reverse strand processing
                ref = get_reverse_complement(ref)
                alt = get_reverse_complement(alt)
                flank40 = get_reverse_complement(str(fasta_dict[chr_pos][pos-21:pos+20].seq).upper())
                counts = {
                    'C': count_pattern(flank40, 'C'),
                    'TC': count_pattern(flank40, 'TC'),
                    'TCA': count_pattern(flank40, 'TCA'),
                    'TCT': count_pattern(flank40, 'TCT'),
                    'YTCA': count_pattern(flank40, '[CT]TCA'),
                    'RTCA': count_pattern(flank40, '[GA]TCA')
                }

                seq = get_reverse_complement(str(fasta_dict[chr_pos][pos-3:pos+2].seq).upper())
                bases = list(seq)
                seq = bases[1] + "x" + bases[3]
                tetraseq3 = bases[1] + "x" + bases[3] + bases[4]
                tetraseq5 = bases[0] + bases[1] + "x" + bases[3]
                fullstring = f"{bases[1]}[{ref}>{alt}]{bases[3]}"
                strand = "-1"

            # Add to output string
            args.output.write(f"{chr_pos}\t{pos}\t{tetraseq5}\t{tetraseq3}\t{seq}\t{ref}>{alt}\t{fullstring}\t{strand}\t{flank40}\t{counts['C']}\t{counts['TC']}\t{counts['TCA']}\t{counts['TCT']}\t{counts['YTCA']}\t{counts['RTCA']}\n")

    except Exception as e:
        print(f"Error processing variants file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
