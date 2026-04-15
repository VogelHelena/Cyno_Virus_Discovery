#!/usr/bin/env python3
import argparse

def process_vcf_line(line):
    # Split line into columns
    cols = line.strip().split('\t')
    
    # Get AF value from INFO column
    info_col = cols[7]
    af_value = float([x.split('=')[1] for x in info_col.split(';') if x.startswith('AF=')][0])
    
    # Get DP4 values
    dp4_str = [x.split('=')[1] for x in info_col.split(';') if x.startswith('DP4=')]
    dp4_values = [int(x) for x in dp4_str[0].split(',')] if dp4_str else []
    
    # Process columns based on AF value
    if af_value > 0.5:
        # Swap columns 4 and 5 (REF and ALT)
        cols[3], cols[4] = cols[4], cols[3]
        af_value = 1 - af_value
        
        # Swap first/second with third/fourth values in DP4
        if dp4_values:
            dp4_values[0], dp4_values[1], dp4_values[2], dp4_values[3] = dp4_values[2], dp4_values[3], dp4_values[0], dp4_values[1]
            
            # Update DP4 string in INFO column
            dp4_part = f'DP4={",".join(map(str, dp4_values))}'
            af_part = f'AF={af_value:.6f}'
            info_parts = info_col.split(';')
            for i, part in enumerate(info_parts):
                if part.startswith('DP4='):
                    info_parts[i] = dp4_part
                elif part.startswith('AF='):
                    info_parts[i] = af_part
            cols[7] = ';'.join(info_parts)
    
    return '\t'.join(cols)

def main():
    parser = argparse.ArgumentParser(description='Process VCF file and transform variants based on AF value')
    parser.add_argument('input', help='Input VCF file to process')
    parser.add_argument('-o', '--output', help='Output file for results (default: stdout)', default=None)
    args = parser.parse_args()
    
    # Output results
    if args.output:
        try:
            with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
                for line in infile:
                    # Preserve header lines
                    if line.startswith('#'):
                        outfile.write(line)
                        continue
                    
                    processed_line = process_vcf_line(line)
                    outfile.write(processed_line)
        
            print(f"Successfully processed {args.input} and wrote results to {args.output}")
        
        except FileNotFoundError:
            print(f"Error: Input file '{args.input}' not found")
        except IOError as e:
            print(f"Error: Unable to read/write files - {str(e)}")
        except Exception as e:
            print(f"Error: An unexpected error occurred - {str(e)}")
    else:
        try:
            with open(args.input, 'r') as infile:
                for line in infile:
                    # Preserve header lines
                    if line.startswith('#'):
                        print(line)
                        continue
                    
                    processed_line = process_vcf_line(line)
                    print(processed_line)
        except FileNotFoundError:
            print(f"Error: Input file '{args.input}' not found")
        except Exception as e:
            print(f"Error: An unexpected error occurred - {str(e)}")

if __name__ == '__main__':
    main()
