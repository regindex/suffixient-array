#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, os, glob, math, platform, re, random

Description = """
Pipeline for filtering non DNA characters from the dataset
"""
base_path = os.path.dirname(os.path.abspath(__file__))

random.seed(0)

def replace_n_with_random(dna_sequence):
    nucleotides = "ATCG"
    return "".join(random.choice(nucleotides) if base == "N" else base for base in dna_sequence)

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_file', nargs='?', help='the input file', type=str)
    parser.add_argument('--replace_N',  help='replace N characters with a random nucleotide', action='store_true')
    parser.add_argument('--trim_N',  help='remove all N characters', action='store_true')
    parser.add_argument('--trim_nonDNA',  help='remove all non DNA characters', action='store_true')
    parser.add_argument('--trim_newline',  help='remove newline characters', action='store_true')
    parser.add_argument('--output_file',  help='output file path', type=str, required=True)
    args = parser.parse_args()

    with open(args.input_file,'r') as text:
        T = text.read()

    if args.trim_newline:
        T = T.replace("\n", "")
    elif args.replace_N:
        T = replace_n_with_random(T)
    elif args.trim_N:
        T = T.replace("N", "")
    elif args.trim_nonDNA:
        valid_dna_chars = {'A','T','C','G'}
        T = ''.join([char for char in T if char in valid_dna_chars])

    with open(args.output_file,'w+') as output:
        output.write(T)

##########################
if __name__ == '__main__':
    main()