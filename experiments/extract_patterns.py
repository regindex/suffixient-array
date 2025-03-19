#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, os, glob, math, platform, re, random

Description = """
Pipeline for extracting patterns out of a text
"""
base_path = os.path.dirname(os.path.abspath(__file__))

random.seed(0)

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input_file', nargs='?', help='the input file',    type=str)
    parser.add_argument('--pattern_number', help='number of patterns',     type=int, required=True)
    parser.add_argument('--pattern_length', help='length of the patterns', type=int, required=True)
    parser.add_argument('--output_file',    help='output file path',       type=str, required=True)
    args = parser.parse_args()

    with open(args.input_file,'r') as text:
        T = text.read()

    with open(args.output_file,'w+') as output:
        for i in range(args.pattern_number):
            output.write(">pattern"+str(i+1)+"\n")
            ri = random.randint(0, len(T)-args.pattern_length)
            output.write(T[ri:ri+args.pattern_length]+"\n")

##########################
if __name__ == '__main__':
    main()