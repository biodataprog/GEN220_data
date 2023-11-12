#!/usr/bin/env python3

import os
import gzip

import argparse
import re
from Bio import SeqIO

parser = argparse.ArgumentParser(prog="compute_AA_stats.py", 
                                 description="A program to calculate some summary stats about proteins")
parser.add_argument('-i', "--input", required=True, help="input protein file", type=str)
parser.add_argument('-o', "--output", help="The name of file to save", type=str)
parser.add_argument("--format", default="fasta", 
                    choices = ['fasta', 'genbank'], 
                    help="The type of output file", type=str)
args = parser.parse_args()

infile = args.input
outfile = ""
if args.output:
    outfile = args.output
else:
    outfile = os.path.basename(infile).split(".")[0] + "_stats.txt"
    
if infile.endswith(".gz"):
    infh = gzip.open(infile, "rt")
else:
    infh = open(infile, "rt")

print(f"Input file: {infile}, output file: {outfile}")
if not infh:
    print("Had trouble opening input file {infile}" )
    exit()

lengths = []
AA_counts = {}
with open(outfile,"wt") as outfh:
    outfh.write(f"Processing data from: {os.path.basename(infile)}\n")
    for seq_record in SeqIO.parse(infh, args.format):
        # print(seq_record.id, len(seq_record))
        lengths.append(len(seq_record))
        for AA in seq_record.seq:
            if AA not in AA_counts:
                AA_counts[AA] = 0
            AA_counts[AA] += 1
    # print(AA_counts)
    total_AA_number = sum(AA_counts.values())
    outfh.write(f"Total number of AAs is {total_AA_number}\n")
    outfh.write(f'AA\tPercentage\n')
    for AA in sorted(AA_counts.keys()):
        outfh.write(f'{AA}\t{100 * AA_counts[AA]/total_AA_number:.2f}\n')
    outfh.write(f"Avg length: {sum(lengths)/len(lengths):.2f}\n")
