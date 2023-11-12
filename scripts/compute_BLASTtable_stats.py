#!/usr/bin/env python3

import os
import csv
import gzip
import argparse
import re

parser = argparse.ArgumentParser(prog="compute_BLASTtable_stats.py", 
                                 description="A program to calculate some summary stats about best BLAST hits")
parser.add_argument('-i', "--input", required=True, help="input BLAST table file", type=str)
parser.add_argument('-o', "--output", help="The name of file to save", type=str)
parser.add_argument("--filetype", default="tsv", 
                    choices = ['tsv', 'csv'], 
                    help="The type of output file", type=str)
args = parser.parse_args()

infile = args.input
outfile = ""
if args.output:
    outfile = args.output
else:
    outfile = infile.split(".")[0] + "_stats." + args.filetype

print(f"Input file: {infile}, output file: {outfile}")

with gzip.open(infile, "rt") as infh, open(outfile,"wt") as outfh:
    print("Reading BLAST table file...")
    incsv = csv.reader(infh, delimiter="\t")
    best_hits = {}
    for row in incsv: 
        queryacc = re.sub(r'gi\|\d+\|ref\|([^\|]+)\|', '\\1', row[0])
        subjectacc = re.sub(r'gi\|\d+\|ref\|([^\|]+)\|', '\\1', row[1])
        pid = float(row[2])
        bitscore = float(row[-1])
        # print(f'query is {queryacc} {subjectacc} {pid}')
        if (queryacc not in best_hits) or (bitscore > best_hits[queryacc][2]):
            best_hits[queryacc] = [subjectacc, pid, bitscore]
#        print(row)
#        print(best_hits)
    print(f"Writing output file {outfile}...")
    outcsv = None
    if args.filetype == "tsv":
        outcsv = csv.writer(outfh, delimiter="\t")
    elif args.filetype == "csv":
        outcsv = csv.writer(outfh, delimiter=",")
    else:
        print("Error: unknown file type")
        exit(1)
    average_pid = sum([best_hits[queryacc][1] for queryacc in best_hits])/len(best_hits)
    print(f'Average percent identity: {average_pid:.2f}')
    outcsv.writerow(["queryacc", "subjectacc", "pid", "bitscore"])
    for queryacc in best_hits:
        outcsv.writerow([queryacc, best_hits[queryacc][0], best_hits[queryacc][1], best_hits[queryacc][2]])
    print("Done.")
    