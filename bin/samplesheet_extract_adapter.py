#!/usr/bin/env python3

"""
Author: Elisabet Thomsen
Date: 03-10-2019
Description: Take adapter1 and adapter2 sequences out of metadata.csv and write output to stdout in FASTA format.
Usage: <scriptname.py> <samplesheet_inputfile>
"""

import sys

# Get filename from commandline
filename = sys.argv[1]

infile = open(filename,'r')

adapter1 = -1
adapter2 = -1

# Get adapter sequences
for line in infile:
        if line.startswith("Adapter,"):
                line = line[:-1]
                adapter1 = line.split(",")[1]
        if line.startswith("AdapterRead2"):
                line = line[:-1]
                adapter2 = line.split(',')[1]

infile.close()

assert adapter1 != -1, "Adapter not found in samplesheet: %s" %(filename)
assert adapter2 != -1, "AdapterRead2 not found in samplesheet: %s" %(filename)

adapter_fasta_str = ">Adapter1" "\n" + adapter1 + "\n" + ">Adapter2" + "\n" + adapter2

print(adapter_fasta_str)

