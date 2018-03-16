#!/usr/bin/env python
#############################################################

import argparse

parser = argparse.ArgumentParser(description="Generate FASTQ files with regularly spaced reads from FASTA sequences, and align them to specified genome builds with subread.")
parser.add_argument("-i", "--index", type=str, nargs='+', help="subread genome indexes to align to", required=True)
parser.add_argument("-f", "--fasta", type=str, nargs="+", help="FASTA files to take read sequences from", required=True)
parser.add_argument("-r", "--readlen", type=int, help="read length in base pairs (default: 50)", default=50)
parser.add_argument("-o", "--outdir", type=str, nargs='+', help="directory for output files, one for each index", required=True)
parser.add_argument("-s", "--spacing", type=int, help="spacing between adjacent reads (default: 10)", default=10)
parser.add_argument("-T", "--threads", type=int, help="number of threads to use for subread alignment (default: 6)", default=6)

args = parser.parse_args()
if len(args.outdir)!=len(args.index):
    raise RuntimeError, "need one output directory for each specified index"

#############################################################

from Bio import SeqIO
import os
import subprocess
import tempfile

# Defining a common quality string for Phred +33.
qualstr = "H"*args.readlen 

# Need this to identify ambiguous bases.
import re
pattern = re.compile("(?i)[^ATGC]")

for input_fa in args.fasta:
    fasta_sequences = SeqIO.parse(open(input_fa),'fasta')

    for entry in fasta_sequences:
        name, sequence = entry.id, str(entry.seq)
        short_name = name.split(" ")[0]
        seqlen = len(sequence)
        
        output_handle, output_file = tempfile.mkstemp(dir=".", suffix="_"+short_name+".fastq")
        OHANDLE = os.fdopen(output_handle, "w")
        pos = 0
        endpos = args.readlen
        curstring = sequence[pos:endpos]
        invalid_count = len(re.findall(pattern, curstring))
        
        # Extracting read sequences (ignoring reads with too much ambiguous sequence).
        while endpos < seqlen:
            if invalid_count < args.readlen/5:
                print >> OHANDLE, "@" + str(pos+1) + "\n" + curstring + "\n+\n" + qualstr
            pos += args.spacing
            endpos += args.spacing
            invalid_count -= len(re.findall(pattern, curstring[:args.spacing]))
            curstring = sequence[pos:endpos]
            invalid_count += len(re.findall(pattern, curstring[-args.spacing:]))

        OHANDLE.close()
            
        # Aligning them to specified builds with subread.
        for x, index in enumerate(args.index):
            code = subprocess.call(["subread-align", "-i", index, "-r", output_file, "-o",
                os.path.join(args.outdir[x], short_name+".bam"), "--BAMoutput", 
                "-u", "-H", "-P", str(3), "-T", str(args.threads)])
            if code:
                raise RuntimeError, "subread failed for '"+short_name+"' on '"+index+"'"

        os.remove(output_file)

#############################################################
# End.

