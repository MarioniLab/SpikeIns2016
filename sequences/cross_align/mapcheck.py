#!/usr/bin/env python
#############################################################

import argparse

parser = argparse.ArgumentParser(description="Collate the mapping edit distances of BAM files, for genome-to-genome or transcriptome-to-transcriptome cross-mapping.")
parser.add_argument("-b", "--bam", type=str, nargs="+", help="BAM files to check", required=True)
parser.add_argument("-s", "--source", type=str, help="GTF file containing annotation of source genome", required=True)
parser.add_argument("-a", "--aligned", type=str, help="GTF file containing annotation of aligned genome", required=True)
parser.add_argument("-o", "--out", type=str, help="output file", required=True)
parser.add_argument("-q", "--minqual", type=int, help="minimum MAPQ score to consider an alignment (default: 10)", default=10)

args = parser.parse_args()

#############################################################
# Preparing the GTF objects.

import HTSeq
source_gtf = HTSeq.GFF_Reader(args.source)
source_features = HTSeq.GenomicArrayOfSets("auto", stranded=False)

for feature in source_gtf:
    if feature.type == "exon":
        source_features[feature.iv] += feature.attr["gene_id"]

dest_gtf = HTSeq.GFF_Reader(args.aligned)
dest_features = HTSeq.GenomicArrayOfSets("auto", stranded=False)

for feature in dest_gtf:
    if feature.type == "exon":
        dest_features[feature.iv] += feature.attr["gene_id"]

#############################################################

import os
OHANDLE = open(args.out, mode="w")

import re

for bf in args.bam:
    almnt_file = HTSeq.BAM_Reader(bf)
    cur_chr = re.sub("\\.bam$", "", os.path.basename(bf))
    counts_genome = [0]
    counts_transcriptome = [0]
    counts_transcriptome_edit = [0]
    ref_rlen = -1

    for almnt in almnt_file:
        # Constructing the reference alignment.
        total_rlen = len(almnt.read)
        if ref_rlen < 0:
            ref_rlen = total_rlen
            counts_genome += [0] * total_rlen
            counts_transcriptome += [0] * total_rlen
            counts_transcriptome_edit += [0] * total_rlen
        elif ref_rlen != total_rlen:
            raise RuntimeError, "differing read lengths are not supported"

        ref_pos = int(almnt.read.name)
        s_almnt = HTSeq.Alignment(almnt.read, HTSeq.GenomicInterval(cur_chr, ref_pos, ref_pos+total_rlen, "+"))
        source_ids = set()
        for iv, val in source_features[s_almnt.iv].steps():
            source_ids |= val

        # Checking whether it's actually aligned to the new genome
        if not almnt.aligned or almnt.aQual < args.minqual:
            counts_genome[0] += 1
            if len(source_ids):
                counts_transcriptome[0] += 1
                counts_transcriptome_edit[0] += 1
            continue
        
        # Computing the number of matching bases (edit distance includes indels, so we add them back in)
        edit_distance = almnt.optional_field("NM")
        total_match =- edit_distance
        for cigop in almnt.cigar:
            if cigop.type in "MDI":
                total_match += cigop.size
        counts_genome[total_match] += 1

        # Checking whether they align in both transcriptomes, where we count it as a transcriptomic mismatch
        if len(source_ids): 
            dest_ids = set()
            for iv, val in dest_features[almnt.iv].steps():
                dest_ids |= val

            if len(dest_ids):
                counts_transcriptome[total_match] += 1
            
                # We also use the edit distance - cross-mapping across splice junctions might 
                # have low matches due to soft clipping. Using the total read length to store 
                # the edit distances (i.e., last entry has edit of 0).
                lindex = max(0, total_rlen - edit_distance)
                counts_transcriptome_edit[lindex] += 1
            else:
                counts_transcriptome[0] += 1
                counts_transcriptome_edit[0] += 1

    # Printing diagnostics to file
    print >> OHANDLE, bf + "\tGENOME\tMATCH",
    for y in counts_genome:
        print >> OHANDLE, "\t" + str(y),
    print >> OHANDLE

    print >> OHANDLE, bf + "\tTRANS\tMATCH",
    for y in counts_transcriptome:
        print >> OHANDLE, "\t" + str(y),
    print >> OHANDLE

    print >> OHANDLE, bf + "\tTRANS\tEDIT",
    for y in counts_transcriptome_edit:
        print >> OHANDLE, "\t" + str(y),
    print >> OHANDLE

#############################################################
# Closing handles.

OHANDLE.close()

#############################################################
# End. 
