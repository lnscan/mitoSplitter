#!/usr/bin/env python

import pysam
import argparse
import gzip
parser = argparse.ArgumentParser(description='make fastq from possorted_genome_bam.bam from cellranger')

parser.add_argument('-f', '--bam', required=True, help="cellranger bam")
parser.add_argument('-o', '--out', required=True, help="output fastq name")
parser.add_argument('-c', '--chrom', required = False, help="chrom")
parser.add_argument('-s', '--start', required = False, help="start")
parser.add_argument('-e', '--end', required = False, help="end")
args = parser.parse_args()

assert (not(args.chrom) and not(args.start) and not(args.end)) or (args.chrom and args.start and args.end), "if specifying region, must specify chrom, start, and end"

fn = args.bam#"possorted_genome_bam.bam"#files[0]
bam = pysam.AlignmentFile(fn, "rb")

open_function = lambda f: gzip.open(f,"rt") if f[-3:] == ".gz" else open(f)

if args.chrom:
    bam = bam.fetch(args.chrom, int(args.start), int(args.end))

with open(args.out,'w') as fastq:
    for (index,read) in enumerate(bam):
        if read.is_secondary or read.is_supplementary:
            continue
        if read.seq is None:
            continue
        pos = read.pos
        readname = read.qname
        fastq.write("@"+read.qname+"\n")
        fastq.write(read.seq+"\n")
        fastq.write("+\n")
        fastq.write(read.qual+"\n")

