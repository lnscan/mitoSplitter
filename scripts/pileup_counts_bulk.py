#!/usr/bin/python

###################################################
# summarizes the total number of reads per position
###################################################

import sys
import re
import os
import pysam

bamfile = sys.argv[1] # Filepath to raw bamfile
outpre = sys.argv[2] # Prefix / basename for raw file output
mitofa = sys.argv[3] # Maximum length of mtDNA genome
base_qual = float(sys.argv[4]) # Minimum base quality to be considered for pileup
sample = sys.argv[5] # Sample name to be considered for downstream analyses
alignment_quality = float(sys.argv[6]) # minimum alignment quality required to be considered

# Export Functions
def writeSparseMatrix(mid, vec):
	with open(outpre + "." + mid + ".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec[i] > 0):
				V.write(str(i+1)+" "+str(vec[i])+"\n")


def writeSparseMatrix2(mid, vec1, vec2, mitoseq):
	with open(outpre + "." + mid + ".txt","w") as V:
		for i in range(0,int(maxBP)):
			if(vec1[i] > 0):
				refbase = mitoseq[i]
				V.write(str(i+1)+" "+str(vec1[i])+" "+str(vec2[i])+" "+refbase+"\n")
# Read mito fasta
mitofasta = pysam.FastaFile(mitofa)
maxBP = mitofasta.lengths[0]
mitoseq = mitofasta.fetch(mitofasta.references[0], 0, maxBP)

n = int(maxBP)

# BAQ
# initialize with a pseudo count to avoid dividing by zero
countsA = [0.00000001] * n
countsC = [0.00000001] * n
countsG = [0.00000001] * n
countsT = [0.00000001] * n

qualA = [0.0] * n
qualC = [0.0] * n
qualG = [0.0] * n
qualT = [0.0] * n

bam2 = pysam.AlignmentFile(bamfile, "rb")
for read in bam2:
	seq = read.seq
	quality = read.query_qualities
	align_qual_read = read.mapping_quality
	for qpos, refpos in read.get_aligned_pairs(True):
		if qpos is not None and refpos is not None and align_qual_read > alignment_quality:
			if(seq[qpos] == "A" and quality[qpos] > base_qual):
				qualA[refpos] += quality[qpos]
				countsA[refpos] += 1
			elif(seq[qpos] == "C" and quality[qpos] > base_qual):
				qualC[refpos] += quality[qpos]
				countsC[refpos] += 1
			elif(seq[qpos] == "G" and quality[qpos] > base_qual):
				qualG[refpos] += quality[qpos]
				countsG[refpos] += 1
			elif(seq[qpos] == "T" and quality[qpos] > base_qual):
				qualT[refpos] += quality[qpos]
				countsT[refpos] += 1

meanQualA = [round(x/y,1) for x, y in zip(qualA, countsA)]
meanQualC = [round(x/y,1) for x, y in zip(qualC, countsC)]
meanQualG = [round(x/y,1) for x, y in zip(qualG, countsG)]
meanQualT = [round(x/y,1) for x, y in zip(qualT, countsT)]

countsA = [ int(round(elem)) for elem in countsA ]
countsC = [ int(round(elem)) for elem in countsC ]
countsG = [ int(round(elem)) for elem in countsG ]
countsT = [ int(round(elem)) for elem in countsT ]


# Output Allele Counts
writeSparseMatrix2("A", countsA, meanQualA, mitoseq)
writeSparseMatrix2("C", countsC, meanQualC, mitoseq)
writeSparseMatrix2("G", countsG, meanQualG, mitoseq)
writeSparseMatrix2("T", countsT, meanQualT, mitoseq)

zipped_list = zip(list(countsA),list(countsC),list(countsG),list(countsT))
sums = [sum(item) for item in zipped_list]
writeSparseMatrix("coverage", sums)


