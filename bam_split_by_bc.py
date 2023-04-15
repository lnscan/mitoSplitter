##### This code is used to split bam file by barcode list,
##### extracted from https://www.jianshu.com/p/1cc868654499 and modified by linl.
##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam  -o sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory

### Python 3.6.8
import pysam
import sys
import os

### Input varibles to set
# file to split on
unsplit_file = sys.argv[1]    # input: /path/to/dir/of/data/sorted_tags.bam
# where to place output files
out_dir = sys.argv[2]    # input: /path/to/dir/of/out_data/
# cell barcode tag in bam file
tag = sys.argv[3]    # input: e.g.CB for 10x

# variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile( unsplit_file, "rb")
for read in samfile.fetch( until_eof=True):
    # barcode itr for current read
    CB_itr = read.get_tag(tag)    # input
    # if change in barcode or first line; open new file  
    if( CB_itr!=CB_hold or itr==0):
        # close previous split file, only if not first read in file
        if( itr!=0):
            split_file.close()
        CB_hold = CB_itr
        print('cell barcode: '+CB_hold)
        itr+=1
        out_dir_num = str(int(itr / 10000))
        out_dir_comb = out_dir + out_dir_num + "/"
        split_file = pysam.AlignmentFile(out_dir_comb+"{}.bam".format(CB_hold), "wb", template=samfile)
    # write read with same barcode to file
    split_file.write(read)
split_file.close()
samfile.close()

