#!/bin/sh

## test run of mitoSplitter

sh ../mitoSplitter_pipeline.sh -i sc_mix_sorted.bam -b barcodes.list -o example_result -l bam_file.list -g barcodes_cluster.txt
