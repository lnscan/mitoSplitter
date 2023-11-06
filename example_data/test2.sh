#!/bin/bash

## test run of mitoSplitter, all cells in the example_data are singlets

bash ../mitoSplitter_pipeline.sh -i sc_mix.bam -b barcodes.list -o example_result -l bam_file.list -g barcodes_cluster.txt test.log 2>&1 &
