#!/bin/bash

# Usage: ./wonderdump.sh --split-files --gzip SRRxxxxxxx

if [[ "$1" == "--split-files" && "$2" == "--gzip" ]]; then
    accession="$3"
    echo "Downloading $accession..."
    fasterq-dump --split-files "$accession" -O .
    echo "Compressing..."
    gzip "${accession}_1.fastq" "${accession}_2.fastq"
    echo "Done: $accession"
else
    echo "Usage: $0 --split-files --gzip <SRA_ID>"
fi
