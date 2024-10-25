import argparse
import csv
import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

import collections
import math

def read_telomeres_csv(file_path):
    telomeres = {}
    with open(file_path, mode='r') as file:
        reader = csv.reader(file, delimiter=',')
        next(reader)  # Skip the header if there is one
        for row in reader:
            if 'term' not in row or 'chr_chrM' in row:
                continue
            # print(row)
            _, chromosome, _, _, start, end, *_ = row

            if not start or not end:
                continue

            start, end = int(start) - 1, int(end)  # Adjust indexing if necessary
            if chromosome not in telomeres:
                telomeres[chromosome] = []
            telomeres[chromosome].append((start, end))
    for chromosome in telomeres:
        telomeres[chromosome].sort()
    return telomeres

def get_soft_clipped_length(cigar_tuples, telo_side):
    left_clipped = 0
    right_clipped = 0
    
    if cigar_tuples[0][0] == 4:  # Soft-clipped is denoted by 4
        left_clipped = cigar_tuples[0][1]
    
    if cigar_tuples[-1][0] == 4:
        right_clipped = cigar_tuples[-1][1]
    
    return left_clipped, right_clipped

def main(input_bam_name, input_telomere_name, output_name, margin, is_debug):
    bam_file = pysam.AlignmentFile(input_bam_name, "rb")
    telomeres = read_telomeres_csv(input_telomere_name)
    read_ids = {}
    read_strand = {}
    telomere_lengths = {}
    telomere_sides = {}
    aln_length = {}
    left_clipping = {}
    right_clipping = {}

    for chr, regions in telomeres.items():
        read_ids[chr] = []
        read_strand[chr] = []
        telomere_lengths[chr] = []
        telomere_sides[chr] = []
        aln_length[chr] = []
        left_clipping[chr] = []
        right_clipping[chr] = []

        if is_debug and chr != "chrX":
            continue

        for start, end in regions:
            telo_side = "left" if end < 10000 else "right"

            if is_debug and telo_side != "left":
                continue

            if telo_side == "left":
                end = end + margin
            else:
                start = start - margin

            for read in bam_file.fetch(chr, start, end):
                extension = 0
                if read.is_unmapped:
                    continue
                left_clipped, right_clipped = get_soft_clipped_length(read.cigar, telo_side)

                read_ids[chr].append(read.query_name)
                telomere_sides[chr].append(telo_side)
                read_strand[chr].append(read.flag)
                aln_length[chr].append(read.query_alignment_length)
                left_clipping[chr].append(left_clipped)
                right_clipping[chr].append(right_clipped)
                
        if is_debug:
            break

    with open(output_name, mode='w', newline='') as tsv_file:
        common_fields = ['Chromosome', 'read_id', 'Strand', 'Type', 'aln_len', 'left_clipped', 'right_clipped']        
        tsv_writer = csv.DictWriter(tsv_file, fieldnames=common_fields, delimiter='\t')
        tsv_writer.writeheader()

        for chr in read_ids.keys():
            for x, a, z, b, d, e in zip(read_ids[chr], read_strand[chr], telomere_sides[chr], aln_length[chr], left_clipping[chr], right_clipping[chr]):
                data_row = {'Chromosome': chr, 'read_id': x, 'Strand': a, 'Type': z, 'aln_len': b, 'left_clipped': d, 'right_clipped': e}

                tsv_writer.writerow(data_row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process telomere data from BAM and CSV files.")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file path")
    parser.add_argument("-t", "--input_telomere", required=True, help="Input Telomere CSV file path")
    parser.add_argument("-o", "--output_tsv", required=True, help="Output TSV file path")
    parser.add_argument("-m", "--margin", default=0, help="Define margin at telomere")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debugging mode")
    args = parser.parse_args()

    main(args.input_bam, args.input_telomere, args.output_tsv, int(args.margin), args.debug)

