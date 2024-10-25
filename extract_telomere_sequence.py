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

# https://onestopdataanalysis.com/shannon-entropy/
def estimate_shannon_entropy(dna_sequence):
    m = len(dna_sequence)
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
    
    shannon_entropy_value = 0
    for base in bases:
        n_i = bases[base]
        p_i = n_i / float(m)
        entropy_i = p_i * (math.log(p_i, 2))
        shannon_entropy_value += entropy_i
    
    return round(shannon_entropy_value * -1, 2)

def estimate_composition(dna_sequence, telo_side):
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    dna_sequence = dna_sequence.upper()
    
    if telo_side == "right":
        matches = re.findall(r"TG", dna_sequence)
    else:
        matches = re.findall(r"CA", dna_sequence)

    return len(matches) / (len(dna_sequence) / 2)

def sliding_window_analysis(full_sequence, window_size, telo_side, metric):
    sequence_length = len(full_sequence)
    metric_values = []
    
    if telo_side == "left":
        range_start = sequence_length - window_size
        range_stop = -1
        step = -1
    else:
        range_start = 0
        range_stop = sequence_length - window_size + 1
        step = 1

    for start in range(range_start, range_stop, step):
        window = full_sequence[start:start + window_size]
        if metric == "entropy":
            results_metric = estimate_shannon_entropy(window)
        elif metric == "composition":
            results_metric = estimate_composition(window, telo_side)
        metric_values.append(results_metric)
    
    return metric_values

def trim_soft_clipped_length(full_sequence, window_size, metric_values, metric_thresholds_low, metric_thresholds_high):
    try:
        pos_threshold_break = next(i for i, value in enumerate(metric_values) if value < metric_thresholds_low or value > metric_thresholds_high)
    except StopIteration:
        return len(full_sequence), len(full_sequence)

    if pos_threshold_break == 0:
        return 0, 0
    else:
        return int(round(pos_threshold_break + window_size / 2, 0)), pos_threshold_break

def read_telomeres_csv(file_path):
    telomeres = {}
    with open(file_path, mode='r') as file:
        reader = csv.reader(file, delimiter=',')
        next(reader)  # Skip the header if there is one
        for row in reader:
            if 'term' not in row or 'chr_chrM' in row:
                continue
            print(row)
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
    if telo_side == "left":
        if cigar_tuples[0][0] == 4:  # Soft-clipped is denoted by 4
            return cigar_tuples[0][1]
    else:
        if cigar_tuples[-1][0] == 4:
            return cigar_tuples[-1][1]
    return 0

def main(input_bam_name, input_telomere_name, output_name, is_debug, entropy, composition, save_seq):
    bam_file = pysam.AlignmentFile(input_bam_name, "rb")
    telomeres = read_telomeres_csv(input_telomere_name)
    read_ids = {}
    read_strand = {}
    telomere_lengths = {}
    telomere_sides = {}
    pos_breaks = {}
    length_overhangs = {}
    sequence_overhang = {}
    sequence_telomere = {}
    mapping_qualities = {}

    for chr, regions in telomeres.items():
        read_ids[chr] = []
        read_strand[chr] = []
        telomere_lengths[chr] = []
        telomere_sides[chr] = []
        pos_breaks[chr] = []
        length_overhangs[chr] = []
        sequence_overhang[chr] = []
        sequence_telomere[chr] = []
        mapping_qualities[chr] = []

        if is_debug and chr != "chrX":
            continue

        for start, end in regions:
            telo_side = "left" if end < 10000 else "right"

            if is_debug and telo_side != "left":
                continue

            for read in bam_file.fetch(chr, start, end):
                extension = 0
                if read.is_unmapped:
                    continue
                soft_clipped_length = get_soft_clipped_length(read.cigar, telo_side)
                if read.flag == 0 or read.flag == 16:
                    if telo_side == "left":
                        extension = end - read.reference_start
                        potential_overhang = read.query_sequence[:soft_clipped_length]
                    else:
                        extension = read.reference_end - start
                        potential_overhang = read.query_sequence[len(read.query_sequence)-soft_clipped_length:]

                    pos_interest_ref = end if telo_side == "left" else start

                    pos_interest_seq = -1
                    for alignment_tuple in enumerate(read.get_aligned_pairs()): # with_seq=True
                        if alignment_tuple[1][1] is not None and alignment_tuple[1][0] is not None:
                            if alignment_tuple[1][1] >= pos_interest_ref:
                                pos_interest_seq = alignment_tuple[1][0]
                                break

                    if pos_interest_seq == -1: # No mapping on telomere start
                        continue
                    
                    extracted_seq = read.query_sequence[:pos_interest_seq] if telo_side == "left" else read.query_sequence[pos_interest_seq:]

                    if extension >= 0:
                        if entropy:
                            window_size = 20
                            metric_thresholds = [0.3, 1.4]

                            metric_values = sliding_window_analysis(potential_overhang, window_size, telo_side, "entropy")
                            length_overhang, pos_break = trim_soft_clipped_length(potential_overhang, window_size, metric_values, metric_thresholds[0], metric_thresholds[1])
                        elif composition:
                            window_size = 30
                            metric_thresholds = [0.4, 0.94]

                            metric_values = sliding_window_analysis(potential_overhang, window_size, telo_side, "composition")
                            length_overhang, pos_break = trim_soft_clipped_length(potential_overhang, window_size, metric_values, metric_thresholds[0], metric_thresholds[1])
                        else:
                            length_overhang = len(potential_overhang)

                        if is_debug:
                            print(read.query_name)
                            print(potential_overhang)
                            print(read.cigar)
                            print(extracted_seq)

                        read_ids[chr].append(read.query_name)
                        telomere_lengths[chr].append(extension)
                        telomere_sides[chr].append(telo_side)
                        read_strand[chr].append(read.flag)
                        length_overhangs[chr].append(length_overhang)
                        sequence_overhang[chr].append(potential_overhang if len(potential_overhang) > 0 else -1)
                        sequence_telomere[chr].append(extracted_seq)
                        mapping_qualities[chr].append(read.mapping_quality)

        if is_debug:
            break

    with open(output_name, mode='w', newline='') as tsv_file:
        common_fields = ['Chromosome', 'read_id', 'Strand', 'Type', 'Seq_Full', 'MapQ']
        if is_debug or save_seq:
            fieldnames = common_fields + ['Seq_Overhang']
        else:
            fieldnames = common_fields
        
        tsv_writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
        tsv_writer.writeheader()

        for chr in read_ids.keys():
            for x, a, y, z, c, b, d, e in zip(read_ids[chr], read_strand[chr], telomere_lengths[chr], telomere_sides[chr], length_overhangs[chr], sequence_overhang[chr], sequence_telomere[chr], mapping_qualities[chr]):
                data_row = {'Chromosome': chr, 'read_id': x, 'Strand': a, 'Type': z, 'Seq_Full': d, 'MapQ': e}
                if is_debug or save_seq:
                    data_row['Seq_Overhang'] = b
                tsv_writer.writerow(data_row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process telomere data from BAM and CSV files.")
    parser.add_argument("-i", "--input_bam", required=True, help="Input BAM file path")
    parser.add_argument("-t", "--input_telomere", required=True, help="Input Telomere CSV file path")
    parser.add_argument("-o", "--output_tsv", required=True, help="Output TSV file path")
    parser.add_argument("-e", "--entropy", action="store_true", help="Enable entropy correction")
    parser.add_argument("-c", "--composition", action="store_true", help="Enable composition correction")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debugging mode")
    parser.add_argument("-s", "--save_seq", action="store_true", help="Enable overhang sequence saving")
    args = parser.parse_args()

    main(args.input_bam, args.input_telomere, args.output_tsv, args.debug, args.entropy, args.composition, args.save_seq)

