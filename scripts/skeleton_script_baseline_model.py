#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    Description:
    Baseline impact predictor of SNPs in VEP format.
    Uses raw BLOSUM62 matrix from a text file for scoring.
"""

import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Generates baseline model predictions.')
    parser.add_argument('vep', help='a path to the VEP input file')
    parser.add_argument('blosum', help='a path to the BLOSUM62 input file')
    parser.add_argument('-o', dest='out_path', help='output .tsv file for baseline model scores', required=True)
    return parser.parse_args()

def parse_blosum(path):
    aas = []
    aa_scores = []
    with open(path, "rb") as f:
        for line in f:
            line = line.decode('UTF-8')
            if line.startswith('#') or line.startswith('x'):
                continue
            else:
                aas.append(line.strip('\n').split()[0])
                aa_scores.append(line.strip('\n').split()[1:len(line)])
    blosum_dict = {}
    for aa in aas:
        blosum_dict[aa] = {}
    # Fill in the BLOSUM62 matrix as a dictionary
    for i in range(len(aas)):
        for j in range(len(aas)):
            score = aa_scores[i][j]
            blosum_dict[aas[i]][aas[j]] = score
    return blosum_dict

def parse_vep(path):
    hgvs_ids = []
    ref_aas = []
    mut_aas = []
    with open(path, "rb") as f:
        for line in f:
            line = line.decode('UTF-8').strip('\n')
            if line.startswith('#'):
                continue
            else:
                hgvs_ids.append(line.split('\t')[0])
                vars = line.split('\t')[1]
                # Get reference and mutated amino acids
                ref_aas.append(vars[0])
                mut_aas.append(vars[2])
    return hgvs_ids, ref_aas, mut_aas

def run_baseline(hgvs_ids, ref_aas, mut_aas, blosum_dict):
    scores = []
    # Calculate BLOSUM62 score for each mutation
    for ref_aa, mut_aa in zip(ref_aas, mut_aas):
        score = blosum_dict[ref_aa][mut_aa]
        scores.append(score)
    return scores

def write_data(hgvs_ids, scores, out_filepath):
    # Write results to file
    with open(out_filepath, 'w') as f:
        f.write('# ID\tScore\n')
        for id, score in zip(hgvs_ids, scores):
            f.write(id + '\t' + score + '\n')

def main():
    args = parse_args()
    vep_path = args.vep
    blosum_path = args.blosum
    out_filepath = args.out_path
    out_dir, out_filename = os.path.split(out_filepath)
    if '.tsv' not in out_filename:
        sys.exit(r'ERROR: filename "%s" in the output file path argument should contain .tsv extension!' % out_filename)
    if not os.path.exists(out_dir):
        sys.exit(r'ERROR: output directory "%s" to store baseline model output does not exist!' % out_dir)
    blosum_dict = parse_blosum(blosum_path)
    hgvs_ids, ref_aas, mut_aas = parse_vep(vep_path)
    scores = run_baseline(hgvs_ids, ref_aas, mut_aas, blosum_dict)
    write_data(hgvs_ids, scores, out_filepath)

if __name__ == "__main__":
    main()
