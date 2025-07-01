#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    Description:
    This script draws ROC plots and calculates the AUC score for variant impact predictors.
    Usage:
    python3 skeleton_script_create_roc_plot.py -ibench <benchmark_filepath> -ipred <predictor_filepath> -o <output_filepath> [-color]
"""
import numbers
import os
import sys
import matplotlib
import matplotlib.collections
import matplotlib.pyplot
import argparse
import numpy
matplotlib.use('AGG')

def parse_args():
    parser = argparse.ArgumentParser(description="Draws and saves a ROC plot for one or more impact predictors (SIFT, PolyPhen, BLOSUM62)")
    parser.add_argument("-ipred", "--input_predictor", help="tab-separated file with predictor scores", action='append', required=True)
    parser.add_argument("-ibench", "--input_benchmark", help="tab-separated benchmark classification file", required=True)
    parser.add_argument("-color", "--use_color_roc_plot", help="plot ROC with gradient color", action='store_true', required=False)
    parser.add_argument("-o", "--out_filepath", help="output .png file for ROC plot", required=True)
    return parser.parse_args()

def parse_predictor(filename):
    global type_predictor
    if 'sift' in filename:
        type_predictor = 'sift'
    elif 'polyphen' in filename:
        type_predictor = 'polyphen'
    elif 'baseline' in filename:
        type_predictor = 'BLOSUM'
    else:
        type_predictor = ''
    predictor_dict = {}
    with open(filename, 'r') as f:
        eof = f.seek(0, 2)
        f.seek(0)
        f.readline()
        cur = f.tell()
        if cur == eof:
            sys.exit('ERROR: input predictor file does not contain predictor results!')
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
            if len(arr) != 2:
                print("Warning: line does not have two elements separated by a tab:", line)
            key = arr[0]
            value = float(arr[1])
            predictor_dict[key] = value
    return predictor_dict

def parse_benchmark(filename):
    benchmark_dict = {}
    with open(filename,'r') as f:
        eof = f.seek(0, 2)
        f.seek(0)
        f.readline()
        cur = f.tell()
        if cur == eof:
            sys.exit('ERROR: input benchmark file does not contain benchmark results!')
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")
            if len(arr) < 2:
                print("Warning: line does not have two elements separated by a tab:", line)
            key = arr[0]
            value = arr[1]
            benchmark_dict[key] = value
    return benchmark_dict

def count_total_results(predictor_score_dict, benchmark_dict):
    pathogenic = 0
    benign = 0
    for key, value in predictor_score_dict.items():
        result = benchmark_dict[key]
        if result == 'Pathogenic':
            pathogenic += 1
        elif result == 'Benign':
            benign += 1
    return [benign, pathogenic]

def calculate_coordinates(predictor_score_dict, benchmark_dict, out_filepath):
    score_hgvs_pairs = [(v, k) for k, v in predictor_score_dict.items()]
    if type_predictor != 'polyphen':
        sorted_score_hgvs_pairs = sorted(score_hgvs_pairs)
    else:
        sorted_score_hgvs_pairs = sorted(score_hgvs_pairs, reverse=True)
    coordinate_score = [sorted_score_hgvs_pairs[0][0]]
    tpr = [0.0]
    fpr = [0.0]
    num_benign = 0
    num_pathogenic = 0
    total_benign, total_pathogenic = count_total_results(predictor_score_dict, benchmark_dict)
    index_prebreakpoint_score = []
    previous_score = sorted_score_hgvs_pairs[0][0]
    for i in range(len(sorted_score_hgvs_pairs)):
        score = sorted_score_hgvs_pairs[i][0]
        if previous_score != score:
            index_prebreakpoint_score.append(i - 1)
        previous_score = score
    index_prebreakpoint_score.append(len(sorted_score_hgvs_pairs) - 1)
    for i in range(len(sorted_score_hgvs_pairs)):
        score = sorted_score_hgvs_pairs[i][0]
        hgvs = sorted_score_hgvs_pairs[i][1]
        if benchmark_dict[hgvs] == 'Pathogenic':
            num_pathogenic += 1
        elif benchmark_dict[hgvs] == 'Benign':
            num_benign += 1
        if i in index_prebreakpoint_score:
            tpr.append(num_pathogenic / total_pathogenic)
            fpr.append(num_benign / total_benign)
            coordinate_score.append(score)
    if out_filepath:
        out_dir, out_filename = os.path.split(out_filepath)
        with open(os.path.join(out_dir, out_filename.split('.')[0] + '_xy.tsv'), 'w') as f:
            for a, b in zip(fpr, tpr):
                f.write(str(a) + '\t' + str(b) + '\n')
    return tpr, fpr, coordinate_score

def integrate(fpr, tpr):
    auc = 0.
    last_fpr = fpr[0]
    last_tpr = tpr[0]
    for cur_fpr, cur_tpr in list(zip(fpr, tpr))[1:]:
        delta_x = cur_fpr - last_fpr
        delta_y = cur_tpr - last_tpr
        area_trap = delta_x * last_tpr + (delta_x * delta_y)/2
        auc += area_trap
        last_fpr = cur_fpr
        last_tpr = cur_tpr
    return auc

def roc_plot(tpr, fpr, coordinator_score, out_filepath, color = False):
    auc = integrate(fpr, tpr)
    lw = 1
    figure, axes = matplotlib.pyplot.subplots(1, 1)
    if color:
        if type_predictor != 'polyphen':
            lc = colorline(fpr, tpr, coordinator_score, cmap=matplotlib.pyplot.get_cmap('coolwarm').reversed(), axes=axes)
        else:
            lc = colorline(fpr, tpr, coordinator_score, axes=axes)
        color_bar = figure.colorbar(lc)
        colorbar_legend = type_predictor + ' score'
        color_bar.ax.set_ylabel(colorbar_legend)
    else:
        axes.plot(fpr, tpr)
    axes.plot((0, 1), (0, 1), '--', color='navy', lw=lw, linestyle='--', label='Random')
    axes.set_xlim([-0.008, 1.008])
    axes.set_ylim([-0.008, 1.008])
    axes.set_xlabel('False Positive Rate', fontsize=14)
    axes.set_ylabel('True Positive Rate', fontsize=14)
    if type_predictor == 'polyphen':
        axes.set_title('ROC curve for PolyPhen-2 (AUC = %.3f)' % auc, fontsize=14)
    if type_predictor == 'sift':
        axes.set_title('ROC curve for SIFT (AUC = %.3f)' % auc, fontsize=14)
    if type_predictor == 'BLOSUM':
        axes.set_title('ROC curve for BLOSUM62 (AUC = %.3f)' % auc, fontsize=14)
    matplotlib.pyplot.savefig(out_filepath)

def roc_plot_together(list_tpr, list_fpr, labels, out_filepath):
    lw = 1
    list_color = ['g','r','m']
    figure, axes = matplotlib.pyplot.subplots(1, 1)
    for tpr, fpr, color, label in zip(list_tpr, list_fpr, list_color, labels):
        auc = integrate(fpr, tpr)
        line_label = '{} (AUC= {:.3f})'.format(label, auc)
        axes.plot(fpr, tpr, c=color, label=line_label)
    axes.plot((0, 1), (0, 1), '--', color='navy', lw=lw, linestyle='--', label='Random')
    axes.set_xlim([-0.008, 1.008])
    axes.set_ylim([-0.008, 1.008])
    axes.legend(fontsize=12)
    axes.set_title('ROC curves for different predictors', fontsize=14)
    axes.set_xlabel('False Positive Rate', fontsize=14)
    axes.set_ylabel('True Positive Rate', fontsize=14)
    matplotlib.pyplot.savefig(out_filepath)

def colorline(x, y, z=None, axes=None, cmap=matplotlib.pyplot.get_cmap('coolwarm'), linewidth=3, alpha=1.0, **kwargs):
    def make_segments(x, y):
        points = numpy.array([x, y]).T.reshape(-1, 1, 2)
        segments = numpy.concatenate([points[:-1], points[1:]], axis=1)
        return segments
    if z is None:
        z = numpy.linspace(0.0, 1.0, len(x))
    if isinstance(z, numbers.Real):
        z = numpy.array([z])
    z = numpy.asarray(z)
    segments = make_segments(x, y)
    lc = matplotlib.collections.LineCollection(segments, array=z, cmap=cmap, linewidth=linewidth, alpha=alpha, **kwargs)
    if axes is None:
        axes = matplotlib.pyplot.gca()
    axes.add_collection(lc)
    axes.autoscale()
    return lc

def main():
    args = parse_args()
    predictor_path = args.input_predictor
    benchmark_path = args.input_benchmark
    out_filepath = args.out_filepath
    color = args.use_color_roc_plot
    out_dir, out_filename = os.path.split(out_filepath)
    if '.png' not in out_filename:
        sys.exit(r'ERROR: filename "%s" in the output file path argument should contain .png extension!' % out_filename)
    if not os.path.exists(out_dir):
        sys.exit(r'ERROR: output directory "%s" to store the ROC plot does not exist!' % out_dir)
    if len(predictor_path) == 1:
        predictor_path = predictor_path[0]
        predictor_results = parse_predictor(predictor_path)
        benchmark_results = parse_benchmark(benchmark_path)
        tpr, fpr, coordinate_score = calculate_coordinates(predictor_results, benchmark_results, out_filepath)
        roc_plot(tpr, fpr, coordinate_score, out_filepath, color)
    elif len(predictor_path) != 3:
        sys.exit('ERROR: to plot three predictors (baseline, sift, and polyphen) all together, please input three files by adding -ipred before each file')
    else:
        labels = []
        list_tpr = []
        list_fpr = []
        for predictor in predictor_path:
            predictor_results = parse_predictor(predictor)
            benchmark_results = parse_benchmark(benchmark_path)
            if type_predictor == 'polyphen':
                labels.append('PolyPhen-2')
            elif type_predictor == 'sift':
                labels.append('SIFT')
            elif type_predictor == 'BLOSUM':
                labels.append('BLOSUM62')
            tpr, fpr, coordinate_score = calculate_coordinates(predictor_results, benchmark_results, None)
            list_tpr.append(tpr)
            list_fpr.append(fpr)
        roc_plot_together(list_tpr, list_fpr, labels, out_filepath)

if __name__ == "__main__":
    main()
