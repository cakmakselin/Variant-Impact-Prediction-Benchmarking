#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
    Description:
    This script draws ROC plot with/without gradient color and calculates the AUC score.
    Code partially adapted from https://gist.github.com/podshumok/c1d1c9394335d86255b8

    It should be executed by specifying input files and output file path e.g.:
    python3 skeleton_script_create_roc_plot.py -ibench <benchmark_filepath> -ipred <predictor_filepath>
    -o <output_filepath> -color

    -color is an optional argument: ROC plot for a single predictor (SIFT, PolyPhen-2, or BLOSUM62) will show
    threshold scores with gradient color. If the argument is absent, gradient color will not be shown.

    To plot ROC curves for all three predictors in one figure (without gradient color):
    python3 skeleton_script_create_roc_plot.py -ibench <benchmark_filepath> -ipred <sift_scores_filepath>
    -ipred <polyphen_scores_filepath> -ipred <baseline_scores_filepath>  -o <ROCplot_output_filepath>
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
    """
        Parses inputs from the commandline.
        :return: inputs as a Namespace object
    """
    parser = argparse.ArgumentParser(description="Draws and saves a ROC plot for one of OR all the three "
                                                 "impact predictors (SIFT, PolyPhen, and BLOSUM62) to a file")

    # Arguments
    parser.add_argument("-ipred", "--input_predictor", help="tab-separated file with predictor scores. "
                                                            "This argument is required!",
                        action='append', required=True)
    parser.add_argument("-ibench", "--input_benchmark", help="tab-separated benchmark classification file. "
                                                             "This argument is required!", required=True)
    parser.add_argument("-color", "--use_color_roc_plot", help="plot ROC with gradient color", action
    ='store_true', required=False)
    parser.add_argument("-o", "--out_filepath", help="a path to write the output .png file with a ROC plot. "
                                                     "This argument is required!", required=True)

    return parser.parse_args()

def parse_predictor(filename):
    """
        Parses scores of every HGVS ID out of the predictor input file.
        :param filename: a str with the predictor input file
        :return: a dict with HGVS IDs (keys), and the corresponding predictor scores (values)
    """
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
        # Total bytes in the file (end of file)
        eof = f.seek(0, 2)
        # Go to the beginning of the file again
        f.seek(0)
        # Read the first line (should be the header)
        f.readline()
        # Get the current position of the file pointer
        cur = f.tell()
        # If the file doesn't contain predictor results (the header excluded), exit
        if cur == eof:
            sys.exit('ERROR: input predictor file does not contain predictor results!')

        for line in f:
            line = line.rstrip()
            arr = line.split("\t")

            if len(arr) != 2:
                print("Warning: the following line does not have two elements separated by a tab:\n", line)

            key = (arr[0])
            value = float(arr[1])
            predictor_dict[key] = value

    f.close()
    return predictor_dict

def parse_benchmark(filename):
    """
        Parses every HGVS classification out of the benchmark file.
        :param filename: a str with the benchmark input file
        :return: a dict with HGVS IDs (keys), and corresponding benchmark classifications (values)
    """

    benchmark_dict = {}

    with open(filename,'r') as f:
        # Total bytes in the file (end of file)
        eof = f.seek(0, 2)
        # Go to the beginning of the file again
        f.seek(0)
        # Read the first line (should be the header)
        f.readline()
        # Get the current position of the file pointer
        cur = f.tell()
        # If the file doesn't contain predictor results (the header excluded), exit
        if cur == eof:
            sys.exit('ERROR: input benchmark file does not contain benchmark results!')

        for line in f:
            line = line.rstrip()
            arr = line.split("\t")

            if len(arr) < 2:
                print("Warning: the following line does not have three elements separated by a tab:\n", line)

            key = arr[0]
            value = arr[1]
            benchmark_dict[key] = value

    return benchmark_dict

def count_total_results(predictor_score_dict, benchmark_dict):
    """
        Calculates the total number of positives (P), or pathogenic results, and negatives (N), or benign results.
        :param predictor_score_dict: a dict of all predictor scores
        :param benchmark_dict: a dict of benchmark classifications
        :return: a list of ints for the total number of benign and pathogenic results
    """

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
    """
        Calculates coordinates of x and y based on the predictor scores.
        :param predictor_score_dict: a dictionary with scores produced by parse_predictor()
        :param benchmark_dict: a dictionary with benchmark classifications produced by parse_benchmark()
        :param out_filepath: a str with the output .png file path
        :return: lists of coordinates for the ROC plot (TPR and FPR), and a list of sorted predictor scores
    """

    # Get a list of tuples from predictor_score_dict: (predictor score, HGVS ID)
    score_hgvs_pairs = [(v, k) for k, v in predictor_score_dict.items()]

    sorted_score_hgvs_pairs = score_hgvs_pairs
    #########################
    ### START CODING HERE ###
    #########################
    # You need to sort the scores in the correct order for the ROC plot.
    # Use the following if-statement and replace the question mark with the type of the predictor.
    # It will put the ROC curve at the correct side of the diagonal line.

    # '!=' not was '==' in original code! But makes no sense.
    # We want most deleterious mutations first, as we are starting form the (0, 0) coordinates.
    # Therefore, we want the first values in our list to be the ones that are most likely to be deleterious.
    # These are the highest values in the case of polyphen, and the lower values in the case of sift/baseline
    if type_predictor != 'polyphen':
        sorted_score_hgvs_pairs = sorted(score_hgvs_pairs)
    else:
        sorted_score_hgvs_pairs = sorted(score_hgvs_pairs, reverse=True)

    #########################
    ###  END CODING HERE  ###
    #########################
    # Later, each coordinate in the ROC plot will be associated with a predictor score (a threshold score). Thus, we
    # need a separate list for predictor scores
    coordinate_score = [sorted_score_hgvs_pairs[0][0]]

    # Create lists to store coordinates. Starts in (0,0)
    tpr = [0.0]
    fpr = [0.0]

    # Create variables to keep track of the number of true positives (TPs) and false positives (FPs)
    num_benign = 0
    num_pathogenic = 0

    # Get the total number of positives (P) and negatives (N) (pathogenic and benign)
    total_benign, total_pathogenic = count_total_results(predictor_score_dict, benchmark_dict)

    # Get a list of indices of scores before breakpoints
    index_prebreakpoint_score = []
    previous_score = sorted_score_hgvs_pairs[0][0]
    for i in range(len(sorted_score_hgvs_pairs)):
        score = sorted_score_hgvs_pairs[i][0]
        if previous_score != score:
            # Add index of the score before the breakpoint
            index_prebreakpoint_score.append(i - 1)
        previous_score = score

    # Add index of the last score (for the last coordinate)
    index_prebreakpoint_score.append(len(sorted_score_hgvs_pairs) - 1)

    # Iterate over HGVS IDs of SNPs and corresponding sorted predictor scores
    for i in range(len(sorted_score_hgvs_pairs)):
        score = sorted_score_hgvs_pairs[i][0]
        hgvs = sorted_score_hgvs_pairs[i][1]

        #########################
        ### START CODING HERE ###
        #########################
        # Determine whether the SNP is classified by the benchmark as:
        #    Pathogenic -> actual positive, thus a true positive (x-coordinate)
        #    Benign  -> actual negative, thus a false positive (y-coordinate)
        # The comment above is incorrect right? Because true positive
        # rate is on the y-axis and false positive rate on the x axis
        if benchmark_dict[hgvs] == 'Pathogenic':
            num_pathogenic += 1
        elif benchmark_dict[hgvs] == 'Benign':
            num_benign += 1

        # Increase the respective value of num_benign or num_pathogenic

        # Now, you need to calculate TPR and FPR for unique scores as TP/P and FP/N, respectively,
        # using num_benign, num_pathogenic, total_benign, and total_pathogenic correctly. Append the values
        # to the corresponding lists: tpr is a list of y-coordinates and fpr is a list of x-coordinates.
        # Calculate the rates if HGVS score index i is the index of the score before a breakpoint
        # (use index_prebreakpoint_score). Also, append the score to coordinate_score.
        if i in index_prebreakpoint_score:
            tpr.append(num_pathogenic / total_pathogenic)
            fpr.append(num_benign / total_benign)
            coordinate_score.append(score)
        #########################
        ###  END CODING HERE  ###
        #########################

    if out_filepath:
        out_dir, out_filename = os.path.split(out_filepath)
        # Write coordinates to a .tsv file
        with open(os.path.join(out_dir, out_filename.split('.')[0] + '_xy.tsv'), 'w') as f:
            for a, b in zip(fpr, tpr):
                f.write(str(a) + '\t' + str(b) + '\n')

    return tpr, fpr, coordinate_score

def integrate(fpr, tpr):
    """
        Calculates the Area Under the Curve (AUC) for a given list of coordinates.
        :param fpr: a list of FPRs
        :param tpr: a list of TPRs
        :return: a float with AUC
    """

    auc = 0.
    last_fpr = fpr[0]
    last_tpr = tpr[0]

    for cur_fpr, cur_tpr in list(zip(fpr, tpr))[1:]:
        #########################
        ### START CODING HERE ###
        #########################
        delta_x = cur_fpr - last_fpr
        delta_y = cur_tpr - last_tpr
        area_trap = delta_x * last_tpr + (delta_x * delta_y)/2
        auc += area_trap
        #########################
        ###  END CODING HERE  ###
        #########################
        last_fpr = cur_fpr
        last_tpr = cur_tpr
    return auc

def roc_plot(tpr, fpr, coordinator_score, out_filepath, color = False):
    """
       Draws ROC plot with gradient color.
       :param tpr: a list of TPRs
       :param fpr: a list of FPRs
       :param coordinator_score: a list of predictor scores
       :param out_filepath: a str with the output .png file path
       :param color: boolean (False by default) to enable gradient color plotting
    """

    # Compute AUC
    auc = integrate(fpr, tpr)

    # Draw ROC plot and write it to a file
    lw = 1
    figure, axes = matplotlib.pyplot.subplots(1, 1)
    # matplotlib.pyplot.get_cmap('coolwarm').reversed()
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
    """
       Draws ROC plot for three predictors in one figure without gradient color.
       :param list_tpr: a list of lists with TPRs for each predictor
       :param list_fpr: a list of lists with FPRs for each predictor
       :param labels: a list with labels for each of the three ROC curves
       :param out_filepath: a str with output .png file path
    """

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
    """
        Plots a colored line with coordinates x and y. Optionally, specify colors in the array z. Optionally,
        specify a colormap, a norm function and a line width.
        :param x: a list of x-coordinates
        :param y: a list of y-coordinates
    """

    def make_segments(x, y):
        """
            Creates a list of line segments from x- and y-coordinates in the correct format for LineCollection:
            an array of the form numlines x (points per line) x 2 (x and y) array.
            :param x: a list of x-coordinates
            :param y: a list of y-coordinates
            :return: a list of line segments
        """

        points = numpy.array([x, y]).T.reshape(-1, 1, 2)
        segments = numpy.concatenate([points[:-1], points[1:]], axis=1)

        return segments

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = numpy.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if isinstance(z, numbers.Real):
        z = numpy.array([z])

    z = numpy.asarray(z)

    segments = make_segments(x, y)
    lc = matplotlib.collections.LineCollection(segments, array=z, cmap=cmap, linewidth=linewidth,
                                               alpha=alpha, **kwargs)

    if axes is None:
        axes = matplotlib.pyplot.gca()

    axes.add_collection(lc)
    axes.autoscale()

    return lc


def main():

    # Process arguments
    args = parse_args()
    predictor_path = args.input_predictor
    benchmark_path = args.input_benchmark
    out_filepath = args.out_filepath
    color = args.use_color_roc_plot

    out_dir, out_filename = os.path.split(out_filepath)
    # Check if output filename contains .png extension
    if '.png' not in out_filename:
        sys.exit(r'ERROR: filename "%s" in the output file path argument should contain .png extension!' % out_filename)

    # Check if output directory exists
    if not os.path.exists(out_dir):
        sys.exit(r'ERROR: output directory "%s" to store the ROC plot does not exist! Follow instructions in'
                 r' the manual!' % out_dir)

    # Parse input files, calculate ROC coordinates and plot
    if len(predictor_path) == 1:
        predictor_path = predictor_path[0]
        # Parse predictor and benchmark files
        predictor_results = parse_predictor(predictor_path)
        benchmark_results = parse_benchmark(benchmark_path)
        # Calculate ROC coordinates
        tpr, fpr, coordinate_score = calculate_coordinates(predictor_results, benchmark_results, out_filepath)
        # Draw and save the ROC plot

        roc_plot(tpr, fpr, coordinate_score, out_filepath, color)

    elif len(predictor_path) != 3:
        sys.exit('ERROR: to plot three predictors (baseline, sift, and polyphen) all together, '
              'please input three files by adding -ipred before each file')

    else:
        # Lists to store labels and coordinates for three predictors
        labels = []
        list_tpr = []
        list_fpr = []
        # For each predictor
        for predictor in predictor_path:
            # Parse predictor and benchmark files
            predictor_results = parse_predictor(predictor)
            benchmark_results = parse_benchmark(benchmark_path)
            # Append predictor type to the label list
            if type_predictor == 'polyphen':
                labels.append('PolyPhen-2')
            elif type_predictor == 'sift':
                labels.append('SIFT')
            elif type_predictor == 'BLOSUM':
                labels.append('BLOSUM62')

            # Calculate ROC coordinates
            tpr, fpr, coordinate_score = calculate_coordinates(predictor_results, benchmark_results, None)
            # Append coordinates to lists
            list_tpr.append(tpr)
            list_fpr.append(fpr)
        # Plot ROC curves for three predictors together
        roc_plot_together(list_tpr, list_fpr, labels, out_filepath)


if __name__ == "__main__":
    main()
