"""
    Description:
    Draws a ROC plot from three .tsv files.

    It should be executed by specifying three input predictor .tsv files from one of the ClinVar datasets
    with HGVS IDs and scores (three predictors: SIFT, PolyPhen, and BLOSUM62).

    -ititle or --plot_title (optional): with this argument followed by a title string, a title provided by the user
    will be added to plot. Otherwise, the default will be used.

    Execute with:
    python3 skeleton_script_plot_tsv.py -itsv <ROCplot_filename_sift>_xy.tsv -itsv <ROCplot_filename_polyphen>_xy.tsv
    -itsv <ROCplot_filename_baseline>_xy.tsv -o <ROCplot_filename_all>_xy.tsv>.png

"""


import argparse
import os
import sys
import matplotlib
import matplotlib.pyplot as plt


def parse_args():
    """
        Parses inputs from the commandline.
        :return: inputs as a Namespace object
    """
    parser = argparse.ArgumentParser(description="Draws and saves a ROC plot for three different predictors "
                                                 "(SIFT, Polyphen, and BLOSUM62) in one figure "
                                                 "using three .tsv files from one of the three datasets")
    parser.add_argument("-itsv", "--tsv_file_of_coordinates", help="tab-separated .tsv file with ROC coordinates. "
                                                                   "This arguments is required!",
                        action='append', required=True)
    parser.add_argument("-ititle", "--plot_title", help="title of the ROC plot (a default str will be provided)",
                        required=False)
    parser.add_argument("-o", "--output_png", help="output .png file path. This arguments is required!", required=True)

    return parser.parse_args()

def parse_tsv_file(filename):
    """
    Parses .tsv file
    :param filename: a str with the input .tsv file name
    :return: list of fpr values and list of tpr values
    """

    global type_vep
    if 'sift' in filename:
        type_vep = 'sift'
    elif 'polyphen' in filename:
        type_vep = 'polyphen'
    elif 'baseline' in filename:
        type_vep = 'BLOSUM'
    else:
        type_vep = ''

    fpr = []
    tpr = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip()
            arr = line.split("\t")

            if len(arr) != 2:
                print("Warning: the following line does not have two elements separated by a tab:\n", line)
            fpr.append(float(arr[0]))
            tpr.append(float(arr[1]))
    return fpr, tpr

def integrate(fpr, tpr):
    """
        Calculate the Area Under the Curve (AUC) for a given list of coordinates.
        :param fpr: a list of fpr
        :param tpr: a list of tpr
        :return: a float with AUC
    """
    auc = 0.
    last_fpr = fpr[0]
    last_tpr = tpr[0]
    for cur_fpr, cur_tpr in list(zip(fpr, tpr))[1:]:
        # Trapezoidal area between points
        delta_x = cur_fpr - last_fpr
        delta_y = cur_tpr - last_tpr
        area_trap = delta_x * last_tpr + (delta_x * delta_y)/2
        auc += area_trap
        last_fpr = cur_fpr
        last_tpr = cur_tpr
    return auc

def roc_plot_together(list_fpr, list_tpr, labels, out_filepath, title):
    """
       Draws a ROC plot for three predictors in one figure without gradient color.
       :param list_tpr: a list of lists with TPRs
       :param list_fpr: a list of lists with FPRs
       :param label: a list with labels for each of the three ROC curves
       :param out_filepath: a str with the output .png file path
       :param title: a str with the title of the plot
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
    axes.legend()
    axes.set_title(title)
    axes.set_xlabel('False Positive Rate')
    axes.set_ylabel('True Positive Rate')
    matplotlib.pyplot.savefig(out_filepath)


def main():
    # Process arguments
    args = parse_args()
    tsv_files = args.tsv_file_of_coordinates
    out_filepath = args.output_png
    title = args.plot_title

    if title is None:
        title = 'ROC of three impact predictors'

    out_dir, out_filename = os.path.split(out_filepath)
    # Check if output filename contains .png extension
    if '.png' not in out_filename:
        sys.exit('Error: filename in the output file path argument should contain .png extension!')

    # Check if output directory exists
    if not os.path.exists(out_dir):
        sys.exit(r'Error: output directory %s to store the ROC plot does not exist! You should create one beforehand.'
                 % out_dir)

    if len(tsv_files) != 3:
        sys.exit('Error: to plot baseline, sift, polyphen all together, '
                 'please input three .tsv files by adding -itsv before each file')
    else:
        labels = []
        list_fpr = []
        list_tpr = []
        for coordinate_list in tsv_files:
            fpr, tpr = parse_tsv_file(coordinate_list)
            labels.append(type_vep)
            list_fpr.append(fpr)
            list_tpr.append(tpr)
        roc_plot_together(list_fpr, list_tpr, labels, out_filepath, title)



if __name__ == "__main__":
    main()
