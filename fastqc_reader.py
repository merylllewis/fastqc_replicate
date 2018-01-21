import sys
import argparse
import matplotlib.pyplot as plt
import numpy as np
import math
import os


def fastq_reader(file_name):
    fastq_read = open(file_name, 'r')
    header = fastq_read.readline().strip("\r\n")
    sequences = []
    quality_scores = []

    while header:
        sequences.append(fastq_read.readline().strip("\r\n"))
        fastq_read.readline().strip("\r\n")
        quality_score_characters = fastq_read.readline().strip("\r\n")
        sequence_quality_scores = [int(ord(score) - 33) for score in quality_score_characters]
        quality_scores.append(sequence_quality_scores)

        header = fastq_read.readline().strip("\r\n")

    return sequences, quality_scores


def per_sequence_quality(sequences, quality_scores, file_name):
    average_q_scores = []
    # y axis: number of reads, x axis: quality scores
    for i in range(len(sequences)):
        sum_q_scores = 0
        for j in range(len(quality_scores[i])):
            sum_q_scores += quality_scores[i][j]
        average_q_scores.append(float(sum_q_scores/len(quality_scores[i])))

    plot_per_sequence_quality(average_q_scores, file_name)


def plot_per_sequence_quality(average_q_scores, file_name):
    bins_array = np.linspace(0, 40, 41)
    fig = plt.figure(figsize = (18, 9))  # width, height
    ax = fig.add_subplot(111)
    bin_frequencies, _, _ = plt.hist(average_q_scores, bins = bins_array, histtype = "step", color = "red", linewidth = 1.5)

    # setting horizontal grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # setting x axis ticks
    xtick_positions = bins_array[::2]
    plt.xticks(xtick_positions)  # xticks positions
    xticks_names = np.array(bins_array[::2], dtype = "int")
    ax.set_xticklabels(xticks_names, fontsize = 20)  # vertical xticks

    # Setting x axis limits
    plt.xlim([0, 42])

    # setting range of y axis as ceiling(max bin frequency), auto-scales the y axis into 10 divisions
    y_range = int(math.ceil(np.amax(bin_frequencies) / 100)) * 100

    yticks_names = np.linspace(0, y_range, 11, dtype = int)  # y axis: 10 ticks
    # setting y axis ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 20)

    ax.set_xlabel("Average quality scores)", fontsize = 20)
    ax.set_ylabel("Number of reads ", fontsize = 20)
    ax.set_title("Quality score distribution over all sequences ", fontsize = 20)
    legend_name = os.path.basename(file_name)[:-len("_plots")]

    #TODO: figure out why legend name is not showing up
    ax.legend(legend_name, loc = "upper left", fontsize = 20)

    plot_path = os.path.join(file_name, "per_sequence_quality.png")
    plt.savefig(plot_path, bbox_inches = "tight")



def main(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fastq_file", type = str, action = "store", help = "input fastq file name")
    args = parser.parse_args(arguments[1:])
    sequences, quality_scores = fastq_reader(args.input_fastq_file)

    plot_directory = os.path.join(os.getcwd(), args.input_fastq_file.split(".")[0] + "_plots")
    if not os.path.exists(plot_directory):
        print "Making new directory for plots"
        os.makedirs(plot_directory)
    per_sequence_quality(sequences, quality_scores, plot_directory)



if __name__ == "__main__":
    main(sys.argv)