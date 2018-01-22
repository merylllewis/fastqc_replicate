import os
import numpy as np
import matplotlib.pyplot as plt
import math


def plot_per_sequence_quality(average_q_scores, file_name):
    """
    Plots average quality score distribution. Y axis: Number of reads for a particular quality score,
    X axis: average quality scores

    :param average_q_scores: same as above
    :param file_name: same as above
    :return: nothing

    """

    # Find the highest average quality score - this is used to decide the x-axis limits
    max_avg_quality_score = max(average_q_scores)

    # Create a histogram to get frequency of each average quality score value
    histogram_bins = np.linspace(0, max_avg_quality_score, max_avg_quality_score)
    bin_frequencies, _ = np.histogram(average_q_scores, bins = histogram_bins)
    # _ : throwaway variable. removing it would change size of bin_frequencies

    # Create a new figure and add a subplot
    fig = plt.figure(figsize = (18, 9))  # width, height
    ax = fig.add_subplot(111)

    # Create a line plot of number of reads versus average quality scores
    plt.plot(histogram_bins[:-1], bin_frequencies)

    # Setting x axis grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # Setting x axis ticks
    xtick_positions = histogram_bins[::2]
    plt.xticks(xtick_positions)  # xticks positions
    xticks_names = np.array(histogram_bins[::2], dtype = "int")
    ax.set_xticklabels(xticks_names, fontsize = 16)

    # Setting x axis limits
    plt.xlim([0, max_avg_quality_score])

    # Setting range of y axis as ceiling (max of bin frequency), auto-scales the y axis into 10 divisions
    y_range = int(math.ceil(np.amax(bin_frequencies) / 100)) * 100

    # Setting y axis ticks
    yticks_names = np.linspace(0, y_range, 11, dtype = int)  # y axis: 10 ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 16)

    # Setting labels, title, legend
    ax.set_xlabel("Average quality scores", fontsize = 16)
    ax.set_ylabel("Number of reads ", fontsize = 16)

    ax.set_title("Quality score distribution over all sequences ", fontsize = 16)

    legend_name = str(os.path.basename(file_name)[:-len("_plots")])
    ax.legend([legend_name], loc = "upper left", fontsize = 16)

    # Saving plot to disk
    print "Plotting per sequence quality"
    plot_path = os.path.join(file_name, "per_sequence_quality.png")
    plt.savefig(plot_path, bbox_inches = "tight")

    plt.close()


def plot_per_base_sequence_content(group_by_base, plot_directory):
    """
    Plots the %occurence of each base at each position. X axis: position in a read, Y axis: percentage

    :param group_by_base: dictionary of % occurrence of each base at each position. key: base, values: %of the base at
                          each position in a read
    :param plot_directory: Name of the directory to save plots
    :return: nothing
    """
    color_palette = ["blue", "black", "red", "green", "violet"]    # Setting colors for each base
    number_positions = len(group_by_base["A"])   # number_positions = length of each sequence

    fig = plt.figure(figsize = (18, 9))  # width, height
    ax = fig.add_subplot(111)

    x_range = np.linspace(1, number_positions, number_positions, dtype = "int")

    # Plotting percentage of each base at each position
    plt.plot(x_range, group_by_base["A"], color = color_palette[0])
    plt.plot(x_range, group_by_base["C"], color = color_palette[1])
    plt.plot(x_range, group_by_base["G"], color = color_palette[2])
    plt.plot(x_range, group_by_base["T"], color = color_palette[3])
    plt.plot(x_range, group_by_base["N"], color = color_palette[4])

    # Setting x axis grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # setting x axis ticks
    plt.xticks(x_range[::5])  # xticks positions
    ax.set_xticklabels(x_range[::5], fontsize = 16)

    # Setting x axis limits
    plt.xlim([1, number_positions])

    # setting range of y axis as 100 %
    yticks_names = np.linspace(0, 100, 11, dtype = int)  # y axis: 10 ticks
    # setting y axis ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 16)

    # Setting labels, title and legend
    ax.set_xlabel("Position in a read)", fontsize = 16)
    ax.set_ylabel("Percentage occurrence ", fontsize = 16)
    ax.set_title("Sequence content across all bases ", fontsize = 16)
    legend_name = ["A", "C", "G", "T", "N"]
    ax.legend(legend_name, loc = "upper left", fontsize = 16)

    # Saving to disk
    print "Plotting per base sequence content"
    plot_path = os.path.join(plot_directory, "per_base_sequence_content.png")
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def gc_content(cg_bases_by_position, x_range, plot_directory):
    """
    Plots the % GC content at each position in a read. X axis: position in a read, Y axis: %GC content

    :param cg_bases_by_position: list of the % GC content at each position in reads
    :param x_range: range for the x axis
    :param plot_directory: Name of the directory to save plots
    :return: nothing
    :return:
    """

    fig = plt.figure(figsize = (18, 9))
    ax = fig.add_subplot(111)

    plt.plot(x_range, cg_bases_by_position)
    plot_path = os.path.join(plot_directory, "gc_content.png")

    # Setting x axis ticks
    plt.xticks(x_range[::5])  # xticks positions
    ax.set_xticklabels(x_range[::5], fontsize = 16)

    # Setting x axis grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # Setting x axis limits
    plt.xlim([1, x_range[-1]])

    # setting range of y axis as 100 %
    yticks_names = np.linspace(0, 100, 11, dtype = int)
    # setting y axis ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 16)

    # Setting labels, title and  legend
    ax.set_xlabel("Position in a read (bp)", fontsize = 16)
    ax.set_ylabel("Percentage occurrence ", fontsize = 16)
    ax.set_title("GC content across all positions in reads ", fontsize = 16)
    legend_name = ["GC_content"]
    ax.legend(legend_name, loc = "upper left", fontsize = 16)

    # Saving plot to disk
    print "Plotting %GC content"
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def plot_per_base_sequence_quality(quality_scores_over_position, plot_directory):
    """
    Plots the quality score distribution (boxplot) over each position in the reads. X axis: position in read, Y axis:
    quality scores

    :param quality_scores_over_position: list of list of all quality scores in a position
    :param plot_directory: Name of the directory to save plots
    :return: nothing

    """

    # Create figure for plot and add subplot
    fig = plt.figure(figsize = (18, 9))
    ax = fig.add_subplot(111)

    # Create a boxplot of quality scores over position
    plt.boxplot(quality_scores_over_position, showfliers = False)    # Do not show points outside box plots

    # Setting x axis ticks
    x_range = np.linspace(1, len(quality_scores_over_position), len(quality_scores_over_position), dtype = "int")
    plt.xticks(x_range[::2])  # setting xticks positions
    ax.set_xticklabels(x_range[::2], fontsize = 12)

    # Setting x axis limits
    plt.xlim([0, x_range[-1]])

    # Setting x axis grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # Setting y axis ticks
    max_q_score = max(max(quality_scores_over_position))
    yticks_names = np.linspace(0, max_q_score, max_q_score+1, dtype = int)  # y axis: 10 ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 12)

    # Setting labels, title
    ax.set_xlabel("Position in a read (bp)", fontsize = 16)
    ax.set_ylabel("Quality scores ", fontsize = 16)
    ax.set_title("Quality scores across all bases ", fontsize = 16)

    # Saving plot to disk
    print "Plotting per base sequence quality"
    plot_path = os.path.join(plot_directory, "per_base_sequence_quality.png")
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def plot_per_sequence_gc_content(gc_content_by_sequence, plot_directory):
    """
    Plots Number of reads v/s mean GC content. X axis: Mean GC content, Y axis: Number of reads

    :param gc_content_by_sequence: list of %GC content of all sequences
    :param plot_directory: Name of the directory to save plots
    :return: nothing
    """

    # Create a new figure and add a subplot
    fig = plt.figure(figsize = (18, 9))  # width, height
    ax = fig.add_subplot(111)

    # Calculate frequency of occurrence of each GC content bin
    x_range = np.linspace(0, 100, 101, dtype = "int")   # Each bin = 1%
    bin_frequencies, _ = np.histogram(gc_content_by_sequence, bins = x_range)
    # _ : throwaway variable to take care of size of bin_frequencies

    # Plot line graph of number of reads vs GC content %
    plt.plot(x_range[:-1], bin_frequencies)

    # Setting x tick positions
    plt.xticks(x_range[::2])  # xticks positions
    ax.set_xticklabels(x_range[::2], fontsize = 12)

    # Setting x axis limits
    plt.xlim([1, x_range[-1]])

    # Setting x axis grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # Setting range of y axis as ceiling (max bin frequency), auto-scales the y axis
    y_range = math.ceil(float(np.amax(bin_frequencies)) / 100) * 100

    # Setting y axis ticks
    yticks_names = np.linspace(0, y_range, 11, dtype = int)  # y axis: 11 ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 12)

    # Setting labels, title
    ax.set_xlabel("Mean GC content (%)", fontsize = 16)
    ax.set_ylabel("Number of reads ", fontsize = 16)
    ax.set_title("GC distribution over all sequences ", fontsize = 16)

    # Saving plot to disk
    print "Plotting number of reads v/s mean %GC content"
    plot_path = os.path.join(plot_directory, "per_sequence_gc_content.png")
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def plot_sequence_duplication(duplication_counts, deduplicated_counts, remainder_after_deduplication, plot_directory):
    """
    Plots the percentage distribution of all duplicated reads out of total reads and out of total unique reads
    X axis: number of duplicated reads, Y axes: Percentage of total reads, Percentage of unique reads

    :param duplication_counts: percentage of duplicated reads versus duplication level
    :param deduplicated_counts: percentage of unique reads versus duplication level
    :param remainder_after_deduplication: number of reads that remain after removing duplicated reads
    :param plot_directory: Name of the directory to save plots
    :return: nothing
    """

    # Create a new figure and add a subplot
    fig = plt.figure(figsize=(18, 9))
    ax = fig.add_subplot(111)

    # x-axis is number of duplication levels we are plotting
    x_range = np.linspace(1, len(duplication_counts), len(duplication_counts))

    # Plotting percentage (out of total reads) duplicated reads
    plt.plot(x_range, duplication_counts, 'blue')

    # Plotting percentage duplicated reads (out of total unique reads)
    plt.plot(np.linspace(1, len(deduplicated_counts), len(deduplicated_counts)), deduplicated_counts, 'red')

    # Setting x ticks
    plt.xticks(x_range[::2])  # xticks positions
    ax.set_xticklabels(x_range[::2], fontsize = 12)

    # Setting x axis limits
    plt.xlim([1, x_range[-1]])

    # Setting x axis grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # Setting y axis ticks
    yticks_names = np.linspace(0, 100, 11, dtype = int)  # y axis: 11 ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 12)

    # Setting labels, title, legend
    ax.set_xlabel("Sequence duplication level (Number of occurrences)", fontsize = 16)
    ax.set_ylabel("Percentage ", fontsize = 16)
    ax.set_title("Percentage of sequences remaining after deduplication: " + str(remainder_after_deduplication),
                 fontsize = 16)

    legend_names = ["%Total sequences", "% Deduplicated sequences"]
    ax.legend(legend_names, fontsize = 16)

    # Saving plots to disk
    print "Plotting percentage of duplicated sequences"
    plot_path = os.path.join(plot_directory, "sequence_duplication.png")
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def plot_sequence_length_distribution(length_all_sequences_dict, plot_directory):
    """
    Plots the distribution of lengths of all sequences in the list of all sequences. X axis: sequence length, Y axis:
    Number of reads

    :param length_all_sequences_dict: dictionary of number of reads of a particular length, where particular length is
                                      the key
    :param plot_directory: Name of the directory to save plots
    :return: nothing
    """

    # Calculate max and min sequence length in the distribution - this will define the x-range of the plot
    max_sequence_length = max(length_all_sequences_dict.keys())
    min_sequence_length = min(length_all_sequences_dict.keys())

    # Convert to a list for matplotlib. All sequences with no occurrences are set to 0
    length_all_sequences = np.zeros(max_sequence_length+10)

    for sequence_length in length_all_sequences_dict.keys():
        length_all_sequences[sequence_length] = length_all_sequences_dict[sequence_length]

    # Create a new figure and add a subplot
    fig = plt.figure(figsize = (18, 9))
    ax = fig.add_subplot(111)

    # x-range is defined by the max and min sequence length, along with a buffer to keep the ploy away from the edges
    # of the graph
    buffer = 10
    min_x_value = min_sequence_length - buffer
    max_x_value = max_sequence_length + buffer
    x_range = np.linspace(min_x_value, max_x_value, (max_x_value-min_x_value), dtype = "int")

    # Plot a line graph of frequency of occurrence of each sequence length
    plt.plot(x_range, length_all_sequences[min_x_value:max_x_value+1])

    # Setting x lims and x ticks
    plt.xlim([min_x_value, max_x_value+1])
    plt.xticks(x_range)  # xticks positions
    ax.set_xticklabels(x_range, fontsize = 12)

    # Setting x axis grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # Setting range of y axis as ceiling(maximum length) to auto-scale
    y_range = int(math.ceil(max(length_all_sequences_dict.values()) / 100.0)) * 100

    # Setting y axis ticks
    yticks_names = np.linspace(0, y_range, 11, dtype = int)  # y axis: 10 ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 16)

    # Setting labels, title
    ax.set_xlabel("Length of sequences", fontsize = 16)
    ax.set_ylabel("Number of reads ", fontsize = 16)
    ax.set_title("Distribution of sequence lengths over all sequences", fontsize = 16)

    # Saving to disk
    print "Plotting sequence length distribution"
    plot_path = os.path.join(plot_directory, "sequence_length_distribution.png")
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def write_csv_over_represented_sequences(sorted_sequences_by_counts, total_num_sequences, plot_directory):
    """
    Writes the most over-represented sequences (most occurring sequences that make up > 0.1% of all the sequences)

    :param sorted_sequences_by_counts: dictionary of over-represented sequences (keys) with occurrences (values)
    :param total_num_sequences: Total number of sequences in input FASTQ file
    :param plot_directory: Name of the directory to save plots
    :return: nothing
    """

    print "Writing over-represented sequences to disk"

    # Create a new CSV file
    output_file = open(os.path.join(plot_directory, "over_represented_sequences.csv"), "w")
    output_file.write("Sequences, Counts, Percentage (of total)" + "\n")

    # Write sorted data to a CSV file
    for i in range(len(sorted_sequences_by_counts)):
        output_file.write(str(sorted_sequences_by_counts[i][0]) + ", " + str(sorted_sequences_by_counts[i][1]) + ", " +
                          str(100.0 * sorted_sequences_by_counts[i][1]/total_num_sequences) + "\n")
    output_file.close()
