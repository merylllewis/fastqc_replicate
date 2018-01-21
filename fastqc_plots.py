import os
import numpy as np
import matplotlib.pyplot as plt
import math


def plot_per_sequence_quality(average_q_scores, file_name):
    """
    Plots average quality score distribution. X axis: Number of reads, Y axis: average quality scores
    :param average_q_scores: same as above
    :param file_name: same as above
    :return:
    """
    bins_array = np.linspace(0, 41, 42)
    fig = plt.figure(figsize = (18, 9))  # width, height
    ax = fig.add_subplot(111)
    bin_frequencies, bin_edges = np.histogram(average_q_scores, bins = bins_array)

    plt.plot(bins_array[:-1], bin_frequencies)

    # setting horizontal grid
    ax.xaxis.grid(which = "major", linestyle = '-.', linewidth = 0.5, alpha = 0.5)

    # setting x axis ticks
    xtick_positions = bins_array[::2]
    plt.xticks(xtick_positions)  # xticks positions
    xticks_names = np.array(bins_array[::2], dtype = "int")
    ax.set_xticklabels(xticks_names, fontsize = 16)  # vertical xticks

    # Setting x axis limits
    plt.xlim([0, 42])

    # setting range of y axis as ceiling(max bin frequency), auto-scales the y axis into 10 divisions
    y_range = int(math.ceil(np.amax(bin_frequencies) / 100)) * 100

    yticks_names = np.linspace(0, y_range, 11, dtype = int)  # y axis: 10 ticks
    # setting y axis ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 16)

    ax.set_xlabel("Average quality scores)", fontsize = 16)
    ax.set_ylabel("Number of reads ", fontsize = 16)
    ax.set_title("Quality score distribution over all sequences ", fontsize = 16)
    legend_name = os.path.basename(file_name)[:-len("_plots")]

    #TODO: figure out why legend name is not showing up
    ax.legend(str(legend_name), loc = "upper left", fontsize = 16)

    plot_path = os.path.join(file_name, "per_sequence_quality.png")
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def plot_per_base_sequence_content(group_by_base, plot_directory):
    color_palette = ["blue", "black", "red", "green", "violet"]
    number_positions = len(group_by_base[0])

    fig = plt.figure(figsize = (18, 9))  # width, height
    ax = fig.add_subplot(111)
    x_range = np.linspace(1, number_positions, number_positions, dtype = "int")
    plt.plot(x_range, group_by_base[0], color = color_palette[0])
    plt.plot(x_range, group_by_base[1], color = color_palette[1])
    plt.plot(x_range, group_by_base[2], color = color_palette[2])
    plt.plot(x_range, group_by_base[3], color = color_palette[3])
    plt.plot(x_range, group_by_base[4], color = color_palette[4])

    # setting horizontal grid
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

    ax.set_xlabel("Position in a read)", fontsize = 16)
    ax.set_ylabel("Percentage occurrence ", fontsize = 16)
    ax.set_title("Sequence content across all bases ", fontsize = 16)
    legend_name = ["A", "C", "G", "T", "N"]

    ax.legend(legend_name, loc = "upper left", fontsize = 16)

    plot_path = os.path.join(plot_directory, "per_base_sequence_content.png")
    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()



def gc_content(cg_bases_by_position, x_range, plot_directory):

    fig = plt.figure(figsize = (18,9))
    ax = fig.add_subplot(111)

    plt.plot(x_range, cg_bases_by_position)
    plot_path = os.path.join(plot_directory, "gc_content.png")
    plt.xticks(x_range[::5])  # xticks positions
    ax.set_xticklabels(x_range[::5], fontsize = 16)

    # Setting x axis limits
    plt.xlim([1, x_range[-1]])

    # setting range of y axis as 100 %

    yticks_names = np.linspace(0, 100, 11, dtype = int)  # y axis: 10 ticks
    # setting y axis ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 16)

    ax.set_xlabel("Position in a read)", fontsize = 16)
    ax.set_ylabel("Percentage occurrence ", fontsize = 16)
    ax.set_title("Sequence content across all bases ", fontsize = 16)
    legend_name = ["GC_content"]
    ax.legend(legend_name, loc = "upper left", fontsize = 16)

    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()


def plot_per_base_sequence_quality(quality_by_position, plot_directory):
    fig = plt.figure(figsize = (18, 9))
    ax = fig.add_subplot(111)
    x_range = np.linspace(1, len(quality_by_position), len(quality_by_position), dtype = "int")

    plt.boxplot(quality_by_position, showfliers = False)
    plot_path = os.path.join(plot_directory, "per_base_sequence_quality.png")
    plt.xticks(x_range[::2])  # xticks positions
    ax.set_xticklabels(x_range[::2], fontsize = 12)

    # Setting x axis limits
    plt.xlim([1, x_range[-1]])

    # setting range of y axis as 100 %

    yticks_names = np.linspace(0, 40, 41, dtype = int)  # y axis: 10 ticks
    # setting y axis ticks
    plt.yticks(yticks_names)
    ax.set_yticklabels(yticks_names, fontsize = 12)

    ax.set_xlabel("Position in a read (bp)", fontsize = 16)
    ax.set_ylabel("Quality scores ", fontsize = 16)
    ax.set_title("Quality scores across all bases ", fontsize = 16)

    plt.savefig(plot_path, bbox_inches = "tight")
    plt.close()