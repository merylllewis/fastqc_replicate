#/usr/binb/env python3
"""
fastqc_reader.py replicates the FASTQC tool and calculates various statistics for quality assurance of FASTQ files
FASTQ files are to be parsed as command line arguments
"""
__author__ = "Meryl Lewis"
__email__ = "merylllewis@gmail.com"

import sys
import argparse
import os
import fastqc_statistics as fastqc_stats


def fastq_reader(file_name):
    """
    Reads an input FASTQ file and returns a list of all the reads and their quality (phred) scores
    :param file_name: input fastq file: line 1: header/identifier, line 2: sequence, line 3: "+", line 4: quality scores
                      and so on...
    :return: sequences, quality scores
    """
    fastq_read = open(file_name, 'r')
    header = fastq_read.readline().strip("\r\n")    # reads first line/ header of the file
    sequences = []
    quality_scores = []

    while header:
        sequences.append(fastq_read.readline().strip("\r\n"))
        fastq_read.readline().strip("\r\n")     # "+" line
        quality_score_characters = fastq_read.readline().strip("\r\n")
        sequence_quality_scores = [int(ord(score) - 33) for score in quality_score_characters]
        # 33 is the encoding conversion for new FASTQ files
        quality_scores.append(sequence_quality_scores)

        header = fastq_read.readline().strip("\r\n")

    return sequences, quality_scores


def main(arguments):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fastq_file", type = str, action = "store", help = "input fastq file name")
    args = parser.parse_args(arguments[1:])
    sequences, quality_scores = fastq_reader(args.input_fastq_file)

    plot_directory = os.path.join(os.getcwd(), args.input_fastq_file.split(".")[0] + "_plots")
    if not os.path.exists(plot_directory):
        print "Making new directory for plots"
        os.makedirs(plot_directory)
    fastqc_stats.per_sequence_quality(sequences, quality_scores, plot_directory)
    fastqc_stats.per_base_sequence_content(sequences, plot_directory)
    fastqc_stats.per_base_sequence_quality(quality_scores, plot_directory)


if __name__ == "__main__":
    main(sys.argv)