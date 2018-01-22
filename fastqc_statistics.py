import numpy as np
import fastqc_plots


def per_sequence_quality(sequences, quality_scores, plot_directory):
    """
    Calculates average quality scores for each sequence. Average is the sum of a sequence's quality scores, divided by
    length of sequence

    :param sequences: a list of all sequences
    :param quality_scores: a list of a list of quality scores for each sequence
    :param plot_directory: Name of the directory to save plots
    :return: nothing, calls plotting function for average quality scores

    """

    # Create new list to store average quality scores
    average_q_scores = []

    for i in range(len(sequences)):

        sum_q_scores = 0    # Find sum of quality scores at all positions for this sequence.

        for j in range(len(quality_scores[i])):
            sum_q_scores += quality_scores[i][j]

        # Average score = sum divided by length of sequence
        average_q_scores.append(float(sum_q_scores/len(quality_scores[i])))

    # Calling plotting function for number of reads v/s average quality scores
    fastqc_plots.plot_per_sequence_quality(average_q_scores, plot_directory)


def per_base_sequence_content(sequences, plot_directory):

    sequence_length = len(sequences[0])
    num_sequences = len(sequences)

    percent_occurrence_all_bases = {key: np.zeros(sequence_length) for key in ["A", "C", "G", "T", "N"]}

    for each_sequence in sequences:
        for position in range(len(each_sequence)):
            base_at_position = each_sequence[position]
            # Add percent contribution of this base at this position
            percent_occurrence_all_bases[base_at_position][position] += (100.0 / float(num_sequences))

    gc_per_position = percent_occurrence_all_bases["G"] + percent_occurrence_all_bases["C"]

    # Calling plotting function for percent occurrence of all bases at each position
    fastqc_plots.plot_per_base_sequence_content(percent_occurrence_all_bases, plot_directory)

    # Calling plotting function for percent GC content at each position
    fastqc_plots.gc_content(gc_per_position, np.linspace(1, sequence_length, sequence_length, dtype = "int"),
                            plot_directory)


def per_base_sequence_quality(quality_scores, plot_directory):

    # Invert indexing of quality scores to [position][sequence]
    quality_scores_by_position = np.array(quality_scores)
    quality_scores_by_position = np.transpose(quality_scores_by_position)
    quality_scores_by_position = quality_scores_by_position.tolist()

    # Calling plotting function for sequence quality
    fastqc_plots.plot_per_base_sequence_quality(quality_scores_by_position, plot_directory)


def calculate_gc_content(sequence):
    gc_bases = ["G", "C"]
    total_gc_content = 0

    for base in sequence:
        if base in gc_bases:
            total_gc_content += 1

    # Calculate percent GC content in sequence
    percent_gc_content = 100 * float(total_gc_content)/len(sequence)

    return percent_gc_content


def per_sequence_gc_content(sequences, plot_directory):

    # TODO: Plot theoretical GC content

    gc_content_by_sequence = []

    for each_sequence in sequences:
        gc_content_by_sequence.append(calculate_gc_content(each_sequence))

    # Calling plotting function to plot Number of reads v/s mean %GC content
    fastqc_plots.plot_per_sequence_gc_content(gc_content_by_sequence, plot_directory)


def sequence_duplication(sequences, plot_directory):

    num_occurrence_each_sequence = dict()   # Dictionary to store the number of occurrences of each sequence
    total_num_sequences = len(sequences)

    for sequence in sequences:

        # If the sequence length is over 75, only the first 50 positions are used, according to the original FASTQC
        if len(sequence) > 75:
            sequence_key = sequence[0:50]
        else:
            sequence_key = sequence

        if sequence_key in num_occurrence_each_sequence.keys():
            num_occurrence_each_sequence[sequence_key] += 1
        else:
            num_occurrence_each_sequence[sequence_key] = 1

    # Define all sequences which make up more than 0.1% of all reads as over-represented (according to the original
    # FASTQC) and create a new dictionary with just these sequences and their frequencies
    over_represented_dict = dict((sequence, counts) for sequence, counts in num_occurrence_each_sequence.items() if
                                 counts >= total_num_sequences * 0.1 / 100)

    # Sort dictionary by value in descending order - most over-represented sequences appear first
    sorted_sequences_by_counts = sorted(over_represented_dict.items(), key = lambda item: item[1], reverse = True)

    # Write to CSV
    fastqc_plots.write_csv_over_represented_sequences(sorted_sequences_by_counts, total_num_sequences, plot_directory)

    # Number of duplicated sequences
    num_duplication_levels_computed = 10
    duplication_counts = np.zeros(num_duplication_levels_computed, dtype = "float")

    # Number of unique sequences
    deduplicated_counts = np.zeros(num_duplication_levels_computed, dtype = "float")

    for duplication_level in range(1, num_duplication_levels_computed+1):

        for key in num_occurrence_each_sequence.keys():

            if num_occurrence_each_sequence[key] == duplication_level:
                # Add to dictionary as a percent contribution of all (unique and duplicate) sequences
                duplication_counts[duplication_level-1] += float(duplication_level) * (100.0/total_num_sequences)

                # Add to dictionary as a percent contribution of unique sequences
                deduplicated_counts[duplication_level-1] += 100.0/len(num_occurrence_each_sequence.keys())

    remainder_after_deduplication = 100.0 * float(len(num_occurrence_each_sequence.keys()))/total_num_sequences

    # Calls plotting function for percentage duplicated reads distribution
    fastqc_plots.plot_sequence_duplication(duplication_counts, deduplicated_counts, remainder_after_deduplication,
                                           plot_directory)


def sequence_length_distribution(sequences, plot_directory):
    length_all_sequences_dict = dict()

    for sequence in sequences:
        sequence_length = len(sequence)

        if sequence_length in length_all_sequences_dict.keys():
            length_all_sequences_dict[sequence_length] += 1
        else:
            length_all_sequences_dict[sequence_length] = 1

    fastqc_plots.plot_sequence_length_distribution(length_all_sequences_dict, plot_directory)
