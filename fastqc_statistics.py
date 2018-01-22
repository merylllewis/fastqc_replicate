import numpy as np
import fastqc_plots


def per_sequence_quality(sequences, quality_scores, plot_directory):
    """
    Calculates average quality scores for each sequence. Average is the sum of a sequence's quality scores/
    length of sequence
    :param sequences: a list of all sequences
    :param quality_scores: a list of a list of quality scores for each sequence
    :param plot_directory: Name of the directory to save plots
    :return: nothing, calls plotting function for average quality scores
    """
    average_q_scores = []
    # y axis: number of reads, x axis: quality scores
    for i in range(len(sequences)):
        sum_q_scores = 0
        for j in range(len(quality_scores[i])):
            sum_q_scores += quality_scores[i][j]
        average_q_scores.append(float(sum_q_scores/len(quality_scores[i])))

    fastqc_plots.plot_per_sequence_quality(average_q_scores, plot_directory)


def per_base_sequence_content(sequences, plot_directory):
    occurrence_bases = {key: np.zeros(len(sequences[0])) for key in ["A", "C", "G", "T", "N"]}
    for each_sequence in sequences:

        for position in range(len(each_sequence)):
                occurrence_bases[each_sequence[position]][position] += 1
    group_by_base = []
    for base in occurrence_bases.keys():
        group_by_base.append(occurrence_bases[base])

    number_positions = len(occurrence_bases["A"])

    for position in range(0, number_positions):

        cumulative = 0

        for base in occurrence_bases.keys():
            cumulative += occurrence_bases[base][position]

        for base in range(0, len(occurrence_bases)):
            group_by_base[base][position] /= float(cumulative) / 100

    fastqc_plots.plot_per_base_sequence_content(group_by_base, plot_directory)

    gc_per_position = np.zeros(number_positions)
    for i in range(number_positions):
        gc_per_position[i] = group_by_base[1][i] + group_by_base[2][i]

    fastqc_plots.gc_content(gc_per_position, np.linspace(1, number_positions, number_positions, dtype = "int"), plot_directory)


def per_base_sequence_quality(quality_scores, plot_directory):

    quality_by_position = []

    for i in range(0, len(quality_scores[0])):
        new_list = []

        for j in range(0, len(quality_scores)):
            new_list.append(quality_scores[j][i])

        quality_by_position.append(new_list)

    fastqc_plots.plot_per_base_sequence_quality(quality_by_position, plot_directory)


def gc_content(sequence):
    gc_bases = ["G", "C"]
    total_gc_content = 0

    for base in sequence:
        if base in gc_bases:
            total_gc_content += 1
    total_gc_content = 100 * float(total_gc_content)/len(sequence)

    return total_gc_content


def per_sequence_gc_content(sequences, plot_directory):
    # TODO: figure out theoretical gc content
    gc_content_by_sequence = []
    for each_sequence in sequences:
        gc_content_by_sequence.append(gc_content(each_sequence))

    fastqc_plots.plot_per_sequence_gc_content(gc_content_by_sequence, plot_directory)


def sequence_duplication(sequences, plot_directory):
    counts_each_sequence = dict()
    total_num_sequences = len(sequences)
    for sequence in sequences:

        if len(sequence) > 75:
            sequence_key = sequence[0:50]
        else:
            sequence_key = sequence

        if sequence_key in counts_each_sequence.keys():
            counts_each_sequence[sequence_key] += 1
        else:
            counts_each_sequence[sequence_key] = 1

    duplication_level_counts = np.zeros(10, dtype = "float")
    deduplicated_counts = np.zeros(10, dtype = "float")

    for duplication_level in range(1, 11):
        for key in counts_each_sequence.keys():
            if counts_each_sequence[key] == duplication_level:
                duplication_level_counts[duplication_level-1] += float(duplication_level)*(100.0/total_num_sequences)
                deduplicated_counts[duplication_level-1] += 100.0/len(counts_each_sequence.keys())

    remainder_after_deduplication = 100.0 * float(len(counts_each_sequence.keys()))/total_num_sequences

    fastqc_plots.plot_sequence_duplication(duplication_level_counts, deduplicated_counts, remainder_after_deduplication,
                                           plot_directory)



