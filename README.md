#fastqc_replicate
fastqc_replicate is a program that replicates the FASTQC quality control program (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). It is written in Python and generates a plots directory for all analysis plots. It takes an input FASTQ file and performs the following quality control steps (based on FASTQC):
1. Per sequence quality scores
2. Per base sequence content + Per base N content
3. Per Sequence GC content
4. Sequence length distribution
5. Duplicated sequences
6. Over represented sequences
7. Per base GC content (More on: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

There are 3 modules: fastqc_reader.py (main function), fastqc_statistics.py (calculates all the statistics) and fastqc_plots.py (plots all the statistics from fastqc_statistics.py)

To run: python fastqc_reader.py "input fastq file (.fastq or .fq)"

Example: Ran a sample FASTQ file: SRA file: SRR1972739, Zika Virus data, first 10000 lines
To fetch data: fastq-dump -X 10000 --split-files SRR1972739.sra
More info: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA257197%20
Results from fastqc_replicate:
https://docs.google.com/document/d/1kzrqJ0Bcq8eS1bqBFqxpQNUDRlJEr3GuFzZpc_EZwfk/edit?usp=sharing
https://docs.google.com/document/d/1r1c4HkbLzOodMBlnJoaS2PpxKRKt-7Nf3SCWlLQBAo4/edit?usp=sharing
