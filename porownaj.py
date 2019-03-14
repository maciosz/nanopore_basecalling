#!/mnt/space/genxone/pysam_biopython_p3_venv/bin/python3 
import os
import sys
import copy
import re
import collections
import argparse
import numpy as np
from Bio import Align
from blast2sam import cigar

aligner = Align.PairwiseAligner()
aligner.open_gap_score = -1.0
aligner.extend_gap_score = 0.1
OPERATIONS = {'insertions':'I',
              'deletions':'D',
              'mismatches':'X',
              'matches':'='}
OPERATIONS_REVERSE = {}
for key, value in OPERATIONS.items():
    OPERATIONS_REVERSE[value] = key

SKIPPED = 0
SKIPPED_TRAILING = 0

def parse_cigar_string(cigar_string):
    global SKIPPED
    global SKIPPED_TRAILING
    operations = re.findall("[0-9]+[IDXM=]", cigar_string)
    operations = [(int(i[:-1]), i[-1]) for i in operations]
    #operations = [(int(i[:-1]), i[-1]) for i in operations if i[:-1] != "0"]
    # paczac na kod po prostu on zaczyna od operacji =
    # i jesli dalej ma inna to zapisuje 0=. Glupie ale malo problematyczne.
    while operations[0][0] == 0 or operations[0][1] in ('I', 'D'):
        operations = operations[1:]
        SKIPPED += 1
    while operations[-1][1] in ('I', 'D'):
        operations = operations[:-1]
        SKIPPED_TRAILING += 1
    return operations

def summarise_operations(operations):
    symbols = OPERATIONS.values()
    operations_by_symbols = {}
    summary = {}
    alignment_length = float(sum([i[0] for i in operations]))
    for symbol in symbols:
        operations_by_symbols[symbol] = [i[0] for i in operations if i[1] == symbol]
    for symbol, counts in operations_by_symbols.items():
        current_summary = {}
        current_summary['number'] = len(counts) 
        current_summary['summaric_length'] = sum(counts)
        current_summary['percent_of_length'] = sum(counts) / alignment_length * 100
        if len(counts) != 0:
            current_summary['mean_length'] = sum(counts) / len(counts)
            current_summary['min'] = min(counts)
            current_summary['max'] = max(counts)
            current_summary['median'] = np.median(counts)
        summary[symbol] = current_summary.copy()
    return summary

def get_best_alignment(read1, read2, cutoff = 50):
    alignments = aligner.align(read1, read2)
    score = alignments.score
    min_indels = len(read1)
    best_alignment = None
    #print("%d alignments" % len(alignments))
    counter = 0
    for alignment in alignments:
        splitted = alignment.__str__().strip().split('\n')
        query = splitted[0]
        subject = splitted[-1]
        query_clipped = query.strip('-')
        subject_clipped = subject.strip('-')
        indels = len(re.findall("-+", query_clipped + " " + subject_clipped))
        if indels <= min_indels:
            min_indels = indels
            best_alignment = (query, subject)
        counter += 1
        if counter >= cutoff:
            break
    #print(best_alignment)
    return best_alignment, score

def save_indel_sequences(alignment, output_prefix, context, read_name):
    query, subject = alignment
    query_deletions = open(output_prefix + "_basecalled_around_deletions.fa", "a+")
    subject_deletions = open(output_prefix + "_real_on_deletions.fa", "a+")
    query_insertions = open(output_prefix + "_basecalled_on_insertions.fa", "a+")
    subject_insertions = open(output_prefix + "_real_around_insertions.fa", "a+")
    save_deletion_sequences((query, subject),
                            (query_deletions, subject_deletions),
                            context, read_name, "deletion")
    save_deletion_sequences((subject, query),
                            (subject_insertions, query_insertions),
                            context, read_name, "insertion")

def save_deletion_sequences(reads, outputs, context, read_name, indel_name, threshold = 0):
    query, subject = reads
    query_fragments, subject_fragments = outputs
    deletions = re.finditer("[ACTG](-+)[ACTG]", query)
    nr = 0
    for deletion in deletions:
        nr += 1
        start, end = deletion.span()
        start += 1
        end -= 1
        if start < threshold or (len(query) - end) < threshold:
            continue
        left_context = max(0, start - context)
        right_context = min(len(query), end + context)
        sep = "-N-"
        query_fragment = query[left_context:start] + sep + query[end:right_context]
        subject_fragment = subject[left_context:start] + sep + \
                           subject[start:end] + sep + \
                           subject[end:right_context]
        query_fragments.write(">%s%d_%s length:%d start-end:%d-%d start-end_with_context:%d-%d\n"
                              % (indel_name, nr, read_name, end-start, start, end, left_context, right_context))
        query_fragments.write("%s\n" % query_fragment)
        subject_fragments.write(">%s%d_%s length:%d start_end:%d-%d start_end_with_context:%d-%d\n"
                              % (indel_name, nr, read_name, end-start, start, end, left_context, right_context))
        subject_fragments.write("%s\n" % subject_fragment)

def compare_reads(alignment):
    query, subject = alignment
    #query, subject = get_best_alignment(read1, read2)
    #query, _, subject = alignment.__str__().split()
    try:
        cigar_string = cigar(subject, query, 0, len(query), len(query))
        #print(cigar_string)
    except IndexError:
        print("Index Error while constructing cigar:")
        print(subject)
        print(query)
        sys.exit()
    operations = parse_cigar_string(cigar_string)
    summary = summarise_operations(operations)
    return summary

def open_read(filename):
    return open(filename).readlines()[1].strip()

def summarise_comparisons(comparisons):
    operation_summary_template = {'summaric_number': 0,
                                  'mean_lengths': [],
                                  'summaric_length': 0,
                                  'minimum_length': 1000,
                                  'maximum_length': 0,
                                  'percents': []}
    summary = {}
    scores = [i['score'] for i in comparisons]
    summary['scores'] = {'mean': np.mean(scores),
                         'median': np.median(scores),
                         'min': min(scores),
                         'max': max(scores)}
    for key in OPERATIONS.keys():
        summary[key] = copy.deepcopy(operation_summary_template)
    for operation_name, operation_symbol in OPERATIONS.items():
        for comparison in comparisons:
            operation_summary = comparison[operation_symbol]
            #operation_name = OPERATIONS_REVERSE[operation_symbol]
            summary[operation_name]['summaric_number'] += operation_summary['number']
            summary[operation_name]['summaric_length'] += operation_summary['summaric_length']
            summary[operation_name]['percents'].append(operation_summary['percent_of_length'])
            if operation_summary['number'] != 0:
                summary[operation_name]['mean_lengths'].append(operation_summary['mean_length'])
                minimum = operation_summary['min']
                maximum = operation_summary['max']
                if minimum < summary[operation_name]['minimum_length']:
                    summary[operation_name]['minimum_length'] = minimum
                if maximum > summary[operation_name]['maximum_length']:
                    summary[operation_name]['maximum_length'] = maximum
        summary[operation_name]['mean_length'] = summary[operation_name]['summaric_length'] / summary[operation_name]['summaric_number']
        summary[operation_name]['mean_of_mean_lengths'] = sum(summary[operation_name]['mean_lengths']) / len(comparisons)
        _ = summary[operation_name].pop('mean_lengths', None)
        summary[operation_name]['mean_number'] = summary[operation_name]['summaric_number'] / len(comparisons)
        summary[operation_name]['mean_percent'] = np.mean(summary[operation_name]['percents'])
        summary[operation_name]['median_percent'] = np.median(summary[operation_name]['percents'])
        summary[operation_name]['min_percent'] = min(summary[operation_name]['percents'])
        summary[operation_name]['max_percent'] = max(summary[operation_name]['percents'])
        _ = summary[operation_name].pop('percents', None)
    return summary

def print_summary(global_summary):
    statistics = ["summaric_number", "mean_number", "summaric_length", "mean_length",
                  "mean_of_mean_lengths", "minimum_length", "maximum_length"]
    print("\t".join(["operation"] + statistics))
    for operation, summary in global_summary.items():
        print("\t".join([operation] +
                        ["{:,.2f}" % summary[i] if type(summary[i]) == float else "{:,}".format(summary[i]) for i in statistics]))

def save_summary(global_summary, output):
    output.write("#READS_SUMMARY\n")
    output.write("\t".join(["number",
                            "basecalled_summaric_length", "basecalled_mean_length",
                            "true_summaric_length", "true_mean_length"]))
    output.write("\n")
    output.write("{:,}\t{:,}\t{:,.2f}\t{:,}\t{:,.2f}\n".format(global_summary['reads']['number'],
                                             global_summary['reads']['basecalled_summaric_length'],
                                             global_summary['reads']['basecalled_mean_length'],
                                             global_summary['reads']['true_summaric_length'],
                                             global_summary['reads']['true_mean_length']))
    output.write("#\n#ALIGNMENT_SCORES_SUMMARY\n")
    statistics = ["mean", "median", "min", "max"]
    output.write("\t".join(statistics))
    output.write("\n")
    output.write("\t".join(["{:,.2f}".format(global_summary['scores'][i]) for i in statistics]))
    output.write("\n")
    output.write("#\n#BASECALLING_SUMMARY\n")
    statistics = ["summaric_number", "mean_number", "summaric_length", "mean_length",
                  "mean_of_mean_lengths", "minimum_length", "maximum_length",
                  "mean_percent", "median_percent", "min_percent", "max_percent"]
    output.write("\t".join(["operation"] + statistics))
    output.write("\n")
    for operation in ("matches", "mismatches", "insertions", "deletions"):
        summary = global_summary[operation]
        output.write("\t".join([operation] + ["{:,}".format(summary[i]) if type(summary[i]) == int else "{:,.2f}".format(summary[i]) for i in statistics]))
        output.write("\n")

def analyse_basecalling(basecalling_directory, true_directory,
                        output_prefix, context):
    comparisons = []
    basecalled_reads = os.listdir(basecalling_directory + "/result/")
    sys.stderr.write("%d reads to parse\n" % len(basecalled_reads))
    counter = 0
    read_summary = collections.defaultdict(int)
    #summary = collections.defaultdict(dict)
    for read_name in basecalled_reads:
        #print(read_name)
        try:
            true_read = open_read(true_directory + read_name)
            basecalled_read = open_read(basecalling_directory + "/result/" + read_name)
        except FileNotFoundError:
            continue
        #print("comparing...")
        read_summary['number'] += 1
        read_summary['basecalled_summaric_length'] += len(basecalled_read)
        read_summary['true_summaric_length'] += len(true_read)
        best_alignment, score = get_best_alignment(basecalled_read, true_read)
        save_indel_sequences(best_alignment, output_prefix, context, read_name)
        comparison = compare_reads(best_alignment)
        comparison['score'] = score
        comparisons.append(comparison)
        #print(comparison)
        counter += 1
        if counter % 500 == 0:
            sys.stderr.write("%d reads parsed\n" % counter)
        # DO TESTOW:
        #if counter > 50:
        #    break
    summary = summarise_comparisons(comparisons)
    read_summary['basecalled_mean_length'] = read_summary['basecalled_summaric_length'] / read_summary['number']
    read_summary['true_mean_length'] = read_summary['true_summaric_length'] / read_summary['number']
    summary['reads'] = read_summary
    return summary

def argument_parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', dest='basecalling', action='store',
                        help='basecalling directory (with result directory in it)')
    parser.add_argument('-t', dest='true', action='store',
                        help='directory with true sequences of reads,'
                             ' each read in different .fastq file')
    parser.add_argument('-o', dest='output', action='store', default='output.tab',
                        help='name of the output file (defaults to output.tab)')
    parser.add_argument('-p', dest='prefix', action='store', default='output',
                        help='prefix for the output files (defaults to output)')
    parser.add_argument('-c', dest='context', action='store', default=10, type=int,
                        help='context of indels - how many nucleotides around an indel to save')
    return parser.parse_args()

def main():
    arguments = argument_parsing()
    basecalling_directory = arguments.basecalling
    true_directory = arguments.true
    output = open(arguments.output, 'w')
    summary = analyse_basecalling(basecalling_directory, true_directory,
                                  arguments.prefix, arguments.context)
    #print_summary(summary)
    save_summary(summary, output)
    output.close()
    sys.stderr.write("%d skipped\n" % SKIPPED)
    sys.stderr.write("%d trailing skipped\n" % SKIPPED_TRAILING)
    sys.stderr.write("Summary written to %s, enjoy.\n" % arguments.output)

if __name__=='__main__':
    main()
