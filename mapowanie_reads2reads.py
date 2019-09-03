import sys
import re
import numpy as np

# operations defined in CIGAR strings
# 0 M match/mismatch
# 1 I insertion
# 2 D deletion
# 3 N intron
# 4 S soft clip
# 5 H hard clip
# 6 P padding (?)
# 7 = match
# 8 X mismatch
OPERATIONS_QUERY_CONSUMING = set(('M', 'I', 'S', '=', 'X'))
OPERATIONS_REFERENCE_CONSUMING = set(('M', 'D', 'N', '=', 'X'))
OPERATIONS_INSERTION_LIKE = OPERATIONS_QUERY_CONSUMING - OPERATIONS_REFERENCE_CONSUMING
OPERATIONS_DELETION_LIKE = OPERATIONS_REFERENCE_CONSUMING - OPERATIONS_QUERY_CONSUMING

sam_file = open(sys.argv[1])
all_reads = set()
unmapped = set()
mismapped = set()
mapped_correctly = set()
mapped_perfectly = set()
deletions = []
number_of_deletions = []
insertions = []
number_of_insertions = []
matched = []
for line in sam_file:
    if line.startswith("@"):
        continue
    line = line.strip().split()
    query, target, cigar = line[0], line[2], line[5]
    all_reads.add(query)
    if target == "*":
        unmapped.add(query)
    elif query != target:
        mismapped.add(query)
    else:
        mapped_correctly.add(query)
        if len(re.findall('[0-9][M=]', cigar)) == 1:
            mapped_perfectly.add(query)
        else:
            query_insertions = re.findall('([0-9]+)I', cigar)
            query_insertions = [int(i) for i in query_insertions]
            query_deletions = re.findall('([0-9]+)D', cigar)
            query_deletions = [int(i) for i in query_deletions]
            query_matches = re.findall('([0-9]+)[XM=]', cigar)
            deletions.extend(query_deletions)
            number_of_deletions.append(len(query_deletions))
            insertions.extend(query_insertions)
            number_of_insertions.append(len(query_insertions))
            matched.append(sum(int(i) for i in query_matches))

print "All reads: %d" % len(all_reads)
print "Unmapped: %d" % len(unmapped)
print "Mismapped: %d" % len(mismapped)
print "Mapped correctly: %d" % len(mapped_correctly)
print "\tIncluding mapped perfectly (without indels): %d" % len(mapped_perfectly)

print "\tMean insertion length: %f" % np.mean(insertions)
print "\tMedian insertion length: %d" % np.median(insertions)
print "\tMax insertion length: %d" % np.max(insertions)
print "\tMean and median number of insertions per read: %f %f" % (np.mean(number_of_insertions), np.median(number_of_insertions))
print "\tMean deletion length: %f" % np.mean(deletions)
print "\tMedian deletion length: %d" % np.median(deletions)
print "\tMax deletion length: %d" % np.max(deletions)
print "\tMean and median number of deletions per read: %f %f" % (np.mean(number_of_deletions), np.median(number_of_deletions))
print "\tMean and median number of matched or mismatched bases per read: %f %f" % (np.mean(matched), np.median(matched))

def get_list_of_operations(cigar, operation='M'):
    lengths = re.findall('([0-9]+)' + operation, cigar)
    return [int(i) for i in lengths]

def split_cigar(cigar):
    cigar = re.sub('[A-Z][0-9]', lambda x: x.group()[0] + " " + x.group()[1:], cigar)
    splitted = cigar.split()
    return [(int(i[:-1]), i[-1]) for i in splitted]

def get_region_coords(splitted_cigar, which):
    start, end = 0, 0
    position_on_query = 1
    deletions = 0
    insertions = 0
    counter = 0
    for how_many, what in splitted_cigar:
        if counter >= which:
            start = position_on_query
            end = start + how_many
            break
        coutner += 1
        if what in OPERATIONS_QUERY_CONSUMING:
            position_on_query += how_many
        #if what in OPERATIONS_REFERENCE_CONSUMING:
        #    position_on_target += how_many
        if what in OPERATIONS_DELETION_LIKE:
            deletions += how_many
        if what in OPERATIONS_INSERTION_LIKE:
            insertions += how_many
    return start, end
    #return coordinate + start + deletions - insertions - 1


sam_file = open(sys.argv[1])
for line in sam_file:
    if line.startswith("@"):
        continue
    line = line.strip().split()
    query, target, cigar = line[0], line[2], line[5]
    if target == "*":
        continue
    query_insertions = get_list_of_operations(cigar, 'I')
    query_deletions = get_list_of_operations(cigar, 'D')
    query_any_matches = get_list_of_operations(cigar, '[M=X]')
    query_exact_matches = get_list_of_operations(cigar, '=')
    query_mismatches = get_list_of_operations(cigar, 'X')
    cigar = split_cigar(cigar)
    position_on_query = 1
    insertion_coords = []
    insertions = 0
    deletions = 0
    for how_many, what in cigar:
        if what == 'I':
            start = position_on_query
            end = start + how_many + 1
            insertion_coords.append((start, end))
        if what in OPERATIONS_QUERY_CONSUMING:
            position_on_query += how_many
        #if what in OPERATIONS_REFERENCE_CONSUMING:
        #    position_on_target += how_many
        if what in OPERATIONS_DELETION_LIKE:
            deletions += how_many
        if what in OPERATIONS_INSERTION_LIKE:
            insertions += how_many
    insertion_sequences = []
    sequences_before_insertion = []
    sequences_after_insertion = []
    read_sequence = line[9]
    for start, end in insertion_coords:
        insertion_sequences.append(read_sequence[start:end])
        #if insertion_sequences[-1] == '':
        #    print start, end
        insertion_length = end - start
        start_before_insertion = max(1, start - insertion_length*2)
        sequences_before_insertion.append(read_sequence[start_before_insertion : start])
        end_after_insertion = min(len(read_sequence), end + insertion_length*2)
        sequences_after_insertion.append(read_sequence[end : end_after_insertion])

"""
for i in xrange(len(insertion_sequences)):
    if insertion_sequences[i] in sequences_before_insertion[i] or insertion_sequences[i] in sequences_after_insertion[i]:
        print "Yup:", i
        print insertion_sequences[i], sequences_after_insertion[i]
    else:
        print "Nope:", i, len(insertion_sequences[i])
"""


