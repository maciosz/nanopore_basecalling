#!/usr/bin/env python3.6
import sys
import pysam
from Bio import SeqIO

bam = pysam.AlignmentFile(sys.argv[1])
# to bardziej do testowania
#referencja = SeqIO.parse("/home/maciosz/dane_nanopor/brca_database.fa", "fasta")
# referencja powinna byc z wariantami wykrytymi w danym pacjencie
referencja = SeqIO.parse(sys.argv[2], "fasta")
brca1 = next(referencja).seq
brca2 = next(referencja).seq
referencje = {"hg38_ct_UserTrack_3545_BRCA1_plus_1000_each_side_range_chr17_43043295_43126483": brca1,
        "hg38_knownGene_uc001uub.2_BRCA2_plus_1000_each_side_range_chr13_32314474_32401266": brca2}

#counter = 0
for read in bam:
    #counter += 1
    start, end = read.reference_start, read.reference_end
    gen = read.reference_name
    if gen not in referencje.keys():
        continue
    referencja = referencje[gen]
    sekwencja = referencja[start:end]
    if read.is_reverse:
        sekwencja = sekwencja.reverse_complement()
    sekwencja = str(sekwencja)
    read_name = read.query_name
    fastq = open(read_name +".fastq", "w")
    fastq.write("@" + read_name + "\n")
    fastq.write(sekwencja + "\n")
    fastq.write("+\n")
    fastq.write("I" * len(sekwencja) + "\n") # najwyzszy score wg Phred+33
    fastq.close()
    #do testowania
    #if counter > 10:
    #    break
