import argparse
import numpy as np
import pandas as pd
from scipy.stats import norm, invgauss
from Bio import SeqIO
 
"""
For a given sequence of nucleotides,
 using template describing expected signal obtained from ich 5-mer,
 generate random signal.
 In other words,
 simulate sequencing given sequence with nanopore.
"""

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        default="/home/maciosz/fast_and_spurious/eno1_transcript/eno1.fa",
                        help='input sequence in fasta file')
    parser.add_argument('-t', '--template',
                        default="/home/maciosz/fast_and_spurious/template_median68pA.model",
                        help='template')
    parser.add_argument('-l', '--length', type=int, default=20,
                        help='How many values for each nucleotide?'
                             ' If --std is set, this is the mean of the normal distr.,'
                             ' from which lengths will be drawn.'
                             ' If not, it\'s constant.')
    parser.add_argument('-s', '--std', type=int, default=0,
                        help='The standard deviation of the normal distr.,'
                             ' from which lengths will be drawn.'
                             ' Defaults to zero (no randomisation).')
    parser.add_argument('-o', '--output', action='store', default='output.txt',
                        help='name of the output file (defaults to output.txt)')
    return parser.parse_args()


def generate_values(line, length_mean, length_std):
    """
    Generate values for a given kmer,
    according to description provided in line.
    The number of values is randomised;
    it is drawn from normal distribution with mean length_mean
    and standard deviation length_std
    (length_std = 0 means no randomisation).

    line - line of a template (row from pandas.DataFrame)

    length_mean - mean of a Gauss distribution from which length of a sequence will be drawn

    length_std - standard deviation of see above.

    Description of columns in original template:

    kmer The kmer being modelled.
    level_mean The mean of a Gaussian distribution representing the current observed for this kmer.
    level_stdv The standard deviation of the above Gaussian distribution of observed currents for this kmer.
    sd_mean The mean of an inverse Gaussian distribution representing the noise observed for this kmer.
    sd_stdv The standard deviation of the above inverse Gaussian distribution of noise observed for this kmer.
    ig_lambda The lambda parameter for the above inverse Gaussian distribution of noise observed for this kmer.
    (so sd_stdv == sqrt(sd_mean ^ 3 / ig_lambda). See Wikipedia for more info.
    weight Used internally for model training purposes.
    """
    signal_var = norm(line["level_mean"], line["level_stdv"])
    mu, lmbda = line["sd_mean"], line["ig_lambda"]
    noise_var = invgauss(mu / lmbda, scale = lmbda)
    n = norm.rvs(loc=length_mean, scale=length_std)
    #n = max(1, int(n))
    n = max(0, int(n))
    values = signal_var.rvs(n) + noise_var.rvs(n)
    return values

def main():
    arguments = parse_arguments()
    template = pd.read_csv(arguments.template, delimiter="\t")
    eno2 = SeqIO.read(arguments.input, "fasta")
    kmery = []
    for i in range(len(eno2) - 5):
        kmery.append(eno2[i:i+5])
    kmery.append((eno2[-5:]))
    kmery = [str(kmer.seq) for kmer in kmery]
    signal = []
    for kmer in kmery:
        line = template.loc[template["kmer"] == kmer, :]
        values = generate_values(line, arguments.length, arguments.std)
        signal.extend(values)
    output = open(arguments.output, "w")
    for kmer in signal:
        output.write(str(kmer))
        output.write("\n")
    output.close()

"""
kmer	level_mean	level_stdv	sd_mean	sd_stdv	ig_lambda	weight
AAAAA	101.068534	2.712467	2.158905	0.844810	14.098812	18928.990707
AAAAC	97.538095	2.712467	2.199744	0.868895	14.098790	17729.724615
AAAAG	98.617195	2.712467	2.623405	1.131637	14.098787	25036.063106
AAAAT	96.610107	2.712467	1.928888	0.713460	14.098788	8936.195540
"""

if __name__ == '__main__':
    main()
