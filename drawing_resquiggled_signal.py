#!/usr/bin/python
# coding: utf-8
import sys
import h5py
import argparse
import pandas
import numpy as np
import matplotlib.pyplot as plt

def mad(arr): 
    med = np.median(arr) 
    return np.median(np.abs(arr - med)) 

def argument_parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input', action='store',
                        help='fast5 file with resquiggled signal')
    parser.add_argument('-o', dest='output', action='store', default=None,
                        help='name of the output file (by default figure is shown, not saved)')
    parser.add_argument('-s', dest='start', type=int, default=0,
                        help='plot from this nucleotide (zero-based, inclusive)')
    parser.add_argument('-e', dest='end', type=int, default=10,
                        help='plot to this nucleotide (zero-based, inclusive)')
    parser.add_argument('--events-adress',
                        default= 'Analyses/RawGenomeCorrected_000/BaseCalled_template/Events',
                        help='where the events are stored')
    return parser.parse_args()

def main():
    arguments = argument_parsing()
    read = h5py.File(arguments.input, "r")
    start, end = arguments.start, arguments.end
    events = read[arguments.events_adress][start : end+2]
    starts = [i[2] for i in events]
    start_sygnal = read[arguments.events_adress].attrs['read_start_rel_to_raw']
    #???
    end_sygnal = starts[-1] + start_sygnal
    read_id = list(read["Raw/Reads"].items())[0][0]
    sygnal_all = read['Raw/Reads/' + read_id + '/Signal']
    sygnal_shift = np.median(sygnal_all)
    sygnal_scale = mad(sygnal_all)
    sygnal = sygnal_all[start_sygnal + starts[0] : end_sygnal]
    sygnal = np.array(sygnal, dtype='float64')
    sygnal -= sygnal_shift
    sygnal /= sygnal_scale
    starts = [i - starts[0] for i in starts]
    mids = []
    for i in range(len(starts)-1):
        mid = starts[i+1] - starts[i]
        mid /= 2
        mid += starts[i]
        mids.append(mid)
    #mids = [mids[i] + starts[i] for i in range(len(mids))]
    sekwencja = [str(i[4]) for i in events]
    sekwencja = ''.join(sekwencja)
    mean_values = [i[0] for i in events][:-1]
    mean_sygnal = sum(sygnal) / len(sygnal)
    plt.plot(sygnal)
    plt.plot(starts, [mean_sygnal] * len(starts), '|')
    plt.hlines(y=mean_values, xmin=starts[:-1], xmax=starts[1:])
    for i in range(len(mids)):
        plt.text(mids[i], mean_sygnal, sekwencja[i])
    if not arguments.output:
        plt.show()
    else:
        plt.savefig(arguments.output)

if __name__ == '__main__':
    main()

