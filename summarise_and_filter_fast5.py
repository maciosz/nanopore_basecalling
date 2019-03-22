#!/usr/bin/env python3
import os
import sys
import argparse
import shutil
import numpy as np
import warnings
warnings.filterwarnings(action="ignore", category=FutureWarning)
import h5py

FLAGS = None

def get_infiles():
    if FLAGS.directory:
        files = os.listdir(FLAGS.directory)
        files = [FLAGS.directory + "/" + name for name in files if name.endswith(".fast5")]
    elif FLAGS.fast5:
        files = FLAGS.fast5
    return files

def save_read(filename):
    #name = os.path.basename(read.filename)
    name = os.path.basename(filename)
    if FLAGS.suffix:
        basename = ''.join(name.split('.')[:-1])
        name = basename + FLAGS.suffix + ".fast5"
    name = FLAGS.output_dir + name
    if os.path.isfile(name):
        sys.exit("File %s exists" % name)
    #output = h5py.File(name, "w")
    shutil.copy(filename, name)

def get_signal_address(read):
    reads = read.get("Raw/Reads")
    if reads:
        read_id = list(reads.keys())[0]
    else:
        sys.exit("Raw/Reads group not present. I don't know what to do. Help.")
    signal_address = "Raw/Reads/" + read_id + "/Signal"
    if not read.get(signal_address):
        sys.exit("%s group not present. You might want to edit get_signal_address function. Or check your file."
                 % signal_address)
    return signal_address

def get_length(read):
    signal_address = get_signal_address(read)
    return len(read[signal_address])

def print_length_summary(lengths):
    if len(lengths) == 0:
        sys.exit("No lengths found.")
    print("## Length summary")
    #print(f"Mean: {np.mean(lengths):,.2f}")
    print("Mean: {:,.2f}".format(np.mean(lengths)))
    print("Quintiles: {:,.2f}\t{:,.2f}\t{:,.2f}".format(
          *[np.quantile(lengths, i) for i in (0.25, 0.5, 0.75)]))
    print(("Deciles:" + "\t{:,.2f}"*9).format(
          *[np.percentile(lengths, i) for i in range(10, 91, 10)]))
    print("Min, max: {:,}\t{:,}".format(min(lengths), max(lengths)))

def set_action(args):
    actions = {'filter': ['filter', 'f'],
               'summarise': ['summarise', 'summarize', 'summary', 's']}
    actions_rev = {}
    for action, abbrevs in actions.items():
        for abbrev in abbrevs:
            actions_rev[abbrev] = action
    #actions = {'filter': 'filter',
    #           'f': 'filter',
    #           'summarise': 'summarise',
    #           'summarize': 'summarise',
    #           's' : 'summarise'}
    #args.action = actions.get(args.action)
    args.action = actions_rev.get(args.action)

def argument_parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='action', default='filter',
                        help='action to perform:'
                             ' filter (f) / summarise (s)')
    parser.add_argument('-d', dest='directory', default=None,
                        help='directory with fast5 files')
    parser.add_argument('-f', dest='fast5', nargs='+', default=None,
                        help='fast5 files')
    parser.add_argument('-o', dest='output_dir', default='.',
                        help='directory to save output files;'
                        ' defaults to current directory')
    parser.add_argument('-s', dest='suffix', default=None,
                        help='suffix for the output files (defaults to None)')
    parser.add_argument('-t', dest='threshold', type=int, default=0,
                        help='save reads with length equal or above this threshold.'
                             ' Defaults to zero (no filtering, just copy everything).')
    args = parser.parse_args()
    if not bool(args.directory) ^ bool(args.fast5):
        sys.exit("Provide either directory or fast5 file(s). Exactly one of them. Not both. Not none. One.")
    args.output_dir += "/"
    set_action(args)
    if not args.action:
        sys.exit("Provide a valid action description: filter / summarise / summary / f / s."
                 " summarize is also acceptable but I will judge you silently.")
    return args

def main():
    global FLAGS
    FLAGS = argument_parsing()
    files = get_infiles()
    lengths = []
    for filename in files:
        fast5 = h5py.File(filename, "r")
        length = get_length(fast5)
        if FLAGS.action == "filter":
            if length >= FLAGS.threshold:
                save_read(filename)
        elif FLAGS.action == "summarise":
            lengths.append(length)
    if FLAGS.action == "summarise":
       print_length_summary(lengths) 

if __name__=='__main__':
    main()
