#!/usr/bin/python3
import os
import argparse
import numpy as np

PENALTY1 = -0.2
PENALTY2 = -0.2
OPEN_PENALTY1 = -20.
OPEN_PENALTY2 = -20.
MATCH_REWARD = 0.1
# match reward to stala dodawana do roznicy miedzy sygnalami
# zeby punktowac dopasowywanie do siebie sygnalu, nawet jesli jakos sie roznia
# nie ma znaczenia w przypadku nieilosciowych sekwencji
# to NIE jest nagroda za dopasowanie takich samych symboli

# to powinno byc znormalizowane tak zeby uwzgledniac srednie roznice miedzy porownywanymi sygnalami
# albo moze znormalizowac najpierw sygnaly?
# w kazdym razie zeby nie bylo ze zawsze sie oplaca zrobic jeden duzy gap obejmujacy cala sekwencje

# just some flags
MOVE_IN_BOTH = 0
MOVE_IN_1 = 1
MOVE_IN_2 = 2
 
def argument_parsing():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--seq1',
                        help='first sequence')
    parser.add_argument('-2', '--seq2',
                        help='second sequence')
    parser.add_argument('-n', '--normalise', action='store_true',
                        help='should the signal be normalised?')
    parser.add_argument('-l', '--local', action='store_true',
                        help='perform local alignment (global by default)')
    parser.add_argument('-o', dest='output', action='store', default='output.tab',
                        help='name of the output file (defaults to output.tab)')
    return parser.parse_args()

def normalise(signal):
    mean_value = np.mean(signal)
    sd = np.std(signal)
    signal = [(i - mean_value) / sd for i in signal]
    return signal

def score(char1, char2, gap1_opened=False, gap2_opened=False):
    if char1 == "-":
        if gap1_opened:
            return PENALTY1
        else:
            return OPEN_PENALTY1
    elif char2 == "-":
        if gap2_opened:
            return PENALTY2
        else:
            return OPEN_PENALTY2
    try:
        # na pewno -abs()?
        #return -abs(float(char1) - float(char2))
        return -abs(char1 - char2) + MATCH_REWARD
    except TypeError:
        return int(char1 == char2)

def read_file(filename):
    values = [i.strip() for i in open(filename)]
    try:
        return [float(i) for i in values]
    except:
        return values

def read_sequence(sequence):
    if os.path.isfile(sequence):
        return read_file(sequence)
    else:
        return list(sequence)

def align(signal1, signal2):
    # skoro i tak na biezaco ustalam trace to nie potrzebuje trzymac calej tablicy
    # wystarcza mi chyba dwa rzedy
    n1, n2 = len(signal1), len(signal2)
    scores = np.zeros((n1 + 1, n2 + 1))
    traces = np.zeros((n1 + 1, n2 + 1))
    gap1_opened = False
    gap2_opened = False
    for row in range(1, n1+1):
        for column in range(1, n2+1):
            char1 = signal1[row-1]
            char2 = signal2[column-1]
            move_in_1 = scores[row - 1, column] + score(char1, "-", gap2_opened = gap2_opened)
            move_in_2 = scores[row, column - 1] + score("-", char2, gap1_opened = gap1_opened)
            move_in_both = scores[row - 1, column - 1] + score(char1, char2)
            max_score = max(move_in_1, move_in_2, move_in_both)
            if max_score == move_in_1 and max_score != move_in_both:
                gap1_opened = False
                gap2_opened = True
                traces[row, column] = MOVE_IN_1
            elif max_score == move_in_2 and max_score != move_in_both:
                gap1_opened = True
                gap2_opened = False
                traces[row, column] = MOVE_IN_2
            else:
                gap1_opened = False
                gap2_opened = False
                traces[row, column] = MOVE_IN_BOTH
            scores[row, column] = max_score
    return scores, traces
    # cos jest nie tak, obczaj tablice alajmetu dla uliniowienia abc z adef

def move(x, y, direction):
    if direction == MOVE_IN_BOTH:
        return x-1, y-1, ('.', '|', '.')
    elif direction == MOVE_IN_1:
        return x-1, y, ('.', ' ', '-')
    elif direction == MOVE_IN_2:
        return x, y-1, ('-', ' ', '.')

def get_alignment(trace, end=None):
    n1, n2 = len(trace), len(trace[0])
    if end == None:
        #global alignment
        end = -1
    direction = trace[-1][end] # czy [end][-1]?? czy end to dwie wspolrzedne?
    row, column, chars = move(n1-1, n2-1, direction)
    seq1, symbols, seq2 = chars
    not_end = True
    while not_end:
        direction = trace[row][column]
        row, column, chars = move(row, column, direction)
        char1, symbol, char2 = chars
        seq1 += char1
        seq2 += char2
        symbols += symbol
        if row == 0 or column == 0:
            not_end = False
    seq1 = list(seq1)
    seq2 = list(seq2)
    symbols = list(symbols)
    seq1.reverse()
    seq2.reverse()
    symbols.reverse()
    return seq1, symbols, seq2
        

"""
    for row in range(n1-1, 0, -1):
        for column in range(n2-1, 0, -1):
            if trace[row][column] == MOVE_IN_BOTH:
                symbols += '|'
                seq1 += '.'
                seq2 += '.'
            elif trace[row][column] == MOVE_IN_1:
                symbols += ' '
                seq1 += '.'
                seq2 += '-'
            elif trace[row][column] == MOVE_IN_2:
                symbols += ' '
                seq1 += '-'
                seq2 += '.'
    seq1.reverse()
    seq2.reverse()
    symbols.reverse()
    return seq1, symbols, seq2
"""

def pretty_print_alignment(alignment, seq1, seq2):
    #print("pretty print alignment:")
    #print(seq1)
    #print(seq2)
    #print(alignment)
    seq1_aligned, symbols, seq2_aligned = alignment
    seq1_output, symbols_output, seq2_output  = [], [], []
    for nr, symbol in enumerate(symbols):
        #print(nr, end=", ")
        if seq1_aligned[nr] == '.':
            try:
                seq1_output.append(seq1.pop(0))
            except:
                seq1_output.append("N")
        elif seq1_aligned[nr] == '-':
            seq1_output.append('-')
        if seq2_aligned[nr] == '.':
            try:
                seq2_output.append(seq2.pop(0))
            except:
                seq2_output.append("N")
        elif seq2_aligned[nr] == '-':
            seq2_output.append('-')
    #print(seq1)
    #print(seq2)
    seq1_output = ''.join([str(i) for i in seq1_output])
    seq2_output = ''.join([str(i) for i in seq2_output])
    symbols = ''.join(symbols)
    # to powinno byc bardziej 'uliniowione' wzgledem siebie
    print(seq1_output)
    print(symbols)
    print(seq2_output)

def pretty_plot_alignment(alignment, seq1, seq2):
    # ale co tu zamierzalas zaimplementowac?
    # moze plotowanie dwoch sygnalow, odpowiednio przeskalowane zeby byly zalajnowane
    # to by byl spoko ficzer
    pass

def find_optimal_end(scores):
    pass
        
def main():
    arguments = argument_parsing()
    seq1, seq2 = read_sequence(arguments.seq1), read_sequence(arguments.seq2)
    if arguments.normalise:
        seq1, seq2 = normalise(seq1), normalise(seq2)
    #print(seq1)
    #print(seq2)
    scores, traces = align(seq1, seq2)
    #print(traces)
    end = None
    if arguments.local:
        end = find_optimal_end(scores)
    alignment = get_alignment(traces, end)
    print(''.join(alignment[1]))
    #pretty_print_alignment(alignment, seq1, seq2)
    print(scores[-1][-1])

if __name__ == '__main__':
    main()
