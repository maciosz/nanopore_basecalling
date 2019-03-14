# coding: utf-8
import sys
import h5py
import pandas
import numpy as np
import matplotlib.pyplot as plt

infile = sys.argv[1]

start_seq, end_seq = 0, 10
if len(sys.argv) == 4:
    start_seq, end_seq = [int(i) for i in sys.argv[2:]]

read = h5py.File(infile, "r")
events = read['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'][start_seq : end_seq]

# calosc
#plt.plot(read['Raw/Reads/Read_289/Signal'])
#plt.show()

# od zresquigglowanego poczatku
#plt.plot(read['Raw/Reads/Read_289/Signal'][1863:])
#plt.show()

starts = [i[2] for i in events]
start_sygnal = read['Analyses/RawGenomeCorrected_000/BaseCalled_template/Events'].attrs['read_start_rel_to_raw']

#???
end_sygnal = starts[-1] + start_sygnal
read_id = list(read["Raw/Reads"].items())[0][0]
sygnal_all = read['Raw/Reads/' + read_id + '/Signal']

def mad(arr): 
    med = np.median(arr) 
    return np.median(np.abs(arr - med)) 

sygnal_shift = np.median(sygnal_all)
sygnal_scale = mad(sygnal_all)
sygnal = sygnal_all[start_sygnal + starts[0] : end_sygnal]
sygnal = np.array(sygnal, dtype='float64')
sygnal -= sygnal_shift
sygnal /= sygnal_scale

# poczatkowy fragment
#plt.plot(sygnal)
#plt.show()

starts = [i - starts[0] for i in starts]

mids = []
for i in range(len(starts)-1):
    mid = starts[i+1] - starts[i]
    mid /= 2
    mid += starts[i]
    mids.append(mid)
    
#mids = [mids[i] + starts[i] for i in range(len(mids))]
print events
sekwencja = [str(i[4]) for i in events]
sekwencja = ''.join(sekwencja)
mean_values = [i[0] for i in events][:-1]

mean_sygnal = sum(sygnal) / len(sygnal)


# poczatkowy fragment z naniesionym resquigglowaniem
plt.plot(sygnal)
plt.plot(starts, [mean_sygnal] * len(starts), 'o')
plt.hlines(y=mean_values, xmin=starts[:-1], xmax=starts[1:])
for i in range(len(mids)):
    plt.text(mids[i], mean_sygnal, sekwencja[i])
    # mozna by dorysowac 'schodki' ze srednich w danych fragmentach
plt.show()

# po tym przykladzie widac ze cos jest nie tak:
# python drawing_resquiggled_signal.py tombo_testing/data/fast5/read1.fast5 15 30
# czemu? wyglada ok
# nie liczac tego ze ciezko znalezc jakis pattern w tym jaki poziom odpowiada jakiemu nukleotydowi...
# ale moze na poziomie calych piatek jest jakis pattern
