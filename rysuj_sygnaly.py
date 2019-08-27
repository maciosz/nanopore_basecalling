#!/usr/bin/python3
# coding: utf-8
import os
import collections
import h5py
import matplotlib.pyplot as plt

sciezka = "/home/maciosz/dane_nanopor/genXone_fast5/20180705_1547_Sample_genxone_poznan_brca_s181,182,185-202_brca_kp175_MN106_SLSK108/20180705_1547_brca_s181,182,185-202_brca_kp175_MN106_SLSK108/reads/0/"
pliki_fast5 = os.listdir(sciezka)
processed_reads = collections.defaultdict(int)
#print(len(pliki_fast5))

for plik in pliki_fast5[:100]:
    read = h5py.File(sciezka + plik, "r")
    plik = plik.split('_')
    nr_readu = plik[plik.index("read") + 1]
    processed_reads[nr_readu] += 1
    try:
        signal = read['Raw']['Reads']["Read_" + nr_readu]['Signal']
    except:
        print(plik)
        read.visit(print)
        continue
    plt.plot(signal)
    nazwa = nr_readu + "_signal_" + str(processed_reads[nr_readu])
    plt.savefig(nazwa + ".png")
    plt.clf()
    plt.plot(signal[:int(len(signal)/10)])
    plt.savefig(nazwa + "_cut.png")
    plt.clf()


"""
events = read['Analyses']['RawGenomeCorrected_000']['BaseCalled_template']['Events']
sekwencja = ''.join([chr(i[4][0]) for i in events])
"""
