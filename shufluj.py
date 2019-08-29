# coding: utf-8
import random
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

#sygnal = open("sygnal4").readlines()
sygnal = open(infile).readlines()
random.shuffle(sygnal)

#output = open("random_sygnal4", "w")
output = open(outfile, "w")
for line in sygnal:
    output.write(line)
    
output.close()
#get_ipython().run_line_magic('save', 'shufluj 1-7')
