import os
import sys
from Bio import SeqIO
import string

 
file = sys.argv[1]
reads = open(file, 'rU')

for sequence in SeqIO.parse(reads, "fasta"):
        count = sequence.description.rsplit(' ', 1)[1].split("=")[1]
        fout = open('%s.%s' %(file,'ununiqued'), 'a')
        for num in range(1,int(count)):
        	SeqIO.write(sequence, fout, 'fasta')