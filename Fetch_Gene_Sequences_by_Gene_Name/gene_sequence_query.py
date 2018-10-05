
import sys
from Bio import Entrez
import numpy as np

Output = open(sys.argv[2], 'w+')
# Output = open("test_out.fasta", 'w+')
Entrez.email = "sdsmith@iastate.edu"

for line in open(sys.argv[1]):
# for line in open("gene_list.txt"):
    item = line
    search_string = item+"[Gene]"

    handle = Entrez.esearch(db="nucleotide", term=search_string)
    record = Entrez.read(handle)
    ids = record['IdList']

    seq_id = ids[0] #you must implement an if to deal with <0 or >1 cases
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    record = handle.read()
    seq = record.split('\n',1)
    out = ">%s from %s \n" %(item.rstrip('\n'), seq[0][1:]) + record.split('\n',1)[1].rstrip('\n') + "\n"
    Output.write(out)
Output.close



