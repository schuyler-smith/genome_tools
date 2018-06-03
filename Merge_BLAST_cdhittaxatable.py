import os
import sys
from subprocess import call

# R_file = sys.argv[1]
# I_file = sys.argv[2]

blast_file = "blast_out.blast"
taxa_file = "cdhit_taxa_table_w_repseq.txt"
database = "GreenGene.fa"

blast_line = "B"
taxa_line = "T"

call(["sort","-k1",blast_file], stdout=open("blast_out_sorted.blast", "w"))
call(["grep",">",database], stdout=open("database_taxa.txt", "w"))

Flag = False
with open("blast_out_sorted.blast", 'r') as B, open(taxa_file, 'r') as T:
	for B_line in B:
		OTU = B_line.split("\t")[0]
		while True:
			try:
				if OTU == T_line.split("\t")[0]:
					print I_line.strip()
					I_line = I.next()
					Flag = True
					break
					else:
						if Flag == True:
							print I_line.strip()
							n += 1
							if n == 3:
								Flag = False
						I_line = I.next()
						continue
				except:
					print "something went wrong"