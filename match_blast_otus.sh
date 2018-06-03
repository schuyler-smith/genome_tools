
BLASTDB='ITSdb.findley.fasta'




head -n 1 cdhit_taxa_table_w_repseq.txt > cdhit_taxa_table.txt

awk {'print $1"\t"$2'} cdhit_taxa_table_w_repseq.txt | tail -n +2 > otu_list.txt

awk {'print $1'} otu_list.txt | while read -r line
	do grep -w $line out.blast
done > blast_output.txt

awk {'print $2'} blast_output.txt | while read -r line;
	do grep $line $BLASTDB | awk {'print $2'}
done > blast_taxa.txt

paste <(cat otu_list.txt ) <(cat blast_taxa.txt ) >> cdhit_taxa_table.txt
sed -i 's/Root;//g' cdhit_taxa_table.txt
sed -i 's/;/\t/g' cdhit_taxa_table.txt

# removes third column?
#sed -ri 's/(\s+)?\S+//3' cdhit_taxa_table.txt

#removes duplicate rows based on column
#awk '!seen[$0]++' otu_table.txt > otu_table.unique.txt