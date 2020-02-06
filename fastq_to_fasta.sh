awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' reads.fastq > reads.fa


-or-

sed -n '1~4s/^@/>/p;2~4p' in.fastq > out.fasta
