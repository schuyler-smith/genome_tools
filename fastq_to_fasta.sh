awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' reads.fastq > reads.fa


-or-

for file in $(ls fastq); do file_name="$(basename -- $file)"; file_name=$(echo "$file_name" | cut -f 1 -d '.'); sed -n '1~4s/^@/>/p;2~4p' fastq/${file_name}*.fastq > fasta/${file_name}.fasta; done
