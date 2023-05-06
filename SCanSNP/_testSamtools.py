import pandas as pd
import subprocess





print('wait')


samtools mpileup --no-BAQ --min-MQ 20 --ff 3844 --no-output-ins --no-output-ends \
--ignore-RG -l <(cut -f1,2 /tmp/temp.loci.tsv) /data/possorted_genome_bam.bam


samtools mpileup --no-BAQ --min-MQ 20 --ff 3844 -Q 3  --output-extra CB --no-output-ins --no-output-ends \
--ignore-RG -l <(cut -f1,2 /tmp/temp.loci.tsv) /data/possorted_genome_bam.bam | cut -f 1,2,5,7



sed -E 's/^(([^<>\t]*\t){4}[^<>\t]*)[<>\t]([^\t]*)(\t[^\t]*)?/\1@\3/g'



How can i reshape this
chr1    209830854       GG>G      CGTCAAATCCTAGCCT-1,ATAGAGAAGCCATCCG-1,TGAGGAGTCTTCTGTA-1,CTCCATGAGGACTGGT-1
chr2    22       <C      TTTGTTGTCGCCATAA-1,AGGCATTCAATTGAAG-1

Into this
chr1    209830854       G      CGTCAAATCCTAGCCT-1
chr1    209830854       G      ATAGAGAAGCCATCCG-1
chr1    209830854       >      TGAGGAGTCTTCTGTA-1
chr1    209830854       G      CTCCATGAGGACTGGT-1
chr2    22       A      TTTGTTGTCGCCATAA-1,AGGCATTCAATTGAAG-1
chr2    22       <      AGGCATTCAATTGAAG-1



sed -E 's/([A-Z]+)([><]*)([A-Z]*)/\1\t\2\t\3/g' test2.tsv | \
awk -F'\t' '{split($4, a, ","); split($3, b, ""); for (i=1; i<=length(a); i++) print $1"\t"$2"\t"b[i]"\t"a[i]}' test2.tsv | \
awk '$ ~ /[[:alpha:]]/' | grep -f barcodes.tsv


time awk -F '\t' '{split($4, a, ","); split($3, b, ""); for (i=1; i<=length(a); i++) print $1"\t"$2"\t"b[i]"\t"a[i]}' testChunk.tsv |\
awk '$3 ~ /[[:alpha:]]/'| grep -f barcodes.tsv | wc -l 


##Final version
cut -f1,2 /tmp/temp.loci.tsv > /tmp/temp.loci.parsed.tsv
time samtools mpileup --no-BAQ --min-MQ 20 --ff 3844 -Q 3  --output-extra CB --no-output-ins --no-output-ends --ignore-RG -l /tmp/temp.loci.parsed.tsv /data/possorted_genome_bam.bam \
| cut -f 1,2,5,7 |  awk -F '\t' '{split($4, a, ","); split($3, b, ""); for (i=1; i<=length(a); i++) print $1"\t"$2"\t"b[i]"\t"a[i]}' | \
awk -v OFS='\t' '$3 ~ /[[:alpha:]]/ { $3 = toupper($3); print }'| grep -f /tmp/barcodes.tsv 




time samtools view -F 3844 --region-file /tmp/temp.loci.parsed.tsv --min-MQ --min-MQ  --tag-file "CB:/tmp/barcodes.tsv"  \
 /data/possorted_genome_bam.bam >/dev/null


time samtools view -@ 5 -h -F 3844 --keep-tag CB -r --region-file /tmp/temp.loci.parsed.tsv --tag-file "CB:/tmp/barcodes.tsv" -b \
 /data/possorted_genome_bam.bam  > /data/possorted_genome_bam.filtered.bam && samtools index /data/possorted_genome_bam.filtered.bam



samtools view -@ 5 -h -F 3844 --keep-tag RG --region-file /tmp/temp.loci.parsed.tsv --tag-file "CB:/tmp/barcodes.tsv" \
 /data/possorted_genome_bam.bam

 
cat /tmp/temp.loci.parsed.tsv | sed  's/\t/@/g'