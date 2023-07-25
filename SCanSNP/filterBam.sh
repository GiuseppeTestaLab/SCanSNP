
samtools view -h -D CB:/data/barcodes.tsv --keep-tag CB,RG /data/possorted_genome_bam.bam 

samtools view -h -D CB:/data/barcodes.tsv  /data/possorted_genome_bam.bam | head -n 300