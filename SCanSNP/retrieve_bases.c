// retrieve_bases.c
#include <htslib/sam.h>
#include <stdio.h>


void retrieve_bases(char* bam_file, int* positions, int n_positions) {
	samFile *bam;
	bam_hdr_t *header;
	bam1_t *read = bam_init1();
	const char* chrom = "chr1"; // change this to the chromosome you want to retrieve bases from
	
	bam = sam_open(bam_file, "r");
	header = sam_hdr_read(bam);
	
	for (int i = 0; i < n_positions; i++) {
		hts_itr_t *iter = sam_itr_querys(header, chrom, positions[i]-1, positions[i]);
		while (sam_itr_next(bam, iter, read) >= 0) {
			uint8_t *seq = bam_get_seq(read);
			int32_t pos = read->core.pos;
			if (pos == positions[i]-1) {
				char base = seq_nt16_str[bam_seqi(seq, 0)];
				printf("%d\t%c\n", pos+1, base);
			}
		}
		hts_itr_destroy(iter);
	}

	bam_destroy1(read);
	sam_close(bam);
	bam_hdr_destroy(header);
}
