#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>



char* extract_bases_from_bam(const char* chrom, int pos, const char* bam_fname) {
    sam_hdr_t *hdr;
    samFile* bam_file = sam_open(bam_fname, "r"); // open BAM file
    if (bam_file == NULL) {
        fprintf(stderr, "Error opening BAM file: %s\n", bam_fname);
        return NULL;
    }
    hts_idx_t* idx = sam_index_load(bam_file, bam_fname); // load BAM index
    if (idx == NULL) {
        fprintf(stderr, "Error loading BAM index: %s\n", bam_fname);
        sam_close(bam_file);
        return NULL;
    }
    hts_itr_t* itr = sam_itr_queryi(idx, HTS_IDX_START, pos, pos+1); // query BAM file
    if (itr == NULL) {
        fprintf(stderr, "Error querying BAM file: %s:%d\n", chrom, pos+1);
        hts_idx_destroy(idx);
        sam_close(bam_file);
        return NULL;
    }
    char* bases = NULL; // initialize bases string
    int num_bases = 0; // initialize number of bases
    hdr = sam_hdr_read(bam_file); // read BAM header
    bam1_t* rec = bam_init1(); // initialize SAM record
    int num_records = 0; // initialize number of records
    while (sam_itr_next(bam_file, itr, rec) >= 0) { // loop through records
        num_records++;
        if (strcmp(hdr->target_name[rec->core.tid], chrom) == 0 // check chromosome
            && rec->core.pos <= pos && rec->core.pos + bam_cigar2rlen(rec->core.n_cigar, bam_get_cigar(rec)) > pos) {
            const uint32_t *cigar = bam_get_cigar(rec); // check position
            const uint8_t* seq = bam_get_seq(rec); // get sequence
            uint8_t base_num = bam_seqi(bam_get_seq(rec), pos); // get base at position
            if (base_num != 15) { // check for N or other ambiguous base
                char c = seq_nt16_str[bam_seqi(seq, pos)]; // convert to character
                bases = realloc(bases, num_bases + 2); // grow bases string
                bases[num_bases++] = c; // add base to string
            }
        }
    }
    printf("Number of records retrieved: %d\n", num_records);
    if (num_bases == 0) {
        printf("No bases extracted\n");
    }
    else {
        bases[num_bases] = '\0'; // terminate bases string
        printf("Bases extracted: %s\n", bases);
    }
    bam_destroy1(rec); // free SAM record
    sam_hdr_destroy(hdr); // free BAM header
    hts_itr_destroy(itr);
    return bases;
}



int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage: %s <chromosome> <position> <BAM file>\n", argv[0]);
        return 1;
    }
    const char* chrom = argv[1];
    int pos = atoi(argv[2]);
    const char* bam_fname = argv[3];
    char* bases = extract_bases_from_bam(chrom, pos, bam_fname);
    if (bases == NULL) {
        printf("Error extracting bases\n");
        return 1;
    }
    printf("Bases at %s:%d: %s\n", chrom, pos+1, bases);
    free(bases);
    return 0;
}

//gcc -o extract_bases extract_bases.c -I/usr/include/ -L/lib/x86_64-linux-gnu/pkgconfig/ -lhts

