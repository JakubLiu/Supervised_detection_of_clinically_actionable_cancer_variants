#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>


void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s <input.bam> <output.bam> <Q_threshold>\n", prog);
    exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {


    if (argc != 4) usage(argv[0]);    // print error message if the command line arguments are wrong

    const char *in_bam  = argv[1];   //input bam
    const char *out_bam = argv[2];   // output bam
    int Q_THRESHOLD = atoi(argv[3]);


    samFile *in = sam_open(in_bam, "r");  // open the input bam file for reading
    if (!in) {
        fprintf(stderr, "Error: cannot open input BAM %s\n", in_bam);
        return EXIT_FAILURE;
    }

    samFile *out = sam_open(out_bam, "wb");  // open the output bamfile for writing
    if (!out) {
        fprintf(stderr, "Error: cannot open output BAM %s\n", out_bam);
        return EXIT_FAILURE;
    }

    bam_hdr_t *hdr = sam_hdr_read(in);  // read the bam header from the input file
    sam_hdr_write(out, hdr);            // write the bam header from the input file to the output file

    bam1_t *b = bam_init1();    // this will hold the alignment record


    while (sam_read1(in, hdr, b) >= 0) {   // loop over the records in the bam file

        uint8_t *seq  = bam_get_seq(b);
        uint8_t *qual = bam_get_qual(b);
        int readlen = b->core.l_qseq;

        /* QUAL=255 means missing */
        if (qual[0] != 255) {
            for (int i = 0; i < readlen; i++) {      // loop over the bases in a given read
                if (qual[i] < Q_THRESHOLD) {       // if the basecalling quality is below the given threshold, mask it as N
                    bam_set_seqi(seq, i, seq_nt16_table['N']);
                }
            }
        }

        sam_write1(out, hdr, b);
    }

    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(in);
    sam_close(out);

    return EXIT_SUCCESS;
}
