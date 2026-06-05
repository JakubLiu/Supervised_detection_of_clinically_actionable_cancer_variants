#!/usr/bin/bash

# usage ./error_rate_reporter.sh reference.fa regions.bed bamfile.bam output.txt

#        bed file structure:

#        1	10001	10001
#        1	10002	10002
#        1	10003	10003
#        1	10004	10004
#        1	10005	10005
#        1	10006	10006
#        1	10007	10007

reference_genome="$1"
bed_file="$2"
bam_file="$3"
output="$4"


bam-readcount -f "$reference_genome" \
            -l "$bed_file" \
            "$bam_file" | \
            awk 'BEGIN {
    OFS="\t";
    print "chr","pos","ref","depth","A","C","G","T","N","error_rate"
}

{
    chr=$1;
    pos=$2;
    ref=$3;
    depth=$4;

    A=C=G=T=N=0;

    for(i=5;i<=NF;i++){
        split($i,a,":");
        base=a[1];
        count=a[2];

        if(base=="A") A=count;
        else if(base=="C") C=count;
        else if(base=="G") G=count;
        else if(base=="T") T=count;
        else if(base=="N") N=count;
    }

    # reference count
    ref_count = 0;
    if(ref=="A") ref_count=A;
    else if(ref=="C") ref_count=C;
    else if(ref=="G") ref_count=G;
    else if(ref=="T") ref_count=T;
    else if(ref=="N") ref_count=N;

    non_ref = A + C + G + T + N - ref_count;

    if(depth > 0){
        error_rate = non_ref / depth;
    } else {
        error_rate = 0;
    }

    print chr,pos,ref,depth,A,C,G,T,N,error_rate;

}' > "$output"

