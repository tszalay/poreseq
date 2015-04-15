if [ $# -lt 2 ]; then
    echo "lastalign [ref fasta] [src fasta]"
    return
fi

# index the ref fasta file
samtools faidx $1
# build a database of it, same name
lastdb $1 $1
# align the fasta file to 
lastal -s 2 -T 0 -Q 0 -a 1 $1 $2 | last-map-probs | maf-convert sam | samtools view -t $1.fai -S -b - | samtools sort - $2
samtools index $2.bam
