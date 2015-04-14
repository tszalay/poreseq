if [ $# -lt 1 ]; then
    echo "selfalign [input fasta]"
    return
fi

# index the fasta file
samtools faidx $1
# build a database of it, same name
lastdb $1 $1
# align the fasta file to 
lastal -p nanopore.mat $1 $1 | maf-convert sam | samtools view -t $1.fai -S -b - | samtools sort - $1
samtools index $1.bam
