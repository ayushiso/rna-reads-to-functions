#!/bin/bash

# remove old logs
rm -r logs
mkdir logs

# generate fastqc reports of raw reads
output=init_fastqc_reports/
mkdir -p $output
for file in RNAseq_fastqs/*.fastq.gz; do
    fastqc ${file} -t 16 -o $output &>> logs/fastqc.log
done

mkdir -p mqc_init

multiqc $output -o mqc_init --force

# filter and conduct another round of quality control

reports=final_fastqc_reports/
mkdir -p filt_reads
mkdir -p $reports

for file in RNAseq_fastqs/*.fastq.gz; do
    fname="$(basename "$file")"
    cutadapt -m 10 -q 20 -o filt_reads/$fname $file &>> logs/cutadapt.log
    fastqc filt_reads/$file -t 16 -o $reports &>> logs/fastqc_final.log 
done

mkdir -p mqc_final
multiqc $output -o mqc_final --force


# generate STAR genome indices
mkdir -p genome_idx
cd genome_idx 
wget ftp://ftp.ensembl.org/pub/release-99/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz &> dl.log 
wget ftp://ftp.ensembl.org/pub/release-99/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz &>> dl.log
gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
gunzip Saccharomyces_cerevisiae.R64-1-1.99.gtf.gz
STAR --runThreadN 2 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa --sjdbGTFfile Saccharomyces_cerevisiae.R64-1-1.99.gtf --genomeSAindexNbases 10 &> star.log
cd ..

# run STAR and get counts
mkdir -p STAR_output
for file1 in filt_reads/*forward_file.fastq.gz; do
    file2=$(echo $file1 | sed 's/forward/reverse/g')
    prefix=$(echo $file1 | awk -F '/' $'{print $2}' | awk -F '_' $'{print $1 "_" $2 "_" $3}')
    STAR --runThreadN 2 --genomeDir genome_idx --readFilesCommand zcat --readFilesIn $file1 $file2 --outSAMtype BAM Unsorted --outFileNamePrefix STAR_output/$prefix- &>> logs/star.log
    htseq-count -f bam STAR_output/$prefix-Aligned.out.bam genome_idx/Saccharomyces_cerevisiae.R64-1-1.99.gtf > $prefix.counts
done
