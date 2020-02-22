# rna-reads-to-functions
Workflow in Bash, R and Python to go from raw RNA reads to mapped counts, and from mapped counts to differentially expressed genes with basic functional analysis

1. For going from reads to counts, use `reads_to_counts.sh` 
    * requires that you have paired RNASeq reads in RNASeq_fastqs/
    * software requirements: 
       1. FastQC: 0.11.9
       2. MultiQC: 1.8
       3. Cutadapt: 2.6
       4. STAR: 2.7.3a
 
 2. For differential gene analysis from counts, use `diff_exp.R`
     * **HIGHLY RECOMMENDED** to use RStudio for this script
     * requires that you have counts as \*.counts files in RNASeq_counts/
     
 3. For basic ontology annotation of differentially expressed genes, use `ontology.py`
     * basic command: `python ontology.py path\to\upregulated\genes path\to\downregulated\genes onto_type`
     * onto_type can be C, F, or P (GO slim types)
     
