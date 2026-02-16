## Variant calling pipeline used for Pennycress GBS data 
## CMM Oct-Nov 2024

##1. Demultiplex samples. Script demultiplexes paired-end reads using "barcodes.fasta" a file containing information on the barcode associated with each sample. Raw fastqs are stored in ./raw_fastqs/

demultiplex_pennycressGBS.sh

##2. Index reference, align reads, check alignment scores 

align_pennycressGBS.sh

##3. Call variants 

callSNPs_pennycressGBS.sh

##4. Filter variants

filter_pennycressGBS.sh 
