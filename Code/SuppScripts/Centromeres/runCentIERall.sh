#!/bin/bash
## Run the command from the current working directory. Launch qsub from the library directory
#$ -cwd
#$ -pe smp 12
#$ -N centier_no_gff
#$ -m aes
## Job Starts Here

source ~/.bashrc
conda activate centier

export PATH=$PATH:~centier/genometools-1.6.5/bin
export PATH=$PATH:~centier/CentIER

centIER.py ~Pennycress/centromere/Tarvensevar_AK34W_873_v1.0.fa -o ./AK34W_CentIER_results_no_gff

centIER.py ~Pennycress/centromere/Tarvensevar_Ames32873_875_v1.0.fa -o ./Ames_CentIER_results_no_gff

centIER.py ~Pennycress/centromere/Tarvensevar_LorettoMN_876_v1.0.fa -o ./Loretto_CentIER_results_no_gff

centIER.py ~Pennycress/centromere/Tarvensevar_MN106_872_v4.0.fa -o ./MN106_CentIER_results_no_gff

centIER.py ~Pennycress/centromere/Tarvensevar_MN134_877_v1.0.fa -o ./MN134_CentIER_results_no_gff

centIER.py ~Pennycress/centromere/Tarvensevar_PI650286_878_v1.0.fa -o ./PI650286_CentIER_results_no_gff

centIER.py ~Pennycress/centromere/Tarvensevar_Tibet33_879_no_scaffold_v1.0.fa -o ./Tibet_CentIER_results_no_gff
