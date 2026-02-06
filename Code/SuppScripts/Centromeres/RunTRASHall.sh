#!/bin/bash
## Run the command from the current working directory. Launch qsub from the library directory
#$ -cwd
#$ -pe smp 12
#$ -N TRASH
#$ -m aes
## Job Starts Here

source ~/.bashrc
mamba activate TRASH_centromeres
export PATH=$PATH:~/Pennycress/centromere/TRASH/TRASH/
export R_LIBS=~/Pennycress/centromere/TRASH/TRASH/libs


bash ~/Pennycress/centromere/TRASH/TRASH/TRASH_run.sh --simpleplot ~/Pennycress/centromere/Tarvensevar_AK34W_873_v1.0.fa --o ~/Pennycress/centromere/TRASH/AK34W_20241007

bash ~/Pennycress/centromere/TRASH/TRASH/TRASH_run.sh --simpleplot ~/Pennycress/centromere/Tarvensevar_Ames32873_875_v1.0.fa --o ~/Pennycress/centromere/TRASH/Ames32873_20241007

bash ~/Pennycress/centromere/TRASH/TRASH/TRASH_run.sh --simpleplot ~/Pennycress/centromere/Tarvensevar_LorettoMN_876_v1.0.fa --o ~/Pennycress/centromere/TRASH/Loretto_20241007

bash ~/Pennycress/centromere/TRASH/TRASH/TRASH_run.sh --simpleplot ~/Pennycress/centromere/Tarvensevar_MN106_872_v4.0.fa --o ~/Pennycress/centromere/TRASH/MN106_20241007

bash ~/Pennycress/centromere/TRASH/TRASH/TRASH_run.sh --simpleplot ~/Pennycress/centromere/Tarvensevar_MN134_877_v1.0.fa --o ~/Pennycress/centromere/TRASH/MN134_20241007

bash ~/Pennycress/centromere/TRASH/TRASH/TRASH_run.sh --simpleplot ~/Pennycress/centromere/Tarvensevar_PI650286_878_v1.0.fa  --o ~/Pennycress/centromere/TRASH/PI650286_20241007

bash ~/Pennycress/centromere/TRASH/TRASH/TRASH_run.sh --simpleplot ~/Pennycress/centromere/Tarvensevar_Tibet33_879_v1.0.fa --o ~/Pennycress/centromere/TRASH/Tibet33_20241004
