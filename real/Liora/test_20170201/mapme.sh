ispet=1
fastq=($(ls cram/*.cram))
genome=/lustre/jmlab/resources/genomes/subread/mm10_ERCC_SIRV

source ${HOME}/Code/mapping/multi_align.sh
