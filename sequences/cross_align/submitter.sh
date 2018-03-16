# Aligner submission:
rm bam/mm10/align*
bsub -R "rusage[mem=10000]" -n 6 -e "bam/mm10/align.err" -o "bam/mm10/align.out" ./readmaker.py -i ~/lustre/genomes/ERCC ~/lustre/genomes/SIRV -f /lustre/reference_data/genomes/Mus_musculus_GRCm38/chroms/chr*.fa -o bam/mm10/ERCC bam/mm10/SIRV

rm bam/SIRV/align*
bsub -R "rusage[mem=10000]" -n 6 -e "bam/SIRV/align.err" -o "bam/SIRV/align.out" ./readmaker.py -i ~/lustre/genomes/ERCC ~/lustre/genomes/mm10 -f ~/lustre/genomes/sequences/SIRV_150601a.fasta -o bam/SIRV/ERCC bam/SIRV/mm10

rm bam/ERCC/align*
bsub -R "rusage[mem=10000]" -n 6 -e "bam/ERCC/align.err" -o "bam/ERCC/align.out" ./readmaker.py -i ~/lustre/genomes/SIRV ~/lustre/genomes/mm10 -f ~/lustre/genomes/sequences/ERCC92.fa -o bam/ERCC/SIRV bam/ERCC/mm10

# Diagnostic submission:
rm results/log*
bsub -R "rusage[mem=2000]" -n 1 -e "results/log.out" -o "results/log.err" ./mapcheck.py -b bam/mm10/SIRV/*.bam -s ~/lustre/annotation/mm10.gtf -a ~/lustre/annotation/SIRV_C_150601a.gtf -o results/mm10-SIRV.txt
bsub -R "rusage[mem=2000]" -n 1 -e "results/log.out" -o "results/log.err" ./mapcheck.py -b bam/mm10/ERCC/*.bam -s ~/lustre/annotation/mm10.gtf -a ~/lustre/annotation/ERCC92.gtf -o results/mm10-ERCC.txt

bsub -R "rusage[mem=2000]" -n 1 -e "results/log.out" -o "results/log.err" ./mapcheck.py -b bam/SIRV/mm10/*.bam -s ~/lustre/annotation/SIRV_C_150601a.gtf -a ~/lustre/annotation/mm10.gtf -o results/SIRV-mm10.txt
bsub -R "rusage[mem=2000]" -n 1 -e "results/log.out" -o "results/log.err" ./mapcheck.py -b bam/SIRV/ERCC/*.bam -s ~/lustre/annotation/SIRV_C_150601a.gtf -a ~/lustre/annotation/ERCC92.gtf -o results/SIRV-ERCC.txt

bsub -R "rusage[mem=2000]" -n 1 -e "results/log.out" -o "results/log.err" ./mapcheck.py -b bam/ERCC/mm10/*.bam -s ~/lustre/annotation/ERCC92.gtf -a ~/lustre/annotation/mm10.gtf -o results/ERCC-mm10.txt
bsub -R "rusage[mem=2000]" -n 1 -e "results/log.out" -o "results/log.err" ./mapcheck.py -b bam/ERCC/SIRV/*.bam -s ~/lustre/annotation/ERCC92.gtf -a ~/lustre/annotation/SIRV_C_150601a.gtf -o results/ERCC-SIRV.txt
