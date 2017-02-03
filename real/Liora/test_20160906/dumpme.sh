bsub -R "rusage[mem=5000]" -o dump.out -e dump.err bash ~/Code/mapping/sanger_dump.sh PAIRED
