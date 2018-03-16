# Sequence sources

- ERCC sequences were obtained from https://www.thermofisher.com/order/catalog/product/4456739
- SIRV sequences were obtained from https://www.lexogen.com/wp-content/uploads/2015/08/SIRV_Sequences_150601.zip
- CBFB-MYH11-mcherry sequence was supplied by Fernando Calero, CIMR (Jan 21, 2016)
- mm10 sequence was obtained from UCSC (via `/lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa`)

# Genome builds 

This combines mm10 with both the ERCC and SIRV spike-ins:

```sh
#bsub -e log.err -o log.out -R "rusage[mem=16000]" \
subread-buildindex -o mm10_ERCC_SIRV /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa \
    /lustre/jmlab/resources/genomes/sequences/ERCC92.fa \
    /lustre/jmlab/resources/genomes/sequences/SIRV_150601a.fasta 
```

... and an extra induced oncogene:

```sh
#bsub -e log.err -o log.out -R "rusage[mem=16000]" \
subread-buildindex -o mm10_ERCC_SIRV_onco /lustre/reference_data/mib-cri/reference_genomes/mus_musculus/mm10/fasta/mmu.mm10.fa \
    /lustre/jmlab/resources/genomes/sequences/ERCC92.fa \
    /lustre/jmlab/resources/genomes/sequences/SIRV_150601a.fasta \
    CBFB-MYH11-mcherry.fa
```

Also building indices for the spike-ins by themselves, to assess specificity:

```sh
subread-buildindex -o ERCC /lustre/jmlab/resources/genomes/sequences/ERCC92.fa
subread-buildindex -o SIRV /lustre/jmlab/resources/genomes/sequences/SIRV_150601a.fasta 
```
