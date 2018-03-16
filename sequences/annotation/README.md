# Annotation sources

- ERCC92.gtf was taken from http://www.thermofisher.com/order/catalog/product/4456739
- CBFB-MYH11-mcherry.gtf was manually constructed, based on information from Fernando Calero
- SIRV_C_150601a.gtf was taken from https://www.lexogen.com/sirvs/
- mm10.gtf was constructed from http://ftp.ensembl.org/pub/release-82/gtf/mus_musculus/

```sh
zcat Mus_musculus.GRCm38.82.gtf.gz | \
    sed -r "s/^([0-9MXY])/chr\1/" | \
    sed "s/^chrMT/chrM/g" | \
    awk '$3 == "exon"' > mm10.gtf
```
