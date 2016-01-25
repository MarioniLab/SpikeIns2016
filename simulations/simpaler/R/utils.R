findChr <- function(ids, txdb, input="GENEID", output="TXCHROM")
# This gets the chromosome ID for each gene, for filtering mitochondrial genes.
{
	suppressMessages(loc <- select(txdb, keytype=input, columns=output, keys=ids))
        return(loc$TXCHROM[match(ids, loc$GENEID)])
}

