all.out <- list()
all.out[["MAGE-TAB Version"]] <- "1.1"
all.out[["Investigation Title"]] <- "Assessing the reliability of spike-in normalization for analyses of single-cell RNA sequencing data"
all.out[["Experiment Description"]] <- "This study aims to assess the reliability of spike-in normalization for analyzing single-cell RNA sequencing data. This is done by performing mixture experiments where two different sets of spike-in RNA (ERCC and SIRV) are added separately to each cell (416B or trophoblasts), and generating sequencing libraries using a modified version of the Smart-seq2 protocol. The aim is to measure the variance of the log-ratio of the total counts between the two spike-in sets. This will quantify how precisely the spike-in RNA was added to each well. As a control, addition was also performed with a premixed solution of both spike-ins, to quantify the variability in the log-ratios due to the experimental protocol. The same data can also be used to measure the well-to-well variability in the differences in behaviour between the two spike-in sets."

all.out[["Experimental Design"]] <- c("hardware variation design", "normalization testing design", "operator variation design")
all.out[["Experimental Design Term Source REF"]] <- c("EFO", "EFO", "EFO") 
all.out[["Experimental Design Term Accession Number"]] <- c("EFO_0001767", "EFO_0001771", "EFO_0001772")

all.out[["Experimental Factor Name"]] <- c("spike-in addition", "treatment") 
all.out[["Experimental Factor Type"]] <- c("spike-in addition", "treatment") 
all.out[["Experimental Factor Term Source REF"]] <- c("", "")
all.out[["Experimental Factor Term Accession Number"]] <- c("", "")

all.out[["Person Last Name"]] <- "Lun"
all.out[["Person First Name"]] <- "Aaron"
all.out[["Person Mid Initials"]] <- "TL"
all.out[["Person Email"]] <- "aaron.lun@cruk.cam.ac.uk"
all.out[["Person Phone"]] <- ""
all.out[["Person Fax"]] <- ""      
all.out[["Person Address"]] <- "University of Cambridge Li Ka Shing Centre Robinson Way Cambridge CB2 0RE United Kingdom"
all.out[["Person Affiliation"]] <- "Cancer Research UK Cambridge Institute"
all.out[["Person Roles"]] <- "submitter"

all.out[["Protocol Name"]] <- c("Obtaining 416B cells",
                                "Obtaining trophoblasts",
                                "Culturing 416B cells",
                                "Culturing trophoblasts",
                                "Extracting RNA",
                                "Creating libraries",
                                "Sequencing libraries",
                                "Assigning reads to genes",
                                "Reverse transcription")
all.out[["Protocol Type"]] <- c("sample collection protocol",
                                "sample collection protocol",
                                "growth protocol",
                                "growth protocol",
                                "nucleic acid extraction protocol",
                                "nucleic acid library construction protocol",
                                "nucleic acid sequencing protocol",
                                "high throughput sequence alignment protocol",
                                "conversion protocol")
all.out[["Protocol Term Source REF"]] <- c("EFO", 
                                           "EFO",   
                                           "EFO",   
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO",
                                           "EFO")
all.out[["Protocol Term Accession Number"]] <- c("EFO_0005518",
                                                 "EFO_0005518",
                                                 "EFO_0003789",
                                                 "EFO_0003789",
                                                 "EFO_0002944",
                                                 "EFO_0004184",
                                                 "EFO_0004170",
                                                 "EFO_0004917",
                                                 "EFO_0005520")

all.out[["Protocol Description"]] <- c("??? Obtaining 416B cells",
                                       "??? Obtaining trophoblasts",
                                       "??? Culturing 416B cells",
                                       "??? Culturing trophoblasts",
                                       "Cells were sorted into individual wells of a 96-well microtiter plate. Each well contained 2.3 ul of lysis buffer with RNAse inhibitor (Ambion) in a 0.2% (v/v) Triton X-100 solution.",
                                       "Preamplification was performed in a total volume of 27 ul that contained 13.5 ul of HiFi Hotstart ReadyMix (2x; KAPA Biosystems) and 0.1 uM of IS PCR primer (Sigma-Aldrich). After 23 cycles of amplification, samples were cleaned with 80% (v/v) of Ampure beads (Beckman Coulter). Sequencing libraries were prepared using the Nextera XT DNA sample preparation kit (Illumina).",
                                       "Libraries for batch 20160113 were sequenced on an Illumina HiSeq 2500 to obtain 125 bp single-end reads. Libraries for batch 20160325 were sequenced on an Illumina HiSeq 4000 to obtain 50 bp single-end reads. Libraries for batches 20160906 and 20170201 were sequenced on an Illumina HiSeq 4000 to obtain 75 bp paired-end reads.",
                                       "Reads were aligned to the mm10 build of the mouse genome using subread v1.5.0 in RNA-seq mode with unique mapping. The number of reads or read pairs mapped to the exonic regions of each gene was then counted for each library, using the featureCounts function in Rsubread v1.24.1 with Ensembl GRCm38 version 82. Only alignments with mapping quality scores above 10 were considered during counting.",
                                       "Reverse transcription was performed in a final volume of 13.2 ul per well, containing 1 uM of oligo-dT (Sigma-Aldrich), 1.04 mM of each dNTP (ThermoFisher), 100 U of SuperScript II retrotranscriptase (Invitrogen/ThermoFisher), 5 U of RNase inhibitor (Ambion), 5 mM of DTT, 1 M of Betaine (Sigma-Alrich), 6 mM of MgCl$_2$ (Ambion) and 1 uM of TSO primer (Exiqon). Spike-in RNA was added such that each well contained a final volume of 0.1 ul of a 1:3,000,000 dilution of the ERCC RNA Spike-In Mix (Invitrogen/ThermoFisher) and 0.12 ul of a 1:3,000,000 dilution of the Spike-in RNA Variant (SIRV) Control Mix E0 (Lexogen).")

all.out[["Protocol Hardware"]] <- c("",
                                    "", 
                                    "", 
                                    "", 
                                    "",
                                    "",
                                    "Illumina HiSeq 2500, Illumina HiSeq 4000",
                                    "",
                                    "")
all.out[["Protocol Software"]] <- c("", 
                                    "",
                                    "",
                                    "",
                                    "",
                                    "",
                                    "",
                                    "(R)subread",
                                    "")
                                 
all.out[["Term Source Name"]] <- "EFO"
all.out[["Term Source File"]] <- "http://www.ebi.ac.uk/efo/"
all.out[["Term Source Version"]] <- ""
all.out[["SDRF File"]] <- "sdrf.tsv"
all.out[["Public Release Date"]] <- "2017-05-31"
all.out[["Comment[AEExperimentType]"]] <- "RNA-seq of coding RNA"

unlink("idf.tsv")
for (x in names(all.out)) {
    write(file="idf.tsv", paste0(c(x, all.out[[x]]), collapse="\t"), append=TRUE)
}

##########################################################################3
# Merging the two SDRF files together.

sdrf.calero <- read.table("sdrf_Calero.tsv", header=TRUE, check.names=FALSE, comment="", sep="\t", quote="", stringsAsFactors=FALSE)
sdrf.liora <- read.table("sdrf_Liora.tsv", header=TRUE, check.names=FALSE, comment="", sep="\t", quote="", stringsAsFactors=FALSE)
stopifnot(identical(colnames(sdrf.calero), colnames(sdrf.liora)))
write.table(file="sdrf.tsv", rbind(sdrf.calero, sdrf.liora), row.names=FALSE, quote=FALSE, sep="\t")

##########################################################################3
# Computing the spike-in quantity per well.

# ERCC data taken from https://www.thermofisher.com/order/catalog/product/4456740,
# under the link "ERCC Controls Analysis: ERCC RNA Spike-In Control Mixes (English)".
ercc.vol <- 0.1 # in uL
ercc.dil <- 3e6
ercc.data <- read.table("cms_095046.txt", header=TRUE, check.names=FALSE, sep="\t", stringsAsFactors=FALSE)
ercc.id <- ercc.data[,"ERCC ID"]
ercc.quant <- ercc.vol/ercc.dil * ercc.data[,"concentration in Mix 1 (attomoles/ul)"]

# SIRV data from https://www.lexogen.com/sirvs/downloads/,
# under the link (SIRV sequence design overview (XLSX)).
# We compute the total molarity across all transcripts for each gene (from femtomoles -> attomoles)
sirv.vol <- 0.12 # in uL
sirv.dil <- 3e6
sirv.data <- c(SIRV1=8, SIRV2=6, SIRV3=11, SIRV4=7, SIRV5=12, SIRV6=18, SIRV7=7) * 1000
sirv.id <- names(sirv.data)
sirv.quant <- sirv.vol/sirv.dil * sirv.data

# Sanity check:
# library(edgeR)
# full <- readRDS("../Calero/trial_20160113/analysis/full.rds")
# plot(ercc.quant, rowSums(full$counts[full$genes$spike1,])[ercc.id])
# points(sirv.quant, rowSums(full$counts[full$genes$spike2,])[sirv.id], col="red")

spike.dir <- "spike-data"
dir.create(spike.dir, showWarnings=FALSE)
write.table(file=file.path(spike.dir, "spikes.txt"), sep="\t", quote=FALSE, row.names=FALSE,
            data.frame(Name=c(ercc.id, sirv.id), "Attomole/well"=c(ercc.quant, sirv.quant), check.names=FALSE))



