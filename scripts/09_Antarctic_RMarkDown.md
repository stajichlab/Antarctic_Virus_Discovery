R Analysis of Antarctic Metagenome Viral Communities
================
Cassie Ettinger

# Load Libraries

``` r
library(GGally)
library(tidyverse)
library(network)
library(vroom)
library(phyloseq)
library(microbiome)
library(vegan)
library(pairwiseAdonis)
library(patchwork)
library(ggnewscale)
library(ggforce)
library(geosphere)
library(rnaturalearth)
library(sf)
library(rnaturalearthdata)
library(biobroom)
```

# vOTU processing

Following some of Christian Santos Medellin’s scripts:
<https://github.com/cmsantosm/SpatioTemporalViromes/blob/master/Processing/Scripts/votu_filtering.Rmd>

``` r
# run on UCR cluster get file names
dir <- "data/coverage/"
cov.files <- list.files(path = dir, pattern = "*.tsv", full.names = T)
cov.files
sampleid <- list.files(path = dir, pattern = "*.tsv", full.names = F) %>%
    str_remove("_mapped_coverge.tsv")  #whoops! mispelled coverage
sampleid

# open coverage files
cov.list <- lapply(cov.files, read.table, sep = "\t", header = F,
    col.names = c("contig", "start", "end", "coverage"))
names(cov.list) <- sampleid


# Christian's function to get coverage for each vOTUs in
# each bed file and filter out instances in which coverage
# is < 0.75

get_coverage <- function(df) {
    df %>%
        mutate(nbase = end - start) %>%
        group_by(contig) %>%
        mutate(length = sum(nbase)) %>%
        mutate(perc_seq = nbase/length) %>%
        filter(coverage > 0) %>%
        summarise(total_coverage = sum(perc_seq)) %>%
        filter(total_coverage >= 0.75)
}


# Apply the function to all coverage files

cov.list.proc <- lapply(cov.list, get_coverage)
cov.list.proc


# Make the list into a data frame

cov.df <- plyr::ldply(cov.list.proc, function(x) x) %>%
    rename(SampleID = ".id")

cov.df

# saveRDS(cov.df, 'results/processed_data/cov.df.RDS')

# Convert it to a matrix

cov.mtx <- cov.df %>%
    spread(key = SampleID, value = total_coverage) %>%
    as.data.frame()
row.names(cov.mtx) <- cov.mtx$contig
cov.mtx <- cov.mtx[, -1]
cov.mtx

# Remove NAs

cov.mtx <- !is.na(cov.mtx)

# saveRDS(cov.mtx, 'results/processed_data/cov.mtx.RDS')

cov.mtx <- readRDS("results/processed_data/cov.mtx.RDS")

# Get the tpmean table generated from BAMM and keep only
# those instances where coverage >= 0.75 we had a few
# samples missing that we added back in later, hence the
# need to join files here
tpmean1 <- vroom("data/relative_abundance/output_file_tpmean_072622.tsv")
tpmean2 <- vroom("data/relative_abundance/output_file_tpmean_missing_080422.tsv")

tpmean <- full_join(tpmean1, tpmean2)

tpmean
rownames(tpmean) <- tpmean$`#contig`
tmp <- colnames(tpmean)
tmp <- str_remove(tmp, "_sortedIndexed.bam")
tmp <- str_remove(tmp, "data/fastq/bbmap/")
colnames(tpmean) <- tmp
tpmean.75 <- tpmean[match(row.names(cov.mtx), row.names(tpmean)),
    match(colnames(cov.mtx), colnames(tpmean))]
tpmean.75 <- cov.mtx * tpmean.75
rownames(tpmean.75) <- row.names(cov.mtx)
tpmean.75.tidy <- tpmean.75 %>%
    as.data.frame() %>%
    mutate(OTU_ID = row.names(.)) %>%
    gather(key = "SampleID", value = "Abundance", -OTU_ID) %>%
    select(SampleID, everything())


# Save formatted tables

# saveRDS(tpmean.75,
# 'results/processed_data/antvir_tpm_75_mtx.RDS')
# saveRDS(tpmean.75.tidy,
# 'results/processed_data/antvir_tpm_75_tidy.RDS')
```

# VContact2 processing

``` r
# this script will take in the vcontact2 results and
# connect them back to our samples/summarize them

# load in data
otu.tidy <- readRDS("results/processed_data/antvir_tpm_75_tidy.RDS")  #relative abundance data
genome.ov <- read.table("data/vcontact2/genome_by_genome_overview_ictv.csv",
    header = T, sep = ",")  #results from vcontact2
ntwk <- read.table("data/vcontact2/c1_ictv.ntw", header = F,
    sep = " ", col.names = c("OTU1", "OTU2", "Score"))  #network from vcontact2
samples <- read.csv("data/samples.csv")  #file with sample names
ictv_names <- vroom("data/vcontact2/INPHARED_1Aug2022_data.tsv")

# define the source of each vOTU, either Antarctica (AV) or
# from the reference database INPHARED (refseq)
genome.ov <- genome.ov %>%
    mutate(Source = ifelse(str_detect(Genome, paste(samples$SAMPLE,
        collapse = "|")), "AV", "refseq"))

genome.ov <- left_join(genome.ov, ictv_names, by = c(Genome = "Accession"))

genome.ov <- genome.ov %>%
    mutate(Order = ifelse(is.na(Order.y), Order.x, Order.y)) %>%
    mutate(Family = ifelse(is.na(Family.y), Family.x, Family.y)) %>%
    mutate(Genus = ifelse(is.na(Genus.y), Genus.x, Genus.y))


# check that # contigs matches number expected vOTUs
# (76986) genome.ov %>% filter(Source == 'AV') %>%
# group_by(Genome) %>% count()

# Generate a data frame that specifies the composition of
# each viral cluster in terms of source of network nodes
# (ie. does the cluster have sequences from Antarctica?)
# ClstrComp = cluster composition, either from Antartica
# (AV), reference data(INPHARED/refseq) or both
clstr.source <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    mutate(AV = str_detect(Genome, paste(samples$SAMPLE, collapse = "|"))) %>%
    group_by(VC, Size) %>%
    summarise(pAV = sum(AV)/n()) %>%
    mutate(ClstrComp = case_when(pAV == 0 ~ "refseq", pAV ==
        1 ~ "AV", TRUE ~ "both"))

# Let's create a data frame with the consensus viral
# taxonomy for each cluster.  the strategy is to check each
# taxonomic rank and see if there's only one assignment or
# a mix.  Some refseq entrees have unassigned ranks so we
# are ignoring those and deciding the proper classification
# based on the rest of the VC members
clstr.phyla <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(Phylum != "Unassigned") %>%
    group_by(VC, Phylum) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(Phylum = ifelse(Duplicates > 1, "Mixed", as.character(Phylum))) %>%
    group_by(VC, Phylum) %>%
    count() %>%
    select(-n)

clstr.class <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(Class != "Unassigned") %>%
    group_by(VC, Class) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(Class = ifelse(Duplicates > 1, "Mixed", as.character(Class))) %>%
    group_by(VC, Class) %>%
    count() %>%
    select(-n)

clstr.order <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(Order != "Unassigned") %>%
    group_by(VC, Order) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(Order = ifelse(Duplicates > 1, "Mixed", as.character(Order))) %>%
    group_by(VC, Order) %>%
    count() %>%
    select(-n)

clstr.family <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(Family != "Unassigned") %>%
    group_by(VC, Family) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(Family = ifelse(Duplicates > 1, "Mixed", as.character(Family))) %>%
    group_by(VC, Family) %>%
    count() %>%
    select(-n)

clstr.genus <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(Genus != "Unassigned") %>%
    group_by(VC, Genus) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(Genus = ifelse(Duplicates > 1, "Mixed", as.character(Genus))) %>%
    group_by(VC, Genus) %>%
    count() %>%
    select(-n)



# combine
clstr.master <- clstr.source %>%
    left_join(clstr.phyla, by = "VC") %>%
    left_join(clstr.class, by = "VC") %>%
    left_join(clstr.order, by = "VC") %>%
    left_join(clstr.family, by = "VC") %>%
    left_join(clstr.genus, by = "VC") %>%
    mutate(Phylum = ifelse(is.na(Phylum), "Unassigned", Phylum)) %>%
    mutate(Class = ifelse(is.na(Class), "Unassigned", Class)) %>%
    mutate(Order = ifelse(is.na(Order), "Unassigned", Order)) %>%
    mutate(Family = ifelse(is.na(Family), "Unassigned", Family)) %>%
    mutate(Genus = ifelse(is.na(Genus), "Unassigned", Genus))


# make a taxonomy table for phyloseq
genome.AV <- genome.ov[which(genome.ov$Source == "AV"), ]

genome.AV <- genome.AV[-c(1, 3:5, 19:59)]

# join the taxomy back to the main cluster dataset
genome.AV.tax <- left_join(genome.AV, clstr.master, by = "VC")

# edit taxonomy a bit more for use in plots
genome.AV.tax <- genome.AV.tax %>%
    mutate(ClusterStatus = ifelse(is.na(Size.x), "Unclustered",
        "Clustered")) %>%
    mutate(VCStatus = ifelse(VC == "", "Unclustered", VC)) %>%
    mutate(Phylum = ifelse(is.na(Phylum), "Unassigned", Phylum)) %>%
    mutate(Class = ifelse(is.na(Class), "Unassigned", Class)) %>%
    mutate(Order = ifelse(is.na(Order), "Unassigned", Order)) %>%
    mutate(Family = ifelse(is.na(Family), "Unassigned", Family)) %>%
    mutate(Genus = ifelse(is.na(Genus), "Unassigned", Genus))

genome.AV.justtax <- genome.AV.tax[-c(2, 4:17)]

# vcontact2 also makes a network of protein similarity it
# uses this to make VCs we can use this to plot viruses and
# color by taxonomy or host or genome quality first we need
# to import network info into R and save as R-readable
# files get network information - too big these steps done
# on cluster

# nodes <- ggnet2(ntwk[,-3], mode = 'fruchtermanreingold',
# layout.par = list(list=(niter=2000)))$data %>%
# rename('Genome' = 'label') edges <- ntwk %>% mutate(Pair
# = paste(OTU1, OTU2, sep = '.')) %>% gather(key =
# 'Member', value = 'Genome', -Pair, -Score) %>%
# inner_join(nodes, by = 'Genome') #save files for later
# saveRDS(nodes,
# 'results/processed_data/ntwk_nodes_ictv.RDS')
# saveRDS(edges,
# 'results/processed_data/ntwk_edges_ictv.RDS')

nodes <- readRDS("results/processed_data/ntwk_nodes_ictv.RDS")
edges <- readRDS("results/processed_data/ntwk_edges_ictv.RDS")


# save all the files we made above
saveRDS(genome.ov, "results/processed_data/genome_vc_master_ictv.RDS")
write_csv(genome.ov, "results/processed_data/genome_vc_master_ictv.csv")

saveRDS(clstr.master, "results/processed_data/cluster_vc_master_ictv.RDS")
write_csv(clstr.master, "results/processed_data/cluster_vc_master_ictv.csv")

saveRDS(genome.AV.justtax, "results/processed_data/tax_vc_phyloseq_ictv.RDS")
write_csv(genome.AV.justtax, "results/processed_data/tax_vc_phyloseq_ictv.csv")


## ICTV name conversions
clstr.phyla.ictv <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(ICTV_Phylum != "Unassigned") %>%
    group_by(VC, ICTV_Phylum) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(ICTV_Phylum = ifelse(Duplicates > 1, "Mixed", as.character(ICTV_Phylum))) %>%
    group_by(VC, ICTV_Phylum) %>%
    count() %>%
    select(-n)

clstr.class.ictv <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(ICTV_Class != "Unassigned") %>%
    group_by(VC, ICTV_Class) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(ICTV_Class = ifelse(Duplicates > 1, "Mixed", as.character(ICTV_Class))) %>%
    group_by(VC, ICTV_Class) %>%
    count() %>%
    select(-n)

clstr.order.ictv <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(ICTV_Order != "Unassigned") %>%
    group_by(VC, ICTV_Order) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(ICTV_Order = ifelse(Duplicates > 1, "Mixed", as.character(ICTV_Order))) %>%
    group_by(VC, ICTV_Order) %>%
    count() %>%
    select(-n)

clstr.family.ictv <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(ICTV_Family != "Unassigned") %>%
    group_by(VC, ICTV_Family) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(ICTV_Family = ifelse(Duplicates > 1, "Mixed", as.character(ICTV_Family))) %>%
    group_by(VC, ICTV_Family) %>%
    count() %>%
    select(-n)

clstr.genus.ictv <- genome.ov %>%
    filter(VC.Status == "Clustered") %>%
    filter(VC != "nan") %>%
    filter(ICTV_Genus != "Unassigned") %>%
    group_by(VC, ICTV_Genus) %>%
    count() %>%
    group_by(VC) %>%
    mutate(Duplicates = n()) %>%
    mutate(ICTV_Genus = ifelse(Duplicates > 1, "Mixed", as.character(ICTV_Genus))) %>%
    group_by(VC, ICTV_Genus) %>%
    count() %>%
    select(-n)


# combine
clstr.master.ictv <- clstr.source %>%
    left_join(clstr.phyla.ictv, by = "VC") %>%
    left_join(clstr.class.ictv, by = "VC") %>%
    left_join(clstr.order.ictv, by = "VC") %>%
    left_join(clstr.family.ictv, by = "VC") %>%
    left_join(clstr.genus.ictv, by = "VC") %>%
    mutate(ICTV_Phylum = ifelse(is.na(ICTV_Phylum), "Unassigned",
        ICTV_Phylum)) %>%
    mutate(ICTV_Class = ifelse(is.na(ICTV_Class), "Unassigned",
        ICTV_Class)) %>%
    mutate(ICTV_Order = ifelse(is.na(ICTV_Order), "Unassigned",
        ICTV_Order)) %>%
    mutate(ICTV_Family = ifelse(is.na(ICTV_Family), "Unassigned",
        ICTV_Family)) %>%
    mutate(ICTV_Genus = ifelse(is.na(ICTV_Genus), "Unassigned",
        ICTV_Genus))

# make a taxonomy table for phyloseq
genome.AV <- genome.ov[which(genome.ov$Source == "AV"), ]

genome.AV <- genome.AV[-c(1, 3, 4, 5, 19:59)]

genome.AV.tax.ictv <- left_join(genome.AV, clstr.master.ictv,
    by = "VC")

genome.AV.tax.ictv <- genome.AV.tax.ictv %>%
    mutate(ClusterStatus = ifelse(is.na(Size.x), "Unclustered",
        "Clustered")) %>%
    mutate(VCStatus = ifelse(VC == "", "Unclustered", VC)) %>%
    mutate(ICTV_Phylum = ifelse(is.na(ICTV_Phylum), "Unassigned",
        ICTV_Phylum)) %>%
    mutate(ICTV_Class = ifelse(is.na(ICTV_Class), "Unassigned",
        ICTV_Class)) %>%
    mutate(ICTV_Order = ifelse(is.na(ICTV_Order), "Unassigned",
        ICTV_Order)) %>%
    mutate(ICTV_Family = ifelse(is.na(ICTV_Family), "Unassigned",
        ICTV_Family)) %>%
    mutate(ICTV_Genus = ifelse(is.na(ICTV_Genus), "Unassigned",
        ICTV_Genus))

genome.AV.justtax.ictv <- genome.AV.tax.ictv[-c(2:17)]

saveRDS(genome.AV.justtax.ictv, "results/processed_data/tax_vc_phyloseq_ictv_conversion.RDS")
```

# Importing vOTU RA into phyloseq and re-naming vOTUS

``` r
# load data
tpmean.75 <- readRDS("results/processed_data/antvir_tpm_75_mtx.RDS")
meta <- vroom("data/metadata_all_samples.tsv")
tax <- readRDS("results/processed_data/tax_vc_phyloseq_ictv.RDS")

rownames(meta) <- meta$Sample_name

rownames(tpmean.75) <- word(rownames(tpmean.75))
tpmean.75 <- tpmean.75 %>%
    arrange(rownames(tpmean.75))

# remove taxa that don't meet the coverage thresholds
`%!in%` <- Negate(`%in%`)
no.cov.tax <- tax[which(tax$Genome %!in% rownames(tpmean.75)),
    ]
summary(as.factor(no.cov.tax$ClusterStatus))
```

    ##   Clustered Unclustered 
    ##        2977        1094

``` r
# Clustered Unclustered 2977 1094
length(unique(no.cov.tax$Genome))
```

    ## [1] 4071

``` r
# 4071

tax <- tax[which(tax$Genome %in% rownames(tpmean.75)), ]
tax <- tax %>%
    arrange(Genome)

rownames(tax) <- tax$Genome

tax <- as.matrix(tax[-c(1)])

# import into phyloseq
otu_tab = otu_table(tpmean.75, taxa_are_rows = TRUE)
mapping_file = sample_data(meta)
tax_file = tax_table(tax)

rownames(mapping_file) <- meta$Sample_name3

# rename vOTUs and save conversions
vOTU.names <- as.data.frame(rownames(otu_tab))
vOTU.names$votu.id <- paste0("vOTU", 1:nrow(vOTU.names))

# write.csv(vOTU.names,
# 'results/processed_data/vOTU.names.csv')

rownames(otu_tab) <- paste0("vOTU", 1:nrow(otu_tab))
rownames(tax_file) <- paste0("vOTU", 1:nrow(tax_file))
```

# CheckV and CD-HIT processing

``` r
### CheckV

# Connect CheckV results for each viral sequence back to
# vOTUs and to VCs (clustered levels based on similarity)

# combine individual CheckV files and connect to 'Genome'
# name used in VContact2 results
CV_data_path <- "data/checkv_results/"  # path to the data
CV_files <- dir(CV_data_path, pattern = "*.tsv")  # get file names


CV_data <- data_frame(filename = CV_files) %>%
    mutate(file_contents = map(filename, ~read_tsv(file.path(CV_data_path,
        .))))

# CV_data

CV_data_un <- unnest(CV_data, cols = c(file_contents))  #turns list of files into one file

# Need to fix names so can combine with other results later
CV_data_un <- CV_data_un %>%
    mutate_at("filename", str_replace, ".checkv.quality_summary.tsv",
        "")

CheckV_results <- CV_data_un %>%
    mutate(Genome = paste0(filename, ".", contig_id))

saveRDS(CheckV_results, "results/processed_data/checkv_results.RDS")
write_csv(CheckV_results, "results/processed_data/checkv_results.csv")

### CD-HIT

# Get cluster membership information from cd-hit, to
# connect to vOTUs
cdhit <- read_tsv("data/cdhit/cdhit_cluster_info.txt")

checkv_4_cdhit <- CheckV_results %>%
    mutate(Genome = ifelse(Genome %in% cdhit$id, Genome, paste0(Genome,
        "_1")))

checkv_cdhit <- left_join(checkv_4_cdhit, cdhit, by = c(Genome = "id"))

saveRDS(checkv_cdhit, "results/processed_data/checkv_cdhit_results.RDS")
write_csv(checkv_cdhit, "results/processed_data/checkv_cdhit_results.csv")

# subset to get rep seqs for each cluster, combine with
# taxonomy VC and then remove genome ID and extra info then
# combine back based on clstr # to get taxonomy for each
# vOTU

checkv_cdhit_rep <- checkv_cdhit %>%
    filter(clstr_rep == 1)

checkv_cdhit_rep_tax <- left_join(checkv_cdhit_rep, genome.AV.justtax)
checkv_cdhit_rep_tax <- left_join(checkv_cdhit_rep_tax, genome.AV.justtax.ictv)

# we need to add in vOTU IDs for those vOTUs that had
# insufficient coverage/counts for inclusion earlier note
# the naming here is semi-manual as need to know how many
# vOTUs prior (that is where 72915 came from) could change
# to length of something here instead in future

# vOTU.names

no.cov.vOTU.names <- tibble(no.cov.tax)
no.cov.vOTU.names$`rownames(otu_tab)` <- no.cov.tax$Genome
no.cov.vOTU.names$votu.id <- paste0("vOTU", 72915 + 1:nrow(no.cov.vOTU.names))
no.cov.vOTU.names <- no.cov.vOTU.names[(10:11)]

vOTUs.ids.for.join <- full_join(vOTU.names, no.cov.vOTU.names)

checkv_cdhit_rep_tax_votu <- left_join(checkv_cdhit_rep_tax,
    vOTUs.ids.for.join, by = c(Genome = "rownames(otu_tab)"))

checkv_cdhit_rep_tax_cltr <- checkv_cdhit_rep_tax_votu[-c(1:16,
    18:22)]

checkv_cdhit_tax <- left_join(checkv_cdhit, checkv_cdhit_rep_tax_cltr)

saveRDS(checkv_cdhit_tax, "results/processed_data/checkv_cdhit_results_with_tax_ictv.RDS")
write_csv(checkv_cdhit_tax, "results/processed_data/checkv_cdhit_results_with_tax_ictv.csv")
```

# DRAM-v Processing

``` r
## note DRAM files have different viral contig names and
## need a little love to connect back to checkv/vcontact2
## results

### DRAM-v vMAG stats

# Combine dram distill quality info
DD_data_path <- "data/dramv_results/"  # path to the data
DMAG_files <- dir(DD_data_path, pattern = "*.vMAG_stats.tsv")  # get file names


DMAG_data <- data_frame(filename = DMAG_files) %>%
    mutate(file_contents = map(filename, ~read_tsv(file.path(DD_data_path,
        .))))

DMAG_data_un <- unnest(DMAG_data, cols = c(file_contents))

DMAG_data_un <- DMAG_data_un %>%
    mutate_at("filename", str_replace, "_dramv-distill.vMAG_stats.tsv",
        "")

# fix viral genome names
DMAG_data_un <- DMAG_data_un %>%
    mutate_at("...1", str_replace, "-cat_1", "") %>%
    mutate_at("...1", str_replace, "-cat_2", "") %>%
    mutate_at("...1", str_replace, "-cat_3", "") %>%
    mutate_at("...1", str_replace, "__", "||") %>%
    mutate(Genome = paste0(filename, ".", ...1))

saveRDS(DMAG_data_un[-c(2)], "results/processed_data/dram_vmag_stats.RDS")
write_csv(DMAG_data_un[-c(2)], "results/processed_data/dram_vmag_stats.csv")

dram_vmag_stats <- left_join(checkv_cdhit_tax, DMAG_data_un[-c(1:2)])

saveRDS(dram_vmag_stats, "results/processed_data/checkv_dram_vmag_stats_with_tax_ictv.RDS")
write_csv(dram_vmag_stats, "results/processed_data/checkv_dram_vmag_stats_with_tax_ictv.csv")


### DRAM-v AMG predictions

# combine DRAM-v AMG gene predictions

DD_data_path <- "data/dramv_results/"  # path to the data
DAMG_files <- dir(DD_data_path, pattern = "*.amg_summary.tsv")  # get file names

DAMG_data <- data_frame(filename = DAMG_files) %>%
    mutate(file_contents = map(filename, ~read_tsv(file.path(DD_data_path,
        .)) %>%
        mutate(auxiliary_score = as.character(auxiliary_score))))


DAMG_data_un <- unnest(DAMG_data, cols = c(file_contents))

DAMG_data_un <- DAMG_data_un %>%
    mutate_at("filename", str_replace, "_dramv-distill.amg_summary.tsv",
        "")

# fixing viral genome names
DAMG_data_un <- DAMG_data_un %>%
    mutate_at("scaffold", str_replace, "-cat_1", "") %>%
    mutate_at("scaffold", str_replace, "-cat_2", "") %>%
    mutate_at("scaffold", str_replace, "-cat_3", "") %>%
    mutate_at("scaffold", str_replace, "__", "||") %>%
    mutate(Genome = paste0(filename, ".", scaffold))

saveRDS(DAMG_data_un, "results/processed_data/dram_amg_results.RDS")
write_csv(DAMG_data_un, "results/processed_data/dram_amg_results.csv")
```

# VirSorter2 predicted viral categories

``` r
### VirSorter boundaries

# the authors say these scores aren't definitive, but still
# can provide us with some estimates of types of viruses we
# may have (though we expect mostly dsDNA phage)

VS_data_path <- "data/virsorter_results/"  # path to the data
VS_files <- dir(VS_data_path, pattern = "*.final-viral-score.tsv")  # get file names


VS_data <- data_frame(filename = VS_files) %>%
    mutate(file_contents = map(filename, ~read_tsv(file.path(VS_data_path,
        .))))


VS_data_un <- unnest(VS_data, cols = c(file_contents))

VS_data_un <- VS_data_un %>%
    mutate_at("filename", str_replace, ".final-viral-score.tsv",
        "")

VS_data_un <- VS_data_un %>%
    mutate(Genome = paste0(filename, ".", seqname))

saveRDS(VS_data_un, "results/processed_data/virsorter_score_results.RDS")
write_csv(VS_data_un, "results/processed_data/virsorter_score_results.csv")
```

# Host Prediction

``` r
# load in host prediction from blast
host <- vroom("data/ncbi_host_prediction/bactarch.antvirus.blastn.tsv",
    col_names = c("qseqid", "sseqid", "pident", "bitscore", "evalue",
        "length", "sgi", "sacc", "sallseqid", "staxids", "sscinames",
        "stitle"))

# there are some length / coverage and score requirements
# for a host prediction to be made setting requirements
# length of alignment
minlen = 2000
minpercid = 70
minbit = 50
mine = 0.001  #we set evalue at 1e-25 so not an issue here, they all should be lower than this

# now filter host pred results by the above reqs
host.filt <- host[which(host$bitscore >= minbit), ]
host.filt <- host.filt[which(host.filt$length >= minlen), ]
host.filt <- host.filt[which(host.filt$pident >= minpercid),
    ]

# take the dataframe & group it by the grouping variable
# (query genome) and take the top five hits
host.filt.top5 <- host.filt %>%
    group_by(qseqid) %>%
    slice(1:5)

# run on command line (not in R) cut -f 2
# bactarch.antvirus.blastn.tsv | sort | uniq > acc.txt

# cat acc.txt | epost -db nuccore | esummary -db nuccore |
# xtract -pattern DocumentSummary -element Caption,TaxId >
# taxa.txt

# conda activate ete3 cut -f2 taxa.txt | ete3 ncbiquery
# --info > taxonomy.txt

# load in taxonomy string (taxonomy.txt) and the link back
# to accession # (taxa.txt)
tax2name <- vroom("data/ncbi_host_prediction/taxonomy.txt")
acc2tax <- vroom("data/ncbi_host_prediction/taxa.txt", col_names = FALSE)

# note this is imperfect at 'species' level, seems OK at
# some higher levels, BLAST should also be accurate to
# genus level anyway so should not be an issue for most
# host predictions

# fix names of different levels
tax2name <- tax2name %>%
    mutate_at("Named Lineage", str_replace, "Bacteroidetes/Chlorobi group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Terrabacteria group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "delta/epsilon subdivisions,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Cystobacterineae,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Sorangiineae,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Roseiflexineae,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Nannocystineae,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "PVC group,", "") %>%
    mutate_at("Named Lineage", str_replace, "FCB group,", "") %>%
    mutate_at("Named Lineage", str_replace, "Cyanobacteria/Melainabacteria group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Klebsiella/Raoultella group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Chromobacterium group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Erythrobacter/Porphyrobacter group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Massilia group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Sinorhizobium/Ensifer group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Brucella/Ochrobactrum group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Rhizobium/Agrobacterium group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Azotobacter group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Alteromonas/Salinimonas group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Chryseobacterium group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "TACK group,", "") %>%
    mutate_at("Named Lineage", str_replace, "Chlamydia/Chlamydophila group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Rhizobium/Agrobacterium group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Sinorhizobium/Ensifer group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Brucella/Ochrobactrum group,",
        "") %>%
    mutate_at("Named Lineage", str_replace, "Klebsiella/Raoultella group,",
        "")

# note this is imperfect at 'species' level and may still
# be a bit dirty
tax2name.split <- tax2name %>%
    separate("Named Lineage", into = c("root", "domain", "superkingdom",
        "phylum", "class", "order", "family", "genus", "species"),
        sep = ",")

# join back to accession information
tax_acc <- left_join(acc2tax, tax2name.split, by = c(X2 = "# Taxid"))

host.filt.top5 <- host.filt.top5 %>%
    mutate_at("sacc", str_replace, "\\d$", "") %>%
    mutate_at("sacc", str_replace, "\\.", "")

host.filt.top5.taxa <- left_join(host.filt.top5, tax_acc, by = c(sacc = "X1"))

# collect how many different taxonomic groups in the top 5
# hits per each viral genome
host.filt.top5.taxa.info <- host.filt.top5.taxa %>%
    summarise(Pd = n_distinct(phylum), Cd = n_distinct(class),
        Od = n_distinct(order), Fd = n_distinct(family), Gd = n_distinct(genus),
        Sd = n_distinct(species))

# make a new df
virus.names <- unique(host.filt.top5.taxa.info$qseqid)
host.gen <- as.data.frame(virus.names)
rm(virus.names)

# take the dataframe & group it by the grouping variable
# (genome) take only the top hit now
host.filt.taxa.top1 <- host.filt.top5.taxa %>%
    group_by(qseqid) %>%
    slice(1)

# if only one taxonomic classification at each level in
# host.filt.top5.taxa.info, then we take the phylum of the
# top hit from host.filt.taxa.top1, otherwise we assign it
# as 'Mixed' meaning that multiple taxonomic groupings were
# matches at that level for that genome host.gen.tax <-
# host.gen %>% mutate(Phylum =
# ifelse(host.filt.top5.taxa.info[host.filt.top5.taxa.info$qseqid
# == virus.names,]$Pd == 1,
# host.filt.taxa.top1[host.filt.taxa.top1$qseqid ==
# virus.names,]$phylum, 'Mixed')) %>% mutate(Class =
# ifelse(host.filt.top5.taxa.info[host.filt.top5.taxa.info$qseqid
# == virus.names,]$Cd == 1,
# host.filt.taxa.top1[host.filt.taxa.top1$qseqid ==
# virus.names,]$class, 'Mixed')) %>% mutate(Order =
# ifelse(host.filt.top5.taxa.info[host.filt.top5.taxa.info$qseqid
# == virus.names,]$Od == 1,
# host.filt.taxa.top1[host.filt.taxa.top1$qseqid ==
# virus.names,]$order, 'Mixed')) %>% mutate(Family =
# ifelse(host.filt.top5.taxa.info[host.filt.top5.taxa.info$qseqid
# == virus.names,]$Fd == 1,
# host.filt.taxa.top1[host.filt.taxa.top1$qseqid ==
# virus.names,]$family, 'Mixed')) %>% mutate(Genus =
# ifelse(host.filt.top5.taxa.info[host.filt.top5.taxa.info$qseqid
# == virus.names,]$Gd == 1,
# host.filt.taxa.top1[host.filt.taxa.top1$qseqid ==
# virus.names,]$genus, 'Mixed'))

# save write_csv(host.gen.tax,
# 'results/processed_data/viral.host.taxonomy.csv')
# saveRDS(host.gen.tax,
# 'results/processed_data/host.gen.tax.RDS')

host.gen.tax <- readRDS("results/processed_data/host.gen.tax.RDS")
```

``` r
# load in host prediction from blast
host.mag <- vroom("data/mag_host_prediction/antmag.blastn.tab",
    col_names = c("qseqid", "sseqid", "pident", "bitscore", "evalue",
        "length", "sgi", "sacc", "sallseqid", "staxids", "sscinames",
        "stitle"))
```

    ## Rows: 4141021 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (6): qseqid, sseqid, sacc, sallseqid, sscinames, stitle
    ## dbl (6): pident, bitscore, evalue, length, sgi, staxids
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
host.mag.tax <- vroom("data/mag_host_prediction/Taxonomy_summary.csv")
```

    ## New names:
    ## • `` -> `...8`

    ## Warning: One or more parsing issues, call `problems()` on your data frame for details,
    ## e.g.:
    ##   dat <- vroom(...)
    ##   problems(dat)

    ## Rows: 2278 Columns: 8
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (7): user_genome, Phylum, Class, Order, Family, Genus, Species
    ## lgl (1): ...8
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
host.mag.tax <- host.mag.tax[-c(8)]

# there are some length / coverage and score requirements
# for a host prediction to be made setting requirements
# length of alignment
minlen = 2000
minpercid = 70
minbit = 50
mine = 0.001  #we set evalue at 1e-25 so not an issue here, they all should be lower than this

# now filter host pred results by the above reqs
host.mag.filt <- host.mag[which(host.mag$bitscore >= minbit),
    ]
host.mag.filt <- host.mag.filt[which(host.mag.filt$length >=
    minlen), ]
host.mag.filt <- host.mag.filt[which(host.mag.filt$pident >=
    minpercid), ]

# take the dataframe & group it by the grouping variable
# (query genome) and take the top five hits
host.mag.filt.top5 <- host.mag.filt %>%
    group_by(qseqid) %>%
    slice(1:5)


# note this is imperfect at 'species' level, seems OK at
# some higher levels, BLAST should also be accurate to
# genus level anyway so should not be an issue for most
# host predictions

# fix names of different levels with NAs
host.mag.tax[is.na(host.mag.tax)] <- "Unclassified"

host.mag.tax <- host.mag.tax %>%
    mutate(user_genome = str_remove(user_genome, ".fa"))

host.mag.filt.top5 <- host.mag.filt.top5 %>%
    mutate_at("sacc", str_replace, "_NODE.*", "") %>%
    mutate_at("sacc", str_replace, "_scaff.*", "")

host.mag.filt.top5.taxa <- left_join(host.mag.filt.top5, host.mag.tax,
    by = c(sacc = "user_genome"))

# collect how many different taxonomic groups in the top 5
# hits per each viral genome
host.mag.filt.top5.taxa.info <- host.mag.filt.top5.taxa %>%
    summarise(Pd = n_distinct(Phylum), Cd = n_distinct(Class),
        Od = n_distinct(Order), Fd = n_distinct(Family), Gd = n_distinct(Genus),
        Sd = n_distinct(Species))

# make a new df
virus.names.mag <- unique(host.mag.filt.top5.taxa.info$qseqid)
host.gen.mag <- as.data.frame(virus.names.mag)
rm(virus.names.mag)

# take the dataframe & group it by the grouping variable
# (genome) take only the top hit now
host.mag.filt.taxa.top1 <- host.mag.filt.top5.taxa %>%
    group_by(qseqid) %>%
    slice(1)

# if only one taxonomic classification at each level in
# host.filt.top5.taxa.info, then we take the phylum of the
# top hit from host.filt.taxa.top1, otherwise we assign it
# as 'Mixed' #tomeaning that multiple taxonomic groupings
# were matches at that level for that genome
# host.gen.tax.mag <- host.gen.mag %>% mutate(Phylum =
# ifelse(host.mag.filt.top5.taxa.info[host.mag.filt.top5.taxa.info$qseqid
# == virus.names,]$Pd == 1,
# host.mag.filt.taxa.top1[host.mag.filt.taxa.top1$qseqid ==
# virus.names,]$Phylum, 'Mixed')) %>% mutate(Class =
# ifelse(host.mag.filt.top5.taxa.info[host.mag.filt.top5.taxa.info$qseqid
# == virus.names,]$Cd == 1,
# host.mag.filt.taxa.top1[host.mag.filt.taxa.top1$qseqid ==
# virus.names,]$Class, 'Mixed')) %>% mutate(Order =
# ifelse(host.mag.filt.top5.taxa.info[host.mag.filt.top5.taxa.info$qseqid
# == virus.names,]$Od == 1,
# host.mag.filt.taxa.top1[host.mag.filt.taxa.top1$qseqid ==
# virus.names,]$Order, 'Mixed')) %>% mutate(Family =
# ifelse(host.mag.filt.top5.taxa.info[host.mag.filt.top5.taxa.info$qseqid
# == virus.names,]$Fd == 1,
# host.mag.filt.taxa.top1[host.mag.filt.taxa.top1$qseqid ==
# virus.names,]$Family, 'Mixed')) %>% mutate(Genus =
# ifelse(host.mag.filt.top5.taxa.info[host.mag.filt.top5.taxa.info$qseqid
# == virus.names,]$Gd == 1,
# host.mag.filt.taxa.top1[host.mag.filt.taxa.top1$qseqid ==
# virus.names,]$Genus, 'Mixed'))

# remove Mars-44.bin.73

# saveRDS(host.gen.tax.mag,
# 'results/processed_data/host.gen.tax.mag.RDS')
# write_csv(host.gen.tax.mag,
# 'results/processed_data/viral.host.taxonomy.mag.csv')

host.gen.tax.mag <- readRDS("results/processed_data/host.gen.tax.mag.RDS")

host.gen.tax.mag <- left_join(host.gen.tax.mag, select(host.mag.filt.taxa.top1,
    qseqid, sacc), by = c(virus.names.mag = "qseqid"))

host.gen.tax.mag <- host.gen.tax.mag[host.gen.tax.mag$sacc !=
    "Mars-44.bin.73", ]
```

# Networks

``` r
# load in files using vroom
sampledata <- vroom("data/metadata_all_samples.tsv")
```

    ## New names:
    ## Rows: 191 Columns: 21
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (16): SampleID, Sample_name, Sample_name2, Sample_name3, Sequencing, Sit... dbl
    ## (5): num, Latitude, Longitude, Fossil, Year of collection
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `Site.name` -> `Site.name...7`
    ## • `Site.name` -> `Site.name...8`

``` r
checkv_cdhit_tax <- vroom("results/processed_data/checkv_cdhit_results_with_tax_ictv.csv")
```

    ## Rows: 101085 Columns: 36
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (24): filename, contig_id, provirus, checkv_quality, miuvig_quality, com...
    ## dbl (12): contig_length, proviral_length, gene_count, viral_genes, host_gene...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
rep.gene.info <- vroom("results/processed_data/dram_vmag_stats.csv")
```

    ## Rows: 95939 Columns: 16
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr  (2): filename, Genome
    ## dbl (11): VIRSorter category, Gene count, Strand switches, potential AMG cou...
    ## lgl  (3): Circular, Prophage, Transposase present
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# we need a column with values that are shared between both
# files but! we don't have one yet - as filename has
# SampleID.contigs we need to remove the '.contigs' from
# each value in the filename column in order for it to
# match SampleID in sampledata
checkv_cdhit_tax2 <- checkv_cdhit_tax %>%
    mutate_at("filename", str_replace, ".contigs", "")

# Mars samples are '_' not '-'
sampledata2 <- sampledata %>%
    mutate_at("SampleID", str_replace, "Mars_", "Mars-")

# join all columns in both files using SampleID = filename
# which should now both have the same values we will use
# this for some of our cluster filtering below and later to
# plot
combo <- full_join(sampledata2, checkv_cdhit_tax2, by = c(SampleID = "filename"))

# these take a while to load the images

clstrd.nodes <- nodes %>%
    left_join(genome.ov, by = "Genome") %>%
    filter(VC.Status == "Clustered")

# plot by family
genome.ov.AV <- left_join(genome.ov, checkv_cdhit_tax)
```

    ## Joining, by = c("Genome", "VC.Status", "Class", "Phylum", "ICTV_Genus",
    ## "ICTV_Family", "ICTV_Order", "ICTV_Class", "ICTV_Phylum", "Order", "Family",
    ## "Genus")

``` r
# plot by fam - only VCs that clustered with reference data
vir.fam.clust <- genome.ov %>%
    filter(Source == "AV") %>%
    select(Genome, VC, VC.Status) %>%
    inner_join(clstr.master, by = "VC") %>%
    filter(ClstrComp == "both") %>%
    filter(Order != "Unassigned") %>%
    left_join(select(combo, VCStatus, checkv_quality, contig_length),
        by = c(VC = "VCStatus")) %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(is.na(contig_length) | contig_length >= 10000) %>%
    group_by(Order, Family) %>%
    count() %>%
    group_by(Order) %>%
    mutate(nOrd = sum(n)) %>%
    arrange(nOrd, n) %>%
    ungroup() %>%
    mutate(Rank = 1:n())

vir.fam.nodes <- clstrd.nodes %>%
    left_join(select(clstr.master, VC, ClstrComp), by = "VC") %>%
    left_join(combo) %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(is.na(contig_length) | contig_length >= 10000) %>%
    mutate(Family2 = case_when(Source == "AV" ~ "This Study",
        Family %in% vir.fam.clust$Family ~ as.character(Family),
        TRUE ~ "Other")) %>%
    left_join(vir.fam.clust, by = c(Family2 = "Family")) %>%
    filter(Source == "refseq" | ClstrComp == "both")
```

    ## Joining, by = c("Genome", "VC.Status", "Class", "Phylum", "ICTV_Genus",
    ## "ICTV_Family", "ICTV_Order", "ICTV_Class", "ICTV_Phylum", "Order", "Family",
    ## "Genus")

``` r
# Plot
ntwk.vir <- vir.fam.nodes %>%
    ggplot(aes(x, y)) + geom_line(data = filter(edges, Genome %in%
    vir.fam.nodes$Genome), aes(group = Pair), alpha = 0.1, color = "gray25",
    size = 0.5) + geom_point(alpha = 0.8, size = 2, shape = 16,
    aes(color = Family2)) + scale_color_manual(name = "Virus Family",
    values = c(RColorBrewer::brewer.pal(8, "Set2")[c(2:5)], RColorBrewer::brewer.pal(8,
        "Set2")[c(1)], "slateblue4", RColorBrewer::brewer.pal(8,
        "Set2")[c(6)], "gray10", "gray75")) + theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), legend.position = "right") +
    guides(color = guide_legend(override.aes = list(shape = 16,
        size = 4)))

# ntwk.vir




# plot by rock type - only HQ

vir.hab.nodes <- clstrd.nodes %>%
    left_join(select(clstr.master, VC, ClstrComp), by = "VC") %>%
    left_join(select(combo, checkv_quality, contig_length, Rocks_v2,
        VCStatus), by = c(VC = "VCStatus")) %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(is.na(contig_length) | contig_length >= 10000) %>%
    mutate(Rocks_v2 = ifelse(is.na(Rocks_v2), "Unknown", as.character(Rocks_v2)))  #%>%
# filter(Source == 'refseq' | ClstrComp == 'both')

vir.hab.nodes <- vir.hab.nodes[-c(66:67)]


# Plot
ntwk.hab.vir.hq <- vir.hab.nodes %>%
    ggplot(aes(x, y)) + geom_line(data = filter(edges, Genome %in%
    vir.hab.nodes$Genome), aes(group = Pair), alpha = 0.1, color = "gray25",
    size = 0.5) + geom_point(alpha = 0.8, size = 2, shape = 16,
    aes(color = Rocks_v2)) + scale_color_manual(name = "Rock Type",
    values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "gray75",
        "#0072B2", "#D55E00", "#CC79A7")) + theme_minimal() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), legend.position = "right") +
    guides(color = guide_legend(override.aes = list(shape = 16,
        size = 4)))


# ntwk.hab.vir.hq
```

# Generating Plots / Tables

# Figure 1

# Network + Bars

``` r
## rock type ##
meta.sub.rock <- combo %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(contig_length >= 10000) %>%
    filter(ClusterStatus == "Clustered") %>%
    group_by(Rocks_v2, Order) %>%
    mutate(Order = ifelse(Order == "Unassigned", "Unique VC",
        "VC with reference genomes")) %>%
    tally() %>%
    arrange(desc(Rocks_v2)) %>%
    mutate(lab_ypos = n + 0.1 * n + ifelse(n > 40, ifelse(n <
        50, 50, 20), 1)) %>%
    mutate(norm = n/sum(n) * 100)


a_hab_seq <- ggplot(meta.sub.rock, aes(y = Rocks_v2, fill = Rocks_v2,
    color = Rocks_v2, x = n)) + theme_bw() + geom_bar(stat = "identity",
    position = "stack", width = 0.6, orientation = "y") + scale_fill_manual(name = "Habitat",
    values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
        "#D55E00", "#CC79A7", "gray75")) + scale_color_manual(name = "Habitat",
    values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
        "#D55E00", "#CC79A7", "gray75")) + ylab("") + facet_wrap(~Order,
    scales = "free") + xlab("Number of viral sequences") + scale_y_discrete(limits = rev) +
    theme(legend.position = "none")  #+
# geom_text(aes(x = lab_ypos,label = n, group = Rocks_v2),
# fontface = 'bold', size = 4, color = 'black')

## taxonomy ##
tax.meta.vcs <- combo %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(contig_length >= 10000) %>%
    filter(ClusterStatus == "Clustered" & Order != "Unassigned") %>%
    mutate(Family = ifelse(Family == "Mixed", "Unclassified",
        as.character(Family))) %>%
    group_by(Family) %>%
    tally() %>%
    arrange(desc(n)) %>%
    group_by(Family) %>%
    mutate(n2 = sum(n)) %>%
    mutate(lab_ypos = n2 + 3)

e_tax <- ggplot(tax.meta.vcs, aes(x = Family, fill = Family,
    y = n)) + theme_bw() + geom_bar(stat = "identity", position = "stack",
    width = 0.6) + scale_fill_manual(values = c(RColorBrewer::brewer.pal(8,
    "Set2")[c(2:4)], RColorBrewer::brewer.pal(8, "Set2")[c(1)],
    "slateblue4", RColorBrewer::brewer.pal(8, "Set2")[c(6)],
    "gray75")) + ylab("Number of viral sequences") + xlab("") +
    theme(legend.position = "none") + coord_flip() + scale_x_discrete(limits = rev)  # +
# geom_text(aes(y = lab_ypos,label = n2, group = Family),
# fontface = 'bold', size = 4, color = 'black')



(ntwk.hab.vir.hq + a_hab_seq)/(ntwk.vir + e_tax) + plot_annotation(tag_levels = "A")
```

![](09_Antarctic_RMarkDown_files/figure-gfm/fig1-1.png)<!-- -->

``` r
ggsave(filename = "plots/exploratory/ant_network_plus_VC.png",
    plot = last_plot(), device = "png", width = 14, height = 9,
    dpi = 300)
```

# Host/AMG plots

``` r
## Host plot ##

host.gen.tax <- readRDS("results/processed_data/host.gen.tax.RDS")
meta <- vroom("data/metadata_all_samples.tsv")
```

    ## New names:
    ## Rows: 191 Columns: 21
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "\t" chr
    ## (16): SampleID, Sample_name, Sample_name2, Sample_name3, Sequencing, Sit... dbl
    ## (5): num, Latitude, Longitude, Fossil, Year of collection
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `Site.name` -> `Site.name...7`
    ## • `Site.name` -> `Site.name...8`

``` r
# join and rename
host.vots <- left_join(host.gen.tax, vOTU.names, by = c(virus.names = "rownames(otu_tab)"))
host.vots <- host.vots[-c(1)]
host.vots$HostPhylum <- host.vots$Phylum
host.vots$HostClass <- host.vots$Class
host.vots$HostOrder <- host.vots$Order
host.vots$HostFamily <- host.vots$Family
host.vots$HostGenus <- host.vots$Genus

host.vots <- host.vots[-c(1:5)]


host.gen.tax.meta <- left_join(combo, host.vots, by = c(votu.id = "votu.id"))

# fix labels
host.gen.tax.meta <- host.gen.tax.meta %>%
    mutate(HostPhylum = fct_explicit_na(HostPhylum, na_level = "No host prediction"),
        HostClass = fct_explicit_na(HostClass, na_level = "No host prediction"),
        HostOrder = fct_explicit_na(HostOrder, na_level = "No host prediction"),
        HostFamily = fct_explicit_na(HostFamily, na_level = "No host prediction"),
        HostGenus = fct_explicit_na(HostGenus, na_level = "No host prediction")) %>%
    mutate(HostPhylum = ifelse(HostPhylum == "Mixed", "No host prediction",
        as.character(HostPhylum)), HostClass = ifelse(HostClass ==
        "Mixed", "No host prediction", as.character(HostClass)),
        HostOrder = ifelse(HostOrder == "Mixed", "No host prediction",
            as.character(HostOrder)), HostFamily = ifelse(HostFamily ==
            "Mixed", "No host prediction", as.character(HostFamily)),
        HostGenus = ifelse(HostGenus == "Mixed", "No host prediction",
            as.character(HostGenus)))

host.gen.tax.meta.high <- host.gen.tax.meta %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(contig_length >= 10000)

# summarize data
host.gen.tax.meta.plotdata <- host.gen.tax.meta.high %>%
    filter(HostPhylum != "No host prediction") %>%
    group_by(HostPhylum) %>%
    tally() %>%
    arrange(desc(n)) %>%
    group_by(HostPhylum) %>%
    mutate(n2 = sum(n)) %>%
    mutate(lab_ypos = n2 + 3)


host_plot <- ggplot(host.gen.tax.meta.plotdata, aes(x = HostPhylum,
    fill = HostPhylum, y = n)) + theme_bw() + geom_bar(stat = "identity",
    position = "stack", width = 0.6) + scale_fill_viridis_d(option = "B") +
    ylab("Number of viral sequences") + xlab("") + theme(legend.position = "none") +
    coord_flip() + scale_x_discrete(limits = rev)  #+
# geom_text(aes(y = lab_ypos,label = n2, group =
# HostPhylum), fontface = 'bold', size = 4, color =
# 'black')


## AMG plot ##

DAMG_data_un_HQ <- DAMG_data_un %>%
    filter(Genome %in% host.gen.tax.meta.high$Genome)

# this has all annotations for a gene - so can have
# multiple annotations, we need to address or we will be
# over-counting
DAMG_data_un_HQ_counts <- DAMG_data_un_HQ %>%
    count(gene)

DAMG_data_un_HQ <- left_join(DAMG_data_un_HQ, DAMG_data_un_HQ_counts)
```

    ## Joining, by = "gene"

``` r
DAMG_data_un_HQ.singleAnnot <- DAMG_data_un_HQ %>%
    group_by(gene) %>%
    slice(1)

AMGs <- left_join(DAMG_data_un_HQ.singleAnnot, host.gen.tax.meta.high)
```

    ## Joining, by = "Genome"

``` r
amg.plotdata <- AMGs %>%
    mutate(category = fct_explicit_na(category, na_level = "Unannotated"),
        module = fct_explicit_na(module, na_level = "Unannotated"),
        header = ifelse(category == "Transporters", "Transporters",
            as.character(header)), header = fct_explicit_na(header,
            na_level = "Unannotated"), category = ifelse(category ==
            "carbon utilization (Woodcroft)", "carbon utilization",
            as.character(category)), header = ifelse(header ==
            "sugar utilization (woodcroft)", "carbon utilization",
            as.character(header)), category = ifelse(category ==
            "MISC", "Information Systems", as.character(category)),
        category = ifelse(category == "carbon utilization", "Carbon Utilization",
            as.character(category))) %>%
    group_by(category, header, module) %>%
    tally()

amg_general <- ggplot(amg.plotdata, aes(x = category, fill = category,
    y = n)) + theme_bw() + geom_bar(stat = "identity", position = "stack",
    width = 0.6) + coord_flip() + scale_fill_viridis_d(option = "G",
    direction = -1) + ylab("Number of auxillary metabolic genes") +
    xlab("") + theme(legend.position = "none") + scale_x_discrete(limits = rev)
```

``` r
## Host MAG ##

# join and rename
host.vots.mag <- left_join(host.gen.tax.mag, vOTU.names, by = c(virus.names.mag = "rownames(otu_tab)"))
host.vots.mag <- host.vots.mag[-c(1)]
host.vots.mag$MAGHostPhylum <- host.vots.mag$Phylum
host.vots.mag$MAGHostClass <- host.vots.mag$Class
host.vots.mag$MAGHostOrder <- host.vots.mag$Order
host.vots.mag$MAGHostFamily <- host.vots.mag$Family
host.vots.mag$MAGHostGenus <- host.vots.mag$Genus
host.vots.mag$MAGHostBin <- host.vots.mag$sacc

host.vots.mag <- host.vots.mag[-c(1:6)]


host.gen.tax.mag.meta <- left_join(combo, host.vots.mag, by = c(votu.id = "votu.id"))

# fix labels
host.gen.tax.mag.meta <- host.gen.tax.mag.meta %>%
    mutate(MAGHostPhylum = fct_explicit_na(MAGHostPhylum, na_level = "No host prediction"),
        MAGHostClass = fct_explicit_na(MAGHostClass, na_level = "No host prediction"),
        MAGHostOrder = fct_explicit_na(MAGHostOrder, na_level = "No host prediction"),
        MAGHostFamily = fct_explicit_na(MAGHostFamily, na_level = "No host prediction"),
        MAGHostGenus = fct_explicit_na(MAGHostGenus, na_level = "No host prediction"),
        MAGHostBin = fct_explicit_na(MAGHostBin, na_level = "No host prediction")) %>%
    mutate(MAGHostPhylum = ifelse(MAGHostPhylum == "Mixed", "No host prediction",
        as.character(MAGHostPhylum)), MAGHostClass = ifelse(MAGHostClass ==
        "Mixed", "No host prediction", as.character(MAGHostClass)),
        MAGHostOrder = ifelse(MAGHostOrder == "Mixed", "No host prediction",
            as.character(MAGHostOrder)), MAGHostFamily = ifelse(MAGHostFamily ==
            "Mixed", "No host prediction", as.character(MAGHostFamily)),
        MAGHostGenus = ifelse(MAGHostGenus == "Mixed", "No host prediction",
            as.character(MAGHostGenus)))

host.gen.tax.mag.meta.high <- host.gen.tax.mag.meta %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(contig_length >= 10000)

# summarize data
host.gen.tax.mag.meta.high.plotdata <- host.gen.tax.mag.meta.high %>%
    filter(MAGHostPhylum != "No host prediction") %>%
    group_by(MAGHostPhylum) %>%
    tally() %>%
    arrange(desc(n)) %>%
    group_by(MAGHostPhylum) %>%
    mutate(n2 = sum(n)) %>%
    mutate(lab_ypos = n2 + 3)


host_mag_plot <- ggplot(host.gen.tax.mag.meta.high.plotdata,
    aes(x = MAGHostPhylum, fill = MAGHostPhylum, y = n)) + theme_bw() +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_viridis_d(option = "B") + ylab("Number of viral sequences") +
    xlab("") + theme(legend.position = "none") + coord_flip() +
    scale_x_discrete(limits = rev)  #+
# geom_text(aes(y = lab_ypos,label = n2, group =
# HostPhylum), fontface = 'bold', size = 4, color =
# 'black')

host_mag_plot
```

![](09_Antarctic_RMarkDown_files/figure-gfm/maghost_plot-1.png)<!-- -->

``` r
host.gen.tax.mag.meta.high.plotdata$DB <- "(A) MAG-based prediction"
host.gen.tax.mag.meta.high.plotdata$HostPhylum <- host.gen.tax.mag.meta.high.plotdata$MAGHostPhylum
host.gen.tax.meta.plotdata$DB <- "(B) RefSeq-based prediction"

host.plotdata <- full_join(host.gen.tax.meta.plotdata, host.gen.tax.mag.meta.high.plotdata[-c(1)])
```

    ## Joining, by = c("HostPhylum", "n", "n2", "lab_ypos", "DB")

``` r
host.plotdata <- host.plotdata %>%
    mutate(HostPhylum = ifelse(HostPhylum == "Actinobacteria",
        "Actinobacteriota", as.character(HostPhylum)), HostPhylum = ifelse(HostPhylum ==
        "Bacteroidetes", "Bacteroidota", as.character(HostPhylum)),
        HostPhylum = ifelse(HostPhylum == "Planctomycetes", "Planctomycetota",
            as.character(HostPhylum)), HostPhylum = ifelse(HostPhylum ==
            "Acidobacteria", "Acidobacteriota", as.character(HostPhylum)),
        HostPhylum = ifelse(HostPhylum == "Chloroflexi", "Chloroflexota",
            as.character(HostPhylum)), HostPhylum = ifelse(HostPhylum ==
            "Deinococcus-Thermus", "Deinococcota", as.character(HostPhylum)))

ggplot(host.plotdata, aes(x = HostPhylum, fill = HostPhylum,
    y = n)) + theme_bw() + facet_wrap(~DB) + geom_bar(stat = "identity",
    position = "stack", width = 0.6) + scale_fill_viridis_d(option = "B") +
    ylab("Number of viral sequences") + xlab("") + theme(legend.position = "none") +
    coord_flip() + scale_x_discrete(limits = rev)
```

![](09_Antarctic_RMarkDown_files/figure-gfm/combo_supp_mag-1.png)<!-- -->

``` r
ggsave(filename = "plots/figs1.pdf", plot = last_plot(), device = "pdf",
    width = 6, height = 4, dpi = 300)
ggsave(filename = "plots/figs1.png", plot = last_plot(), device = "png",
    width = 6, height = 4, dpi = 300)
```

``` r
# host vs ncbi

host.combo <- full_join(host.gen.tax.meta.high, host.gen.tax.mag.meta.high)
```

    ## Joining, by = c("num", "SampleID", "Sample_name", "Sample_name2",
    ## "Sample_name3", "Sequencing", "Site.name...7", "Site.name...8", "Sample",
    ## "AreaSample", "Type of rocks", "Rocks_v2", "Latitude", "Longitude", "Sun
    ## exposure", "Elevation (m asl)", "Sea distance (km)", "OriginalYear",
    ## "Continent", "Fossil", "Year of collection", "contig_id", "contig_length",
    ## "provirus", "proviral_length", "gene_count", "viral_genes", "host_genes",
    ## "checkv_quality", "miuvig_quality", "completeness", "completeness_method",
    ## "contamination", "kmer_freq", "warnings", "Genome", "clstr", "clstr_size",
    ## "length", "clstr_rep", "clstr_iden", "clstr_cov", "VC.Status", "Phylum",
    ## "Class", "Order", "Family", "Genus", "ClusterStatus", "VCStatus",
    ## "ICTV_Phylum", "ICTV_Class", "ICTV_Order", "ICTV_Family", "ICTV_Genus",
    ## "votu.id")

``` r
host.combo <- host.combo %>%
    mutate(MAGvNCBI = ifelse(HostPhylum == "No host prediction" &
        MAGHostPhylum == "No host prediction", "No host prediction",
        ifelse(HostPhylum != "No host prediction" & MAGHostPhylum !=
            "No host prediction", "NCBI and MAG based prediction",
            ifelse(MAGHostPhylum != "No host prediction", "MAG-based prediction",
                "NCBI-based prediction"))))

# mag based methods should be superior to NCBI methods
summary(as.factor(host.combo$MAGvNCBI))
```

    ##          MAG-based prediction NCBI and MAG based prediction 
    ##                          5668                           968 
    ##         NCBI-based prediction            No host prediction 
    ##                           139                          8021

``` r
host.combo <- host.combo %>%
    mutate(HostPhylum = ifelse(HostPhylum == "Actinobacteria",
        "Actinobacteriota", as.character(HostPhylum)), HostPhylum = ifelse(HostPhylum ==
        "Bacteroidetes", "Bacteroidota", as.character(HostPhylum)),
        HostPhylum = ifelse(HostPhylum == "Planctomycetes", "Planctomycetota",
            as.character(HostPhylum)), HostPhylum = ifelse(HostPhylum ==
            "Acidobacteria", "Acidobacteriota", as.character(HostPhylum)),
        HostPhylum = ifelse(HostPhylum == "Chloroflexi", "Chloroflexota",
            as.character(HostPhylum)), HostPhylum = ifelse(HostPhylum ==
            "Deinococcus-Thermus", "Deinococcota", as.character(HostPhylum)))

host.combo <- host.combo %>%
    mutate(ComboPhylum = ifelse(MAGvNCBI == "MAG-based prediction",
        as.character(MAGHostPhylum), ifelse(MAGvNCBI == "NCBI-based prediction",
            as.character(HostPhylum), ifelse(MAGvNCBI == "NCBI and MAG based prediction",
                ifelse(HostPhylum == MAGHostPhylum, MAGHostPhylum,
                  "Mixed"), as.character(MAGvNCBI)))), ComboPhylum_MAG = ifelse(MAGvNCBI ==
        "MAG-based prediction", as.character(MAGHostPhylum),
        ifelse(MAGvNCBI == "NCBI-based prediction", as.character(HostPhylum),
            ifelse(MAGvNCBI == "NCBI and MAG based prediction",
                as.character(MAGHostPhylum), as.character(MAGvNCBI)))))


host.combo.plotdata <- host.combo %>%
    filter(ComboPhylum_MAG != "No host prediction") %>%
    group_by(ComboPhylum_MAG) %>%
    tally() %>%
    arrange(desc(n)) %>%
    group_by(ComboPhylum_MAG) %>%
    mutate(n2 = sum(n)) %>%
    mutate(lab_ypos = n2 + 3)

ComboPhylum_MAG_plot <- ggplot(host.combo.plotdata, aes(x = ComboPhylum_MAG,
    fill = ComboPhylum_MAG, y = n)) + theme_bw() + geom_bar(stat = "identity",
    position = "stack", width = 0.6) + scale_fill_viridis_d(option = "B") +
    ylab("Number of viral sequences") + xlab("") + theme(legend.position = "none") +
    coord_flip() + scale_x_discrete(limits = rev)  #+
# geom_text(aes(y = lab_ypos,label = n2, group =
# HostPhylum), fontface = 'bold', size = 4, color =
# 'black')

ComboPhylum_MAG_plot
```

![](09_Antarctic_RMarkDown_files/figure-gfm/hostvncbi-1.png)<!-- -->

``` r
ggsave(filename = "plots/exploratory/host_combo_mag.pdf", plot = last_plot(),
    device = "pdf", width = 6, height = 6, dpi = 300)
```

# Combined Fig 1

``` r
# Figure 1
(ntwk.vir + e_tax)/(ntwk.hab.vir.hq + a_hab_seq)/(ComboPhylum_MAG_plot +
    amg_general) + plot_annotation(tag_levels = "A")
```

![](09_Antarctic_RMarkDown_files/figure-gfm/fig1_plot-1.png)<!-- -->

``` r
ggsave(filename = "plots/fig1.pdf", plot = last_plot(), device = "pdf",
    width = 14, height = 8, dpi = 300)
ggsave(filename = "plots/fig1.png", plot = last_plot(), device = "png",
    width = 14, height = 8, dpi = 300)
```

# Figure 2

# vOTU analyses

``` r
# import vOTU and host data into phyloseq
ps <- phyloseq(otu_tab, tax_file, mapping_file)

## Get relative abundance averages for plotting ##

# For some visualizations - show only genomes with a
# quality score from CheckV

votus.high <- unique(na.omit(checkv_cdhit_tax$votu.id[checkv_cdhit_tax$checkv_quality %in%
    c("High-quality", "Complete", "Medium-quality", "Low-quality")]))
votus.long <- checkv_cdhit_tax$votu.id[checkv_cdhit_tax$contig_length >=
    10000]

ps.HQ <- prune_taxa(votus.high, ps)
ps.HQ.long <- prune_taxa(votus.long, ps.HQ)

ps.RA <- transform_sample_counts(ps.HQ.long, function(x) x/sum(x))

df_ps.RA.filt <- psmelt(ps.RA)


## Site ##

grouped <- group_by(df_ps.RA.filt, Site.name...8, ClusterStatus,
    Phylum, Class, Order, Family, VCStatus, VC.Status, OTU)

avgs_grouped <- summarise(grouped, mean = 100 * mean(Abundance),
    sd = 100 * sd(Abundance))

vOTU_avgs_grouped <- avgs_grouped


vOTU_avgs_grouped <- left_join(vOTU_avgs_grouped, select(clstr.source,
    VC, ClstrComp), by = c(VCStatus = "VC"))

vOTU_avgs_grouped <- vOTU_avgs_grouped %>%
    mutate(ClstrComp = ifelse(is.na(ClstrComp), "none", as.character(ClstrComp))) %>%
    mutate(Family2 = ifelse(Family == "Mixed", "Unclassified",
        as.character(Family))) %>%
    mutate(Family2 = ifelse(Family2 == "Unassigned", ifelse(ClusterStatus ==
        "Clustered", "Unclassified", "Unclustered"), as.character(Family2))) %>%
    mutate(Family3 = ifelse(Family == "Unassigned", ifelse(ClusterStatus ==
        "Clustered", "Unclassified", as.character(VC.Status)),
        as.character(Family2))) %>%
    mutate(Family3 = ifelse(str_detect(as.character(Family3),
        "Overlap"), "Overlap", as.character(Family3))) %>%
    mutate(Family4 = ifelse(Family3 == "Unclassified", ifelse(ClstrComp ==
        "AV", "Unique VC", as.character(Family3)), as.character(Family3)))


vOTU_avgs_grouped$Site.name...8 <- factor(vOTU_avgs_grouped$Site.name...8,
    levels = c("Machu Picchu Base", "Dufayel Island", "Helliwell Hills",
        "Crater Circle", "Mt. Burrow", "Archambault Ridge", "Random Hills",
        "Kay Island", "Timber Peak", "Mt. New Zealand", "Mt. Dickason",
        "Mt. Nansen", "Mt. Keinath", "Anderson Ridge", "Vegetation Island",
        "Mt. Larsen", "Inexpressible Island", "Mt. Billing",
        "Widowmaker Pass", "Trio Nunatak", "Unnamed Nunatak",
        "Harrow Peak", "Mt. McGee", "Ricker Hills", "Mt. Bowen",
        "Pudding Butte", "Starr Nunatak", "Richard Nunatak",
        "Schultz Peak", "Battleship Promontory", "Convoy Range",
        "Mt. Elektra", "Siegfried Peak", "Linnaeus Terrace",
        "Finger Mt.", "University Valley", "Knobhead"))

vOTU_avgs_grouped$Family4 <- factor(vOTU_avgs_grouped$Family4,
    levels = c("Helgolandviridae", "Inoviridae", "Myoviridae",
        "Podoviridae", "Schitoviridae", "Siphoviridae", "Unclassified",
        "Unique VC", "Overlap", "Outlier", "Singleton"))

plot_s = ggplot(vOTU_avgs_grouped, aes(x = Site.name...8, y = (mean),
    fill = Family4)) + geom_bar(stat = "identity", position = "stack")

site_bar_ra <- plot_s + theme_bw() + theme(text = element_text(size = 14)) +
    ylab("Mean Relative Abundance") + xlab("") + guides(fill = guide_legend(title = "Family")) +
    scale_fill_viridis_d()  #+  theme(legend.position = 'none') +
# theme(axis.text.x = element_text(angle = -70, hjust = 0,
# vjust = 0.5))

site_bar_ra + coord_flip()
```

![](09_Antarctic_RMarkDown_files/figure-gfm/morephyloseq-1.png)<!-- -->

``` r
ggsave(filename = "plots/exploratory/art_virus_RA_ictv.png",
    plot = last_plot(), device = "png", width = 6, height = 10,
    dpi = 300)
```

``` r
# get number of vOTUs across three regions shared
vOTU_avgs_grouped.shared <- vOTU_avgs_grouped %>%
    left_join(select(meta, Site.name...8, AreaSample)) %>%
    filter(mean != "0") %>%
    group_by(AreaSample, OTU) %>%
    tally() %>%
    group_by(OTU) %>%
    summarize(n_uniq = length(unique(AreaSample)))
```

    ## Joining, by = "Site.name...8"

``` r
# found in one region
length(vOTU_avgs_grouped.shared$OTU[vOTU_avgs_grouped.shared$n_uniq ==
    "1"])
```

    ## [1] 7759

``` r
# found in two regions
length(vOTU_avgs_grouped.shared$OTU[vOTU_avgs_grouped.shared$n_uniq ==
    "2"])
```

    ## [1] 3066

``` r
# found in three regions
length(vOTU_avgs_grouped.shared$OTU[vOTU_avgs_grouped.shared$n_uniq ==
    "3"])
```

    ## [1] 159

``` r
vOTU.3region.info <- combo[combo$votu.id %in% vOTU_avgs_grouped.shared$OTU[vOTU_avgs_grouped.shared$n_uniq ==
    "3"], ]

vOTU.3region.info <- vOTU.3region.info %>%
    left_join(select(clstr.source, VC, ClstrComp), by = c(VCStatus = "VC")) %>%
    mutate(ClstrComp = ifelse(is.na(ClstrComp), "none", as.character(ClstrComp))) %>%
    mutate(Family2 = ifelse(Family == "Mixed", "Unclassified",
        as.character(Family))) %>%
    mutate(Family2 = ifelse(Family2 == "Unassigned", ifelse(VC.Status ==
        "Clustered", "Unclassified", as.character(VC.Status)),
        as.character(Family2))) %>%
    mutate(Family2 = ifelse(str_detect(as.character(Family2),
        "Overlap"), "Overlap", as.character(Family2))) %>%
    mutate(Family2 = ifelse(Family2 == "Unclassified", ifelse(ClstrComp ==
        "AV", "Unique VC", as.character(Family2)), as.character(Family2))) %>%
    left_join(select(host.gen.tax.meta, Genome, HostPhylum, HostClass,
        HostOrder, HostFamily, HostGenus))
```

    ## Joining, by = "Genome"

``` r
# VC status of overlapping sequences
summary(as.factor(vOTU.3region.info$Family2))
```

    ## Clustered/Singleton             Outlier             Overlap           Singleton 
    ##                   4                  24                 123                  19 
    ##           Unique VC 
    ##                 228

``` r
# VC status, for vOTUs
summary(as.factor(vOTU.3region.info$Family2[vOTU.3region.info$clstr_rep ==
    1]))
```

    ## Clustered/Singleton             Outlier             Overlap           Singleton 
    ##                   2                  16                  21                  14 
    ##           Unique VC 
    ##                 106

``` r
# Hosts status of overlapping sequences
summary(as.factor(vOTU.3region.info$HostPhylum))
```

    ##       Acidobacteria      Actinobacteria       Bacteroidetes Deinococcus-Thermus 
    ##                   5                  16                   3                   1 
    ##  No host prediction      Proteobacteria 
    ##                 259                 114

``` r
# Host phyla, vOTUs
summary(as.factor(vOTU.3region.info$HostPhylum[vOTU.3region.info$clstr_rep ==
    1]))
```

    ##       Acidobacteria      Actinobacteria       Bacteroidetes Deinococcus-Thermus 
    ##                   2                   6                   2                   1 
    ##  No host prediction      Proteobacteria 
    ##                  97                  51

``` r
# vOTUs across sites
vOTU_avgs_grouped.shared.site <- vOTU_avgs_grouped %>%
    filter(mean != "0") %>%
    group_by(Site.name...8, OTU) %>%
    tally() %>%
    group_by(OTU) %>%
    summarize(n_uniq = length(unique(Site.name...8)))

# found in only one site
length(vOTU_avgs_grouped.shared.site$OTU[vOTU_avgs_grouped.shared.site$n_uniq ==
    "1"])
```

    ## [1] 6485

``` r
# found more than one site
length(vOTU_avgs_grouped.shared.site$OTU[vOTU_avgs_grouped.shared.site$n_uniq >
    "1"])
```

    ## [1] 4499

# Beta Diversity

``` r
# euclidean dist of hellinger transformation (hellinger
# distance)
ps_hell <- transform(ps.HQ.long, "hellinger")

ps_hell_ord <- ordinate(physeq = ps_hell, method = "PCoA", distance = "euclidean")

ps_hell@sam_data$Site.name...8 <- factor(ps_hell@sam_data$Site.name...8,
    levels = c("Machu Picchu Base", "Dufayel Island", "Helliwell Hills",
        "Crater Circle", "Mt. Burrow", "Archambault Ridge", "Random Hills",
        "Kay Island", "Timber Peak", "Mt. New Zealand", "Mt. Dickason",
        "Mt. Nansen", "Mt. Keinath", "Anderson Ridge", "Vegetation Island",
        "Mt. Larsen", "Inexpressible Island", "Mt. Billing",
        "Widowmaker Pass", "Trio Nunatak", "Unnamed Nunatak",
        "Harrow Peak", "Mt. McGee", "Ricker Hills", "Mt. Bowen",
        "Pudding Butte", "Starr Nunatak", "Richard Nunatak",
        "Schultz Peak", "Battleship Promontory", "Convoy Range",
        "Mt. Elektra", "Siegfried Peak", "Linnaeus Terrace",
        "Finger Mt.", "University Valley", "Knobhead"))

site_pcoa <- plot_ordination(physeq = ps_hell, ordination = ps_hell_ord,
    shape = "AreaSample", color = "Site.name...8") + theme_bw() +
    scale_shape_manual(name = "Geographic Area", values = c(15:17)) +
    geom_point(size = 4) + scale_color_viridis_d(name = "Site",
    option = "B") + xlab("PCoA1 (15.1%)") + ylab("PCoA2 (5.9%)") +
    theme(text = element_text(size = 14))


site_pcoa
```

![](09_Antarctic_RMarkDown_files/figure-gfm/beta_plot-1.png)<!-- -->

``` r
ggsave("plots/exploratory/ant_pcoa_by_site.png", plot = last_plot(),
    device = "png", width = 6, height = 6, dpi = 300)
```

``` r
set.seed(5311)

# beta div stats
Dist.hell <- phyloseq::distance(ps_hell, method = "euclidean",
    type = "samples")

# here we will first test models with one variable to start
# adonis tests for differences in centroid and/or
# dispersion only variables we have without too many
# missing values are geography (site / lat / long), rock
# type, and year

# adonis2(Dist.hell ~ Rocks_v2, as(sample_data(ps_hell),
# 'data.frame'), permutations = 9999) #sig, ~ 4% variation
# adonis2(Dist.hell ~ Site.name...8,
# as(sample_data(ps_hell), 'data.frame'), permutations =
# 9999) #sig, ~ 33% adonis2(Dist.hell ~ Year.of.collection,
# as(sample_data(ps_hell), 'data.frame'), permutations =
# 9999) #sig 2% varation adonis2(Dist.hell ~ Latitude,
# as(sample_data(ps_hell), 'data.frame'), permutations =
# 9999) #sig 2% var adonis2(Dist.hell ~ Longitude,
# as(sample_data(ps_hell), 'data.frame'), permutations =
# 9999) #not sig adonis2(Dist.hell ~ Latitude*Longitude,
# as(sample_data(ps_hell), 'data.frame'), permutations =
# 9999) #sig 2% var each adonis2(Dist.hell ~
# Longitude*Latitude, as(sample_data(ps_hell),
# 'data.frame'), permutations = 9999) # latitude explains
# more 3%, interacts with long, long not sig here

# Sitename = a pseudo combination of lat & long so we
# wouldn't include all three in a model Sitename >> more
# variation explained than lat alone or lat/long combined
# so going to use Site in our next models rock type does
# not interact with year / site, so four possible models,
# given that site & year interact, best models have site
# before year in the model as explains more variation - so
# two models now adonis2(Dist.hell ~ Site.name...8 *
# Year.of.collection + Rocks_v2, as(sample_data(ps_hell),
# 'data.frame'), permutations = 9999) # sig site 33%, year
# <1% not sig, rocks 1%, sitexyear 5% adonis2(Dist.hell ~
# Rocks_v2 +Site.name...8 * Year.of.collection ,
# as(sample_data(ps_hell), 'data.frame'), permutations =
# 9999) # sig rock up 4%, but site down 30% so likely site
# explains that missing bit when first

# margin now to get contribution when other factors
# included vs. sequential check interactions first
# adonis2(formula = Dist.hell ~ Site.name...8 * Rocks_v2,
# data = as(sample_data(ps_hell), 'data.frame'),
# permutations = 9999, by = 'margin') #notsig
# adonis2(formula = Dist.hell ~ Year.of.collection *
# Rocks_v2, data = as(sample_data(ps_hell), 'data.frame'),
# permutations = 9999, by = 'margin') #sig - so probably
# can't use margin model if both in model adonis2(formula =
# Dist.hell ~ Year.of.collection * Site.name...8, data =
# as(sample_data(ps_hell), 'data.frame'), permutations =
# 9999, by = 'margin') #sig

# adonis2(formula = Dist.hell ~ Rocks_v2 + Site.name...8 +
# Year.of.collection, data = as(sample_data(ps_hell),
# 'data.frame'), permutations = 9999, by = 'margin') # rock
# 2.5%, site 29%, year not sig

# using ordistep & capscale to do a less hypothesis driven
# approach dbRDA (capscale), cca/rda != on dissimilarities
# upr <- capscale(Dist.hell ~ ., data =
# select(as(sample_data(ps_hell), 'data.frame'),
# 'Site.name...8', 'AreaSample', 'Rocks_v2',
# 'Year.of.collection')) lwr <- capscale(Dist.hell ~ 1,
# data = as(sample_data(ps_hell), 'data.frame')) model <-
# ordistep(lwr, scope = formula(upr), trace=0) #site +
# areasample + rocks is best model (but should site be
# nested in areasample?)  model2 <- ordiR2step(lwr, upr,
# trace=0) #Site.name...8 only as best model

# year was not in either model & rock was not in the second
# ordistep model! trying best models of just site + rock no
# interaction term with margin between these two

# with margin
adonis2(formula = Dist.hell ~ Site.name...8 + Rocks_v2, data = as(sample_data(ps_hell),
    "data.frame"), permutations = 9999, by = "margin")  #site sig and most explanitory 
```

    ## Permutation test for adonis under reduced model
    ## Marginal effects of terms
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = Dist.hell ~ Site.name...8 + Rocks_v2, data = as(sample_data(ps_hell), "data.frame"), permutations = 9999, by = "margin")
    ##                Df SumOfSqs      R2      F Pr(>F)    
    ## Site.name...8  34   46.761 0.30261 1.9552 0.0001 ***
    ## Rocks_v2        2    1.791 0.01159 1.2732 0.0373 *  
    ## Residual      144  101.295 0.65553                  
    ## Total         181  154.523 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# with terms
adonis2(formula = Dist.hell ~ Site.name...8 + Rocks_v2, data = as(sample_data(ps_hell),
    "data.frame"), permutations = 9999, by = "terms")  #site most important/explains most variation so put first in model
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## adonis2(formula = Dist.hell ~ Site.name...8 + Rocks_v2, data = as(sample_data(ps_hell), "data.frame"), permutations = 9999, by = "terms")
    ##                Df SumOfSqs      R2      F Pr(>F)    
    ## Site.name...8  35   51.437 0.33288 2.0892 0.0001 ***
    ## Rocks_v2        2    1.791 0.01159 1.2732 0.0399 *  
    ## Residual      144  101.295 0.65553                  
    ## Total         181  154.523 1.00000                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# no posthoc adonis test so need to do pair-wise
site.pair.res <- pairwise.adonis(x = Dist.hell, as.factor(as(sample_data(ps_hell),
    "data.frame")[, "Site.name...8"]), perm = 9999, p.adjust.m = "BH")
summary(site.pair.res$p.value < 0.05)
```

    ##    Mode   FALSE    TRUE    NA's 
    ## logical     399     221      10

``` r
write.csv(site.pair.res, "results/stats/site.pairwise.permanova.results.csv")

rock.pair.res <- pairwise.adonis(x = Dist.hell, as.factor(as(sample_data(ps_hell),
    "data.frame")[, "Rocks_v2"]), perm = 9999, p.adjust.m = "BH")
summary(rock.pair.res$p.value < 0.05)
```

    ##    Mode   FALSE    TRUE 
    ## logical       1       5

``` r
write.csv(rock.pair.res, "results/stats/rock.pairwise.permanova.results.csv")

# betadisper tests for dispersion only only can test single
# variables doing this as permanova are sensitive to
# dispersion differences
disp_dist_rock <- betadisper(Dist.hell, as(sample_data(ps_hell),
    "data.frame")[, "Rocks_v2"])
permutest(disp_dist_rock, permutations = 9999, pairwise = TRUE)
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Response: Distances
    ##            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)  
    ## Groups      3 0.10429 0.034763 2.8144   9999 0.0301 *
    ## Residuals 178 2.19862 0.012352                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##                 Basalt/Dolerite   Granite    Quartz Sandstone
    ## Basalt/Dolerite                 0.0043000 0.1288000    0.0002
    ## Granite               0.0161097           0.2339000    0.5259
    ## Quartz                0.1353131 0.2269756              0.2651
    ## Sandstone             0.0097774 0.5217567 0.2574563

``` r
# sig

disp_dist_site <- betadisper(Dist.hell, as(sample_data(ps_hell),
    "data.frame")[, "Site.name...8"])
permutest(disp_dist_site, permutations = 9999, pairwise = TRUE)  #sig
```

    ## 
    ## Permutation test for homogeneity of multivariate dispersions
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Response: Distances
    ##            Df Sum Sq  Mean Sq      F N.Perm Pr(>F)    
    ## Groups     35 4.7359 0.135312 4.3951   9999  1e-04 ***
    ## Residuals 146 4.4949 0.030787                         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Pairwise comparisons:
    ## (Observed p-value below diagonal, permuted p-value above diagonal)
    ##                       Machu Picchu Base Dufayel Island Helliwell Hills
    ## Machu Picchu Base                                           8.7990e-01
    ## Dufayel Island                                                        
    ## Helliwell Hills              8.7677e-01                               
    ## Crater Circle                6.6454e-29                     9.3081e-01
    ## Mt. Burrow                   3.2942e-28                     9.0097e-01
    ## Archambault Ridge            1.1484e-31                     1.1698e-01
    ## Random Hills                 9.8134e-03                     5.2431e-01
    ## Kay Island                   8.7295e-01                     7.7543e-01
    ## Mt. New Zealand              8.7070e-01                     9.2594e-01
    ## Mt. Dickason                 6.7331e-01                     7.1359e-01
    ## Mt. Nansen                   5.0929e-01                     7.6295e-01
    ## Mt. Keinath                                                           
    ## Anderson Ridge               3.4916e-32                     2.1486e-02
    ## Vegetation Island            1.3667e-01                     3.1938e-01
    ## Mt. Larsen                   6.1320e-32                     5.2540e-02
    ## Inexpressible Island         1.2255e-01                     1.3274e-01
    ## Mt. Billing                  8.5562e-01                     8.6658e-01
    ## Widowmaker Pass              7.5377e-01                     5.8270e-01
    ## Trio Nunatak                 8.0077e-01                     9.4890e-01
    ## Unnamed Nunatak              8.8191e-01                     9.9128e-01
    ## Harrow Peak                                                           
    ## Mt. McGee                    7.8802e-31                     3.2271e-01
    ## Ricker Hills                 5.2677e-01                     2.8159e-01
    ## Mt. Bowen                                                             
    ## Pudding Butte                4.7839e-01                     5.6600e-01
    ## Starr Nunatak                3.5912e-01                     5.3172e-01
    ## Richard Nunatak              6.3271e-31                     4.0939e-01
    ## Schultz Peak                 4.2210e-30                     6.7261e-01
    ## Battleship Promontory        1.4100e-01                     1.4029e-01
    ## Convoy Range                 8.9559e-01                     9.7093e-01
    ## Mt. Elektra                  1.2351e-02                     3.1328e-02
    ## Siegfried Peak               2.3853e-02                     4.0649e-01
    ## Linnaeus Terrace             1.4484e-01                     1.0849e-01
    ## Finger Mt.                   3.4195e-01                     4.1068e-01
    ## University Valley            3.7636e-01                     4.0300e-01
    ## Knobhead                                                              
    ##                       Crater Circle Mt. Burrow Archambault Ridge Random Hills
    ## Machu Picchu Base        1.0000e-04 1.0000e-04        1.0000e-04   9.9000e-03
    ## Dufayel Island                                                               
    ## Helliwell Hills          9.3400e-01 9.0330e-01        1.1320e-01   5.4260e-01
    ## Crater Circle                       1.0000e-04        1.0000e-04   1.1200e-02
    ## Mt. Burrow               0.0000e+00                   1.0000e-04   1.0300e-02
    ## Archambault Ridge        0.0000e+00 0.0000e+00                     5.0000e-04
    ## Random Hills             1.3252e-02 1.1192e-02        2.6011e-04             
    ## Kay Island               7.7506e-01 8.2846e-01        5.0867e-02   1.2649e-01
    ## Mt. New Zealand          9.0484e-01 8.8600e-01        2.6079e-01   7.4238e-01
    ## Mt. Dickason             7.1458e-01 6.9169e-01        1.2566e-01   8.6010e-01
    ## Mt. Nansen               5.7707e-01 5.3891e-01        3.3948e-02   5.9177e-01
    ## Mt. Keinath                                                                  
    ## Anderson Ridge           0.0000e+00 0.0000e+00        0.0000e+00   6.6501e-05
    ## Vegetation Island        1.6076e-01 1.4697e-01        6.4633e-03   5.6551e-01
    ## Mt. Larsen               0.0000e+00 0.0000e+00        0.0000e+00   1.2985e-04
    ## Inexpressible Island     1.4217e-01 1.3102e-01        3.9784e-03   3.7698e-01
    ## Mt. Billing              8.8111e-01 8.6705e-01        3.4785e-01   8.6655e-01
    ## Widowmaker Pass          7.1086e-01 7.3435e-01        3.1079e-01   2.8640e-01
    ## Trio Nunatak             8.6897e-01 8.3122e-01        4.0621e-02   4.5983e-01
    ## Unnamed Nunatak          9.3070e-01 9.0377e-01        1.4671e-01   5.7120e-01
    ## Harrow Peak                                                                  
    ## Mt. McGee                3.3969e-31 3.6811e-31        5.4246e-31   8.4465e-04
    ## Ricker Hills             5.0926e-01 5.1886e-01        9.6532e-01   2.6878e-01
    ## Mt. Bowen                                                                    
    ## Pudding Butte            5.4419e-01 5.0732e-01        6.1292e-03   7.0006e-01
    ## Starr Nunatak            4.1292e-01 3.8254e-01        1.2425e-02   8.7450e-01
    ## Richard Nunatak          0.0000e+00 0.0000e+00        0.0000e+00   1.2336e-03
    ## Schultz Peak             0.0000e+00 0.0000e+00        0.0000e+00   3.7107e-03
    ## Battleship Promontory    1.7311e-01 1.5475e-01        3.6140e-04   6.7621e-01
    ## Convoy Range             9.3175e-01 9.1180e-01        2.7357e-01   6.8156e-01
    ## Mt. Elektra              1.7576e-02 1.4477e-02        4.6642e-06   1.9359e-01
    ## Siegfried Peak           3.4763e-02 2.8190e-02        1.0189e-04   9.8123e-01
    ## Linnaeus Terrace         1.7337e-01 1.5714e-01        6.9195e-04   5.6205e-01
    ## Finger Mt.               3.9762e-01 3.6622e-01        3.4168e-03   8.9108e-01
    ## University Valley        4.2382e-01 3.9719e-01        1.1525e-02   9.7804e-01
    ## Knobhead                                                                     
    ##                       Kay Island Mt. New Zealand Mt. Dickason Mt. Nansen
    ## Machu Picchu Base     8.8490e-01      8.7490e-01   6.8880e-01 5.3170e-01
    ## Dufayel Island                                                          
    ## Helliwell Hills       7.9460e-01      9.3390e-01   7.2700e-01 7.7270e-01
    ## Crater Circle         7.9720e-01      9.1120e-01   7.3790e-01 6.0630e-01
    ## Mt. Burrow            8.4550e-01      8.9310e-01   7.0440e-01 5.6740e-01
    ## Archambault Ridge     4.4500e-02      2.6380e-01   1.1330e-01 2.6400e-02
    ## Random Hills          1.1840e-01      7.6310e-01   8.7310e-01 6.2130e-01
    ## Kay Island                            8.0810e-01   5.8640e-01 4.5690e-01
    ## Mt. New Zealand       7.8941e-01                   8.5320e-01 9.0800e-01
    ## Mt. Dickason          5.5469e-01      8.4276e-01              9.3200e-01
    ## Mt. Nansen            4.3353e-01      9.0303e-01   9.2687e-01           
    ## Mt. Keinath                                                             
    ## Anderson Ridge        9.8486e-03      8.3246e-02   3.5091e-02 8.7373e-03
    ## Vegetation Island     1.0365e-01      5.5527e-01   6.3477e-01 3.8946e-01
    ## Mt. Larsen            2.1912e-02      1.5548e-01   6.8549e-02 1.7136e-02
    ## Inexpressible Island  6.3300e-02      3.0384e-01   3.6878e-01 2.3781e-01
    ## Mt. Billing           7.8424e-01      9.3055e-01   9.3646e-01 9.8574e-01
    ## Widowmaker Pass       7.7134e-01      6.1870e-01   4.3375e-01 4.4204e-01
    ## Trio Nunatak          6.7146e-01      9.4372e-01   6.7818e-01 7.4761e-01
    ## Unnamed Nunatak       7.8660e-01      9.3497e-01   7.3485e-01 7.9028e-01
    ## Harrow Peak                                                             
    ## Mt. McGee             2.0242e-01      4.7982e-01   2.7895e-01 1.0041e-01
    ## Ricker Hills          4.5880e-01      2.9239e-01   2.5698e-01 3.2830e-01
    ## Mt. Bowen                                                               
    ## Pudding Butte         3.2852e-01      7.6322e-01   9.4696e-01 9.2720e-01
    ## Starr Nunatak         2.5950e-01      7.3718e-01   9.1361e-01 7.6310e-01
    ## Richard Nunatak       3.0105e-01      5.5134e-01   3.3959e-01 1.3784e-01
    ## Schultz Peak          7.3195e-01      7.3910e-01   5.2378e-01 3.0429e-01
    ## Battleship Promontory 6.2251e-02      3.1503e-01   4.9962e-01 3.6419e-01
    ## Convoy Range          8.1552e-01      9.6754e-01   8.1480e-01 8.5912e-01
    ## Mt. Elektra           3.8770e-03      1.5345e-01   2.1174e-01 7.4667e-02
    ## Siegfried Peak        5.9891e-02      6.6425e-01   8.1819e-01 5.0581e-01
    ## Linnaeus Terrace      6.3323e-02      2.4154e-01   4.0230e-01 3.1266e-01
    ## Finger Mt.            2.1241e-01      6.2623e-01   8.7557e-01 7.3069e-01
    ## University Valley     2.4734e-01      5.9707e-01   7.9932e-01 6.7812e-01
    ## Knobhead                                                                
    ##                       Mt. Keinath Anderson Ridge Vegetation Island Mt. Larsen
    ## Machu Picchu Base                     1.0000e-04        1.2700e-01 1.0000e-04
    ## Dufayel Island                                                               
    ## Helliwell Hills                       2.2800e-02        3.3120e-01 4.9900e-02
    ## Crater Circle                         1.0000e-04        1.5570e-01 1.0000e-04
    ## Mt. Burrow                            1.0000e-04        1.4130e-01 1.0000e-04
    ## Archambault Ridge                     1.0000e-04        5.1000e-03 1.0000e-04
    ## Random Hills                          1.0000e-04        5.9130e-01 2.0000e-04
    ## Kay Island                            8.5000e-03        9.5800e-02 1.7900e-02
    ## Mt. New Zealand                       8.2600e-02        5.8390e-01 1.5280e-01
    ## Mt. Dickason                          3.1600e-02        6.6290e-01 5.9000e-02
    ## Mt. Nansen                            7.2000e-03        4.0620e-01 1.4900e-02
    ## Mt. Keinath                                                                  
    ## Anderson Ridge                                          1.1000e-03 1.0000e-04
    ## Vegetation Island                     1.3621e-03                   2.5000e-03
    ## Mt. Larsen                            0.0000e+00        2.9778e-03           
    ## Inexpressible Island                  4.4117e-04        6.2239e-01 1.3741e-03
    ## Mt. Billing                           1.3462e-01        7.2060e-01 2.3041e-01
    ## Widowmaker Pass                       7.6319e-02        1.6172e-01 1.6235e-01
    ## Trio Nunatak                          2.9506e-03        2.3649e-01 1.1845e-02
    ## Unnamed Nunatak                       3.0588e-02        3.6445e-01 7.0513e-02
    ## Harrow Peak                                                                  
    ## Mt. McGee                             7.0788e-32        2.1413e-02 1.6725e-31
    ## Ricker Hills                          5.5270e-01        1.6195e-01 7.6735e-01
    ## Mt. Bowen                                                                    
    ## Pudding Butte                         1.5429e-04        3.6659e-01 1.0851e-03
    ## Starr Nunatak                         1.7151e-03        5.3479e-01 4.6731e-03
    ## Richard Nunatak                       0.0000e+00        3.0299e-02 0.0000e+00
    ## Schultz Peak                          0.0000e+00        7.3438e-02 0.0000e+00
    ## Battleship Promontory                 3.5326e-06        9.0791e-01 4.0533e-05
    ## Convoy Range                          1.0416e-01        4.9537e-01 1.7362e-01
    ## Mt. Elektra                           4.0258e-08        5.6698e-01 4.5827e-07
    ## Siegfried Peak                        1.0903e-05        4.6590e-01 3.2781e-05
    ## Linnaeus Terrace                      7.9371e-06        8.8907e-01 8.6734e-05
    ## Finger Mt.                            9.9050e-05        5.1125e-01 6.3272e-04
    ## University Valley                     7.6325e-04        6.9111e-01 3.1891e-03
    ## Knobhead                                                                     
    ##                       Inexpressible Island Mt. Billing Widowmaker Pass
    ## Machu Picchu Base               1.2280e-01  8.6450e-01      7.6830e-01
    ## Dufayel Island                                                        
    ## Helliwell Hills                 1.2720e-01  8.6620e-01      5.9950e-01
    ## Crater Circle                   1.3800e-01  8.8200e-01      7.2990e-01
    ## Mt. Burrow                      1.1900e-01  8.6500e-01      7.5700e-01
    ## Archambault Ridge               5.4000e-03  3.3490e-01      3.1790e-01
    ## Random Hills                    3.9160e-01  8.6910e-01      2.9360e-01
    ## Kay Island                      6.1500e-02  7.8870e-01      7.9370e-01
    ## Mt. New Zealand                 3.1610e-01  9.3510e-01      6.4170e-01
    ## Mt. Dickason                    3.8250e-01  9.3690e-01      4.5400e-01
    ## Mt. Nansen                      2.4150e-01  9.8750e-01      4.5600e-01
    ## Mt. Keinath                                                           
    ## Anderson Ridge                  1.2000e-03  1.3650e-01      6.8600e-02
    ## Vegetation Island               6.4550e-01  7.3770e-01      1.5610e-01
    ## Mt. Larsen                      2.8000e-03  2.1850e-01      1.5970e-01
    ## Inexpressible Island                        4.9760e-01      4.8400e-02
    ## Mt. Billing                     4.8292e-01                  6.3400e-01
    ## Widowmaker Pass                 5.5612e-02  6.2475e-01                
    ## Trio Nunatak                    7.2563e-02  8.6045e-01      4.4843e-01
    ## Unnamed Nunatak                 1.5821e-01  8.7347e-01      5.9583e-01
    ## Harrow Peak                                                           
    ## Mt. McGee                       1.7776e-02  5.5009e-01      6.6195e-01
    ## Ricker Hills                    4.6635e-02  2.2975e-01      5.0398e-01
    ## Mt. Bowen                                                             
    ## Pudding Butte                   1.0818e-01  9.2404e-01      1.8102e-01
    ## Starr Nunatak                   2.7098e-01  8.7457e-01      2.4460e-01
    ## Richard Nunatak                 2.6447e-02  6.0932e-01      7.7834e-01
    ## Schultz Peak                    6.7523e-02  7.5647e-01      9.3312e-01
    ## Battleship Promontory           4.0113e-01  5.1833e-01      2.7855e-02
    ## Convoy Range                    2.7450e-01  9.1684e-01      6.6836e-01
    ## Mt. Elektra                     7.9263e-01  3.5662e-01      4.8855e-03
    ## Siegfried Peak                  2.5360e-01  8.2817e-01      1.6059e-01
    ## Linnaeus Terrace                6.2269e-01  4.0861e-01      2.3587e-02
    ## Finger Mt.                      1.8270e-01  8.0748e-01      1.2388e-01
    ## University Valley               3.2824e-01  7.6247e-01      1.5369e-01
    ## Knobhead                                                              
    ##                       Trio Nunatak Unnamed Nunatak Harrow Peak  Mt. McGee
    ## Machu Picchu Base       8.0070e-01      8.8970e-01             1.0000e-04
    ## Dufayel Island                                                           
    ## Helliwell Hills         9.5250e-01      9.9190e-01             3.2640e-01
    ## Crater Circle           8.7220e-01      9.3280e-01             1.0000e-04
    ## Mt. Burrow              8.3510e-01      9.0280e-01             1.0000e-04
    ## Archambault Ridge       3.7300e-02      1.3770e-01             1.0000e-04
    ## Random Hills            4.6950e-01      5.9440e-01             9.0000e-04
    ## Kay Island              6.8220e-01      8.0380e-01             1.9140e-01
    ## Mt. New Zealand         9.4540e-01      9.4040e-01             4.9030e-01
    ## Mt. Dickason            6.8250e-01      7.5400e-01             2.7570e-01
    ## Mt. Nansen              7.5770e-01      8.0360e-01             8.6300e-02
    ## Mt. Keinath                                                              
    ## Anderson Ridge          5.4000e-03      3.1300e-02             1.0000e-04
    ## Vegetation Island       2.4410e-01      3.7470e-01             1.6800e-02
    ## Mt. Larsen              1.4900e-02      6.7700e-02             1.0000e-04
    ## Inexpressible Island    6.4800e-02      1.5800e-01             1.8400e-02
    ## Mt. Billing             8.6600e-01      8.7830e-01             5.4570e-01
    ## Widowmaker Pass         4.6660e-01      6.1520e-01             6.8340e-01
    ## Trio Nunatak                            9.6470e-01             1.8340e-01
    ## Unnamed Nunatak         9.6254e-01                             3.6820e-01
    ## Harrow Peak                                                              
    ## Mt. McGee               1.8869e-01      3.6381e-01                       
    ## Ricker Hills            1.5748e-01      2.8515e-01             7.8552e-01
    ## Mt. Bowen                                                                
    ## Pudding Butte           5.3593e-01      5.9555e-01             5.4732e-02
    ## Starr Nunatak           4.7770e-01      5.7014e-01             5.2879e-02
    ## Richard Nunatak         2.6857e-01      4.4899e-01             8.8614e-30
    ## Schultz Peak            5.5239e-01      6.9627e-01             8.1695e-31
    ## Battleship Promontory   8.6387e-02      1.6343e-01             6.2788e-03
    ## Convoy Range            9.9953e-01      9.7881e-01             4.8849e-01
    ## Mt. Elektra             1.2451e-02      4.4613e-02             1.3452e-04
    ## Siegfried Peak          3.4485e-01      4.5844e-01             6.7200e-04
    ## Linnaeus Terrace        5.9007e-02      1.2490e-01             9.3858e-03
    ## Finger Mt.              3.5382e-01      4.4431e-01             3.1557e-02
    ## University Valley       3.3480e-01      4.3321e-01             6.1526e-02
    ## Knobhead                                                                 
    ##                       Ricker Hills Mt. Bowen Pudding Butte Starr Nunatak
    ## Machu Picchu Base       5.2770e-01              4.6770e-01    3.6550e-01
    ## Dufayel Island                                                          
    ## Helliwell Hills         2.8160e-01              5.7330e-01    5.6220e-01
    ## Crater Circle           5.1640e-01              5.3580e-01    4.2700e-01
    ## Mt. Burrow              5.2200e-01              4.8910e-01    3.8760e-01
    ## Archambault Ridge       9.6650e-01              1.1900e-02    1.3500e-02
    ## Random Hills            2.6910e-01              7.0390e-01    8.8550e-01
    ## Kay Island              4.7520e-01              3.2580e-01    2.6830e-01
    ## Mt. New Zealand         2.9310e-01              7.6650e-01    7.5230e-01
    ## Mt. Dickason            2.6080e-01              9.4520e-01    9.1800e-01
    ## Mt. Nansen              3.2710e-01              9.3210e-01    7.8220e-01
    ## Mt. Keinath                                                             
    ## Anderson Ridge          5.6000e-01              1.0000e-03    2.0000e-03
    ## Vegetation Island       1.6210e-01              3.8130e-01    5.5770e-01
    ## Mt. Larsen              7.6890e-01              3.9000e-03    5.3000e-03
    ## Inexpressible Island    4.0400e-02              1.0870e-01    2.8320e-01
    ## Mt. Billing             2.2680e-01              9.3020e-01    8.8020e-01
    ## Widowmaker Pass         5.2190e-01              1.7970e-01    2.5790e-01
    ## Trio Nunatak            1.5560e-01              5.5220e-01    4.9410e-01
    ## Unnamed Nunatak         2.8960e-01              6.0080e-01    5.8730e-01
    ## Harrow Peak                                                             
    ## Mt. McGee               7.8680e-01              5.6800e-02    4.8000e-02
    ## Ricker Hills                                    4.7800e-02    1.6970e-01
    ## Mt. Bowen                                                               
    ## Pudding Butte           4.7111e-02                            7.9330e-01
    ## Starr Nunatak           1.6864e-01              7.8338e-01              
    ## Richard Nunatak         7.2676e-01              9.1603e-02    7.8591e-02
    ## Schultz Peak            5.9996e-01              2.6925e-01    2.0092e-01
    ## Battleship Promontory   7.1504e-03              2.1166e-01    5.0092e-01
    ## Convoy Range            4.0649e-01              7.0588e-01    6.7881e-01
    ## Mt. Elektra             4.8512e-03              3.2412e-02    1.2342e-01
    ## Siegfried Peak          1.4791e-01              6.2813e-01    8.4843e-01
    ## Linnaeus Terrace        3.9698e-03              1.2889e-01    3.8784e-01
    ## Finger Mt.              4.2221e-02              7.0669e-01    9.8415e-01
    ## University Valley       6.1373e-02              5.9529e-01    8.5891e-01
    ## Knobhead                                                                
    ##                       Richard Nunatak Schultz Peak Battleship Promontory
    ## Machu Picchu Base          1.0000e-04   1.0000e-04            1.3780e-01
    ## Dufayel Island                                                          
    ## Helliwell Hills            4.1870e-01   6.7550e-01            1.3640e-01
    ## Crater Circle              1.0000e-04   1.0000e-04            1.6800e-01
    ## Mt. Burrow                 1.0000e-04   1.0000e-04            1.4590e-01
    ## Archambault Ridge          1.0000e-04   1.0000e-04            2.9000e-03
    ## Random Hills               1.6000e-03   3.3000e-03            6.7600e-01
    ## Kay Island                 3.0960e-01   7.5110e-01            6.7400e-02
    ## Mt. New Zealand            5.6310e-01   7.5020e-01            3.2020e-01
    ## Mt. Dickason               3.4270e-01   5.3610e-01            5.0080e-01
    ## Mt. Nansen                 1.2450e-01   3.0480e-01            3.5740e-01
    ## Mt. Keinath                                                             
    ## Anderson Ridge             1.0000e-04   1.0000e-04            1.0000e-04
    ## Vegetation Island          2.8300e-02   6.4300e-02            9.0740e-01
    ## Mt. Larsen                 1.0000e-04   1.0000e-04            1.0000e-04
    ## Inexpressible Island       2.7200e-02   6.7300e-02            4.0860e-01
    ## Mt. Billing                6.0080e-01   7.5020e-01            5.2510e-01
    ## Widowmaker Pass            8.0020e-01   9.4010e-01            2.8800e-02
    ## Trio Nunatak               2.6270e-01   5.5630e-01            8.4700e-02
    ## Unnamed Nunatak            4.5370e-01   7.0830e-01            1.5970e-01
    ## Harrow Peak                                                             
    ## Mt. McGee                  1.0000e-04   1.0000e-04            9.1000e-03
    ## Ricker Hills               7.2860e-01   6.0390e-01            7.2000e-03
    ## Mt. Bowen                                                               
    ## Pudding Butte              9.8500e-02   2.5800e-01            2.1710e-01
    ## Starr Nunatak              7.6000e-02   2.0340e-01            5.0610e-01
    ## Richard Nunatak                         1.0000e-04            2.0700e-02
    ## Schultz Peak               0.0000e+00                         6.2900e-02
    ## Battleship Promontory      1.2649e-02   5.8656e-02                      
    ## Convoy Range               5.6118e-01   7.5649e-01            2.8535e-01
    ## Mt. Elektra                3.3550e-04   3.0040e-03            3.4964e-01
    ## Siegfried Peak             1.2141e-03   6.3222e-03            5.8472e-01
    ## Linnaeus Terrace           1.7498e-02   6.7321e-02            6.8791e-01
    ## Finger Mt.                 5.4418e-02   1.7716e-01            4.0564e-01
    ## University Valley          9.2662e-02   2.2639e-01            6.2573e-01
    ## Knobhead                                                                
    ##                       Convoy Range Mt. Elektra Siegfried Peak Linnaeus Terrace
    ## Machu Picchu Base       9.0490e-01  1.7900e-02     2.1000e-02       1.3830e-01
    ## Dufayel Island                                                                
    ## Helliwell Hills         9.7140e-01  2.8500e-02     4.1840e-01       1.0330e-01
    ## Crater Circle           9.3670e-01  2.3200e-02     3.3700e-02       1.6820e-01
    ## Mt. Burrow              9.1910e-01  2.1000e-02     2.8600e-02       1.4610e-01
    ## Archambault Ridge       2.7210e-01  1.0000e-04     5.0000e-04       3.3000e-03
    ## Random Hills            7.0220e-01  1.9190e-01     9.8320e-01       5.6210e-01
    ## Kay Island              8.3620e-01  7.0000e-03     5.4100e-02       6.5700e-02
    ## Mt. New Zealand         9.6780e-01  1.5250e-01     6.8120e-01       2.4130e-01
    ## Mt. Dickason            8.2290e-01  2.0710e-01     8.2310e-01       3.9950e-01
    ## Mt. Nansen              8.7260e-01  7.2400e-02     5.2660e-01       3.0030e-01
    ## Mt. Keinath                                                                   
    ## Anderson Ridge          9.0400e-02  1.0000e-04     1.0000e-04       4.0000e-04
    ## Vegetation Island       5.2150e-01  5.7490e-01     4.8810e-01       8.9660e-01
    ## Mt. Larsen              1.6790e-01  2.0000e-04     1.0000e-04       1.1000e-03
    ## Inexpressible Island    2.7430e-01  7.9670e-01     2.6220e-01       6.2490e-01
    ## Mt. Billing             9.2070e-01  3.6080e-01     8.3280e-01       4.1690e-01
    ## Widowmaker Pass         6.9120e-01  5.5000e-03     1.6230e-01       2.4700e-02
    ## Trio Nunatak            9.9950e-01  1.0900e-02     3.4920e-01       5.7900e-02
    ## Unnamed Nunatak         9.8220e-01  3.8200e-02     4.7680e-01       1.2110e-01
    ## Harrow Peak                                                                   
    ## Mt. McGee               5.0500e-01  1.3000e-03     4.0000e-04       1.5200e-02
    ## Ricker Hills            4.1650e-01  4.6000e-03     1.4420e-01       4.6000e-03
    ## Mt. Bowen                                                                     
    ## Pudding Butte           7.1090e-01  3.1000e-02     6.2840e-01       1.3550e-01
    ## Starr Nunatak           6.9680e-01  1.2260e-01     8.5700e-01       3.8580e-01
    ## Richard Nunatak         5.7980e-01  2.3000e-03     2.1000e-03       2.8900e-02
    ## Schultz Peak            7.7550e-01  7.4000e-03     7.0000e-03       7.1100e-02
    ## Battleship Promontory   2.7380e-01  3.5490e-01     5.9610e-01       6.9180e-01
    ## Convoy Range                        1.0730e-01     6.0490e-01       2.2680e-01
    ## Mt. Elektra             1.1076e-01                 9.7600e-02       6.8360e-01
    ## Siegfried Peak          5.7788e-01  1.0244e-01                      4.5870e-01
    ## Linnaeus Terrace        2.3868e-01  6.6894e-01     4.5064e-01                 
    ## Finger Mt.              5.7143e-01  8.1135e-02     8.6918e-01       2.6317e-01
    ## University Valley       5.5993e-01  2.1257e-01     9.6428e-01       4.4608e-01
    ## Knobhead                                                                      
    ##                       Finger Mt. University Valley Knobhead
    ## Machu Picchu Base     3.2370e-01        3.6830e-01         
    ## Dufayel Island                                             
    ## Helliwell Hills       4.1400e-01        4.0700e-01         
    ## Crater Circle         3.8640e-01        4.2600e-01         
    ## Mt. Burrow            3.5400e-01        3.9070e-01         
    ## Archambault Ridge     7.0000e-03        1.4900e-02         
    ## Random Hills          8.9250e-01        9.7840e-01         
    ## Kay Island            2.1530e-01        2.5150e-01         
    ## Mt. New Zealand       6.3930e-01        6.1920e-01         
    ## Mt. Dickason          8.7640e-01        8.0080e-01         
    ## Mt. Nansen            7.2740e-01        6.8700e-01         
    ## Mt. Keinath                                                
    ## Anderson Ridge        1.0000e-03        2.7000e-03         
    ## Vegetation Island     5.2250e-01        7.0100e-01         
    ## Mt. Larsen            2.1000e-03        6.0000e-03         
    ## Inexpressible Island  1.8100e-01        3.4280e-01         
    ## Mt. Billing           8.1290e-01        7.6800e-01         
    ## Widowmaker Pass       1.2460e-01        1.6120e-01         
    ## Trio Nunatak          3.6520e-01        3.4530e-01         
    ## Unnamed Nunatak       4.3980e-01        4.3910e-01         
    ## Harrow Peak                                                
    ## Mt. McGee             3.5800e-02        6.3700e-02         
    ## Ricker Hills          3.7000e-02        5.6800e-02         
    ## Mt. Bowen                                                  
    ## Pudding Butte         7.1240e-01        6.1100e-01         
    ## Starr Nunatak         9.8550e-01        8.6010e-01         
    ## Richard Nunatak       5.8200e-02        9.2800e-02         
    ## Schultz Peak          1.7540e-01        2.1830e-01         
    ## Battleship Promontory 4.1200e-01        6.3520e-01         
    ## Convoy Range          5.6450e-01        5.7320e-01         
    ## Mt. Elektra           7.9100e-02        2.2110e-01         
    ## Siegfried Peak        8.7500e-01        9.6420e-01         
    ## Linnaeus Terrace      2.6840e-01        4.6070e-01         
    ## Finger Mt.                              8.3850e-01         
    ## University Valley     8.3585e-01                           
    ## Knobhead

``` r
# post hoc test to see which pairs driving
tuk_disp_dist_rock <- TukeyHSD(disp_dist_rock)
tuk_disp_dist_site <- TukeyHSD(disp_dist_site)

summary(tidy(tuk_disp_dist_rock)$adj.p.value < 0.05)
```

    ##    Mode   FALSE    TRUE 
    ## logical       5       1

``` r
summary(tidy(tuk_disp_dist_site)$adj.p.value < 0.05)
```

    ##    Mode   FALSE    TRUE 
    ## logical     567      63

``` r
write.csv(tidy(tuk_disp_dist_rock), "results/stats/rock.tukeyposthoc.dispersion.csv")
write.csv(tidy(tuk_disp_dist_site), "results/stats/site.tukeyposthoc.dispersion.csv")
```

# Distance Decay

``` r
# distancy decay
ps.dist <- ps.HQ.long

# since most samples on one coast, removing rock type as
# possible variable
ps.HQ.EC <- subset_samples(ps.dist, AreaSample != "Antarctic Peninsula")
ps.HQ.EC <- subset_samples(ps.HQ.EC, Rocks_v2 == "Sandstone")

# get lat and long
geo = data.frame(ps.HQ.EC@sam_data$Longitude, ps.HQ.EC@sam_data$Latitude)

# defining geographic distance between two locations -
# haversine distance (accounts for spherical earth)
d.geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)

ps_hell <- transform(ps.HQ.EC, "hellinger")

# getting community distances have to multiply by 1/sqrt(2)
# to get scaled 0-1 as current range is 0-sqrt(2)
comm.dist.hell = vegdist(1/sqrt(2) * as.matrix(t(ps_hell@otu_table)),
    method = "euclidean")

# test
mantel.test.hell = vegan::mantel(comm.dist.hell, dist.geo, method = "spearman",
    permutations = 9999, na.rm = TRUE)
mantel.test.hell
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## vegan::mantel(xdis = comm.dist.hell, ydis = dist.geo, method = "spearman",      permutations = 9999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1966 
    ##       Significance: 1e-04 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0458 0.0595 0.0719 0.0890 
    ## Permutation: free
    ## Number of permutations: 9999

``` r
mantel.stats <- data.frame(label = paste("r = ", signif(mantel.test.hell$statistic,
    3), "\nP = ", signif(mantel.test.hell$signif, 3)))
# plot
dist_hell <- ggplot(mapping = aes(x = jitter(dist.geo/1000, amount = 1),
    y = comm.dist.hell)) + theme_bw() + geom_point(shape = 16,
    size = 1, alpha = 0.1, color = "gray25") + geom_smooth(method = "lm",
    color = "#0292e3", se = F) + labs(x = "Geographical Separation (km)",
    y = "Hellinger Distance") + geom_text(data = mantel.stats,
    aes(x = 475, y = 0.15, label = label), hjust = 0, size = 3) +
    theme(text = element_text(size = 14))

dist_hell
```

    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.

    ## `geom_smooth()` using formula 'y ~ x'

![](09_Antarctic_RMarkDown_files/figure-gfm/dist-1.png)<!-- -->

``` r
# what happens if we remove all instances from same
# location (e.g. geographic distance = 0)?

# first transform to matrices to edit
dist.geo.matrix <- as.matrix(dist.geo)
comm.dist.matrix <- as.matrix(comm.dist.hell)

# change all 0s to NAs
dist.geo.matrix[dist.geo.matrix == 0] <- NA
comm.dist.matrix[is.na(dist.geo.matrix)] <- NA

# change back to distances
dist.geo.nonzero <- as.dist(dist.geo.matrix)
comm.dist.nonzero <- as.dist(comm.dist.matrix)

# mantel test with no 0s bc using na.rm = TRUE
mantel.test.nonzero = vegan::mantel(comm.dist.nonzero, dist.geo.nonzero,
    method = "spearman", permutations = 9999, na.rm = TRUE)

mantel.test.nonzero
```

    ## 
    ## Mantel statistic based on Spearman's rank correlation rho 
    ## 
    ## Call:
    ## vegan::mantel(xdis = comm.dist.nonzero, ydis = dist.geo.nonzero,      method = "spearman", permutations = 9999, na.rm = TRUE) 
    ## 
    ## Mantel statistic r: 0.1701 
    ##       Significance: 1e-04 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0495 0.0659 0.0793 0.0953 
    ## Permutation: free
    ## Number of permutations: 9999

# Making a map of samples

``` r
worldmap <- ne_countries(scale = "medium", returnclass = "sf")

antartcia_sf <- worldmap[worldmap$name == "Antarctica", ]
# antartcia <- map_data('world') %>%
# filter(region=='Antarctica')

meta <- meta %>%
    arrange(desc(Latitude))

meta$Site <- factor(meta$Site.name...8, levels = c("Machu Picchu Base",
    "Dufayel Island", "Helliwell Hills", "Crater Circle", "Mt. Burrow",
    "Archambault Ridge", "Random Hills", "Kay Island", "Timber Peak",
    "Mt. New Zealand", "Mt. Dickason", "Mt. Nansen", "Mt. Keinath",
    "Anderson Ridge", "Vegetation Island", "Mt. Larsen", "Inexpressible Island",
    "Mt. Billing", "Widowmaker Pass", "Trio Nunatak", "Unnamed Nunatak",
    "Harrow Peak", "Mt. McGee", "Ricker Hills", "Mt. Bowen",
    "Pudding Butte", "Starr Nunatak", "Richard Nunatak", "Schultz Peak",
    "Battleship Promontory", "Convoy Range", "Mt. Elektra", "Siegfried Peak",
    "Linnaeus Terrace", "Finger Mt.", "University Valley", "Knobhead"))


meta <- meta %>%
    mutate(AreaSample = ifelse(AreaSample == "Northen Victoria Land",
        "Northern Victoria Land", as.character(AreaSample)))

as.big.nocol <- ggplot() + geom_sf(data = antartcia_sf) + coord_sf(expand = FALSE) +
    geom_point(data = meta, aes(x = Longitude, y = Latitude,
        shape = AreaSample, color = AreaSample)) + scale_shape_manual(name = "Geographic Area",
    values = c(15:17)) + scale_color_manual(name = "Geographic Area",
    values = c("black", "#FCFFA4", "#BB3754")) + theme_bw() +
    xlab("Longitude") + ylab("Latitude") + scale_y_continuous(limits = c(-85,
    -60)) + theme(panel.grid.major = element_line(color = gray(0.5),
    linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "#f5f9fa")) +
    theme(text = element_text(size = 14))


as.big.nocol
```

![](09_Antarctic_RMarkDown_files/figure-gfm/map-1.png)<!-- -->

``` r
ggsave(filename = "plots/exploratory/art_map.png", plot = last_plot(),
    device = "png", width = 10, height = 6, dpi = 300)
```

# Figure 2

``` r
addSmallLegend <- function(myPlot, pointSize = 3, textSize = 10,
    spaceLegend = 0.5) {
    myPlot + guides(shape = guide_legend(override.aes = list(size = pointSize)),
        color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = textSize), legend.text = element_text(size = textSize),
            legend.key.size = unit(spaceLegend, "lines"))
}


# flip RA plot
site_flip <- site_bar_ra + coord_flip()

(addSmallLegend(as.big.nocol + theme(legend.position = "top")))/(addSmallLegend(site_flip +
    theme(legend.position = "left")) + ((addSmallLegend(site_pcoa)/dist_hell))) +
    plot_annotation(tag_levels = "A")
```

    ## Warning: Duplicated override.aes is ignored.

    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.
    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.

    ## `geom_smooth()` using formula 'y ~ x'

![](09_Antarctic_RMarkDown_files/figure-gfm/fig2_plots-1.png)<!-- -->

``` r
ggsave(filename = "plots/fig2.pdf", plot = last_plot(), device = "pdf",
    width = 20, height = 10, dpi = 300)
```

    ## Warning: Duplicated override.aes is ignored.

    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.

    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.

    ## `geom_smooth()` using formula 'y ~ x'

``` r
ggsave(filename = "plots/fig2.png", plot = last_plot(), device = "png",
    width = 20, height = 10, dpi = 300)
```

    ## Warning: Duplicated override.aes is ignored.

    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.

    ## Don't know how to automatically pick scale for object of type dist. Defaulting to continuous.

    ## `geom_smooth()` using formula 'y ~ x'

``` r
# without map (site_flip + ( site_pcoa / dist_hell )) +
# plot_annotation(tag_levels = 'A') + plot_layout(guides =
# 'collect')
```

# Data tables

``` r
## Tables to summarize all the info!
write.csv(combo, "results/vOTU_info.csv")

# whole dataset

combo2 <- combo %>%
    filter(!is.na(votu.id))

# total number of unique VCs
length(unique(combo2$VCStatus)) - 1  #Unclustered
```

    ## [1] 10624

``` r
# total number of vOTUs
length(unique(combo2$votu.id))
```

    ## [1] 76984

``` r
# total number of viral seqs
length(unique(combo2$Genome))
```

    ## [1] 101085

``` r
checkv_cdhit_tax.meta.samplesum <- combo %>%
    group_by(SampleID) %>%
    summarise(contig_length_avg = mean(contig_length), gene_count_avg = mean(gene_count),
        num_VCs = length(na.omit(unique(VCStatus))), num_vOTUs = length(na.omit(unique(votu.id))),
        num_seqs = length(na.omit(unique(Genome))))

number_provirus <- combo %>%
    group_by(SampleID, provirus) %>%
    tally() %>%
    filter(provirus == "Yes")

number_provirus$provirus_num <- number_provirus$n

number_check <- combo %>%
    group_by(SampleID, checkv_quality) %>%
    tally() %>%
    spread(checkv_quality, n) %>%
    replace(is.na(.), 0)

checkv_cdhit_tax.meta.samplesum.combo <- checkv_cdhit_tax.meta.samplesum %>%
    left_join(dplyr::select(number_provirus, SampleID, provirus_num)) %>%
    left_join(number_check[-c(7)]) %>%
    left_join(dplyr::select(sampledata2, SampleID, AreaSample,
        Site.name...8, `Year of collection`, Rocks_v2, Latitude,
        Longitude)) %>%
    replace(is.na(.), 0)
```

    ## Joining, by = "SampleID"
    ## Joining, by = "SampleID"
    ## Joining, by = "SampleID"

``` r
# contig length
summary(checkv_cdhit_tax.meta.samplesum.combo$contig_length_avg)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##       0    9415   10038    9900   10878   17350

``` r
# avg gene count per seq
summary(checkv_cdhit_tax.meta.samplesum.combo$gene_count_avg)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    0.00   11.24   12.11   11.90   12.92   20.28

``` r
# avg number viral seq
summary(checkv_cdhit_tax.meta.samplesum.combo$num_seqs)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##     0.0   283.5   452.0   529.2   766.0  1514.0

``` r
# total number of viral seq
sum(checkv_cdhit_tax.meta.samplesum.combo$num_seqs)
```

    ## [1] 101085

``` r
# avg number provirus
summary(checkv_cdhit_tax.meta.samplesum.combo$provirus_num)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    0.00   10.50   20.00   24.99   34.00   96.00

``` r
# total number of provirus
sum(checkv_cdhit_tax.meta.samplesum.combo$provirus_num)
```

    ## [1] 4773

``` r
# total num complete genomes
sum(checkv_cdhit_tax.meta.samplesum.combo$Complete)
```

    ## [1] 77

``` r
# total number high-quality
sum(checkv_cdhit_tax.meta.samplesum.combo$`High-quality`)
```

    ## [1] 366

``` r
# total number medium-quality
sum(checkv_cdhit_tax.meta.samplesum.combo$`Medium-quality`)
```

    ## [1] 1600

``` r
write.csv(checkv_cdhit_tax.meta.samplesum.combo, "results/sample_sum_info.csv")




# filtered dataset
combo.filt <- combo %>%
    filter(!checkv_quality %in% c("Not-determined")) %>%
    filter(is.na(contig_length) | contig_length >= 10000) %>%
    filter(!is.na(votu.id))

# total number of unique VCs
length(unique(combo.filt$VCStatus)) - 1  #Unclustered
```

    ## [1] 2851

``` r
# total number of vOTUs
length(unique(combo.filt$votu.id))
```

    ## [1] 11806

``` r
# total number of viral seqs
length(unique(combo.filt$Genome))
```

    ## [1] 14796

``` r
checkv_cdhit_tax.meta.samplesum.filt <- combo.filt %>%
    group_by(SampleID) %>%
    summarise(contig_length_avg = mean(contig_length), gene_count_avg = mean(gene_count),
        num_VCs = length(na.omit(unique(VCStatus))), num_vOTUs = length(na.omit(unique(votu.id))),
        num_seqs = length(na.omit(unique(Genome))))

number_provirus.filt <- combo.filt %>%
    group_by(SampleID, provirus) %>%
    tally() %>%
    filter(provirus == "Yes")

number_provirus.filt$provirus_num <- number_provirus.filt$n

number_check.filt <- combo.filt %>%
    group_by(SampleID, checkv_quality) %>%
    tally() %>%
    spread(checkv_quality, n) %>%
    replace(is.na(.), 0)

checkv_cdhit_tax.meta.samplesum.combo.filt <- checkv_cdhit_tax.meta.samplesum.filt %>%
    left_join(dplyr::select(number_provirus.filt, SampleID, provirus_num)) %>%
    left_join(number_check.filt) %>%
    left_join(dplyr::select(sampledata2, SampleID, AreaSample,
        Site.name...8, `Year of collection`, Rocks_v2, Latitude,
        Longitude)) %>%
    replace(is.na(.), 0)
```

    ## Joining, by = "SampleID"
    ## Joining, by = "SampleID"
    ## Joining, by = "SampleID"

``` r
# contig length
summary(checkv_cdhit_tax.meta.samplesum.combo.filt$contig_length_avg)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   13297   20293   22208   22610   24234   43442

``` r
# avg gene count per seq
summary(checkv_cdhit_tax.meta.samplesum.combo.filt$gene_count_avg)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   18.21   25.52   27.71   28.09   30.25   42.78

``` r
# avg number viral seq
summary(checkv_cdhit_tax.meta.samplesum.combo.filt$num_seqs)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    4.00   39.00   66.50   80.41  111.00  250.00

``` r
# total number of viral seq
sum(checkv_cdhit_tax.meta.samplesum.combo.filt$num_seqs)
```

    ## [1] 14796

``` r
# avg number provirus
summary(checkv_cdhit_tax.meta.samplesum.combo.filt$provirus_num)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    0.00    5.75   11.00   14.65   19.00   68.00

``` r
# total number of provirus
sum(checkv_cdhit_tax.meta.samplesum.combo.filt$provirus_num)
```

    ## [1] 2695

``` r
# total num complete genomes
sum(checkv_cdhit_tax.meta.samplesum.combo.filt$Complete)
```

    ## [1] 76

``` r
# total number high-quality
sum(checkv_cdhit_tax.meta.samplesum.combo.filt$`High-quality`)
```

    ## [1] 341

``` r
# total number medium-quality
sum(checkv_cdhit_tax.meta.samplesum.combo.filt$`Medium-quality`)
```

    ## [1] 1539

``` r
write.csv(checkv_cdhit_tax.meta.samplesum.combo.filt, "results/sample_sum_info.filt.csv")
```
