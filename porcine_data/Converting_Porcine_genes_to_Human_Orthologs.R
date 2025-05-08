#### CONVERT PORCINE DEGs TO HUMAN ####
setwd("C:/Users/anna-/OneDrive/Documents/Anna/MSc. Reproduction and Pregnancy/Research Project 2/pig/galois_results")
data<- read.csv("all_galois_genes_cluster_5.csv")

#Extract DEG names from the rownames
porcine_gene_names <- data$Gene

#Connect to Ensembl datasets (use old archive to avoid server error)
human <- useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose=TRUE, host="https://feb2021.archive.ensembl.org")
porcine <- useMart(biomart="ensembl", dataset = "sscrofa_gene_ensembl", verbose=TRUE, host="https://feb2021.archive.ensembl.org")

#Set up a query to convert the gene names using the get LDS function:
genesV2 <- getLDS(attributes = c("external_gene_name"),
                  filters ="external_gene_name",
                  values = porcine_gene_names,
                  mart = porcine,
                  attributesL = c("hgnc_symbol"),
                  martL = human,
                  uniqueRows = TRUE)


# Match row names (GeneID) in dge_output with GeneID in genesV2
matched_genes <- intersect(data$Gene, genesV2$Gene.name) #These are the pig genes that have a human ortholog


# Add HGNC.symbol column to dge_output based on matching row names
data$HGNC.symbol <- NA  # make a new column
match_indices <- match(data$Gene, genesV2$Gene.name)
data$HGNC.symbol <- genesV2$HGNC.symbol[match_indices]
data<- data[!is.na(data$HGNC.symbol), ]#omit rows without human ortholog

write.csv(data, "Human_orthologs_of_all_galois_cluster_5.csv")
