#### Finding 1000 random genes to use for hypergraph control ####

library(biomaRt)

#Connect to human database:
human <- useMart(biomart="ensembl", dataset = "hsapiens_gene_ensembl", verbose=TRUE, host="https://feb2021.archive.ensembl.org")

#retrieve all human gene symbols:
all_genes <- getBM(attributes = "hgnc_symbol", mart = human)
all_genes <- all_genes$hgnc_symbol[all_genes$hgnc_symbol != ""]

#set seed to allow me to come back to this if needed.
set.seed(1)
random_genes <- sample(all_genes, 1000)

#Put these genes into a dataframe
random_genes <- data.frame(gene = random_genes)

#save as csv
write.csv(random_genes, paste0(get_wd(),"random_genes.csv"))

