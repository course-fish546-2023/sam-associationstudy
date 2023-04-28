# Libraries
# data.table is a lot faster than traditional data frames
library(data.table)
library(tidyverse)

# List of accession numbers that we're interested in
accession_list <- c("GCF_000006765.1","GCF_000011965.2","GCF_000157115.2","GCF_000012345.1")

new_accession_list <- c()

# Paste the accession IDs together with a URL-safe comma separating them
accession_string <- paste(new_accession_list, collapse = "%2C")

# Paste the string into the API's URL
# Request the genome GFF with annotations
accession_url <- paste0("https://api.ncbi.nlm.nih.gov/datasets/",
                        "v2alpha/genome/accession/", accession_string,
                        "/download?include_annotation_type=GENOME_GFF&",
                        "filename=string_file.zip")
# You can change where this is saved with the "destfile" argument, righ
# now it's just the working directory
download.file(accession_url, destfile = "./data/temp.zip")
unzip("./data/temp.zip",exdir="./data")


# Loop over each entry in accession list, read in the table from the
# associated folder, and extract the gene IDs

gene_presence_list <- lapply(accession_list, function(accession_number_i){
  filename <- paste0("./data/ncbi_dataset/data/", accession_number_i, "/genomic.gff")
  # Apparently the files have different length headers so we need to skip them dynamically
  # from https://stackoverflow.com/q/18920777
  # This makes the reader much slower but I can't figure out how to avoid it otherwise
  lineskip_command <- paste("grep -v '^#'", filename)
  # fread is data.table's fast table reader function
  # we also index for column 3 equaling "gene" at the end here
  tab_i <- fread(cmd=lineskip_command, sep = "\t", fill = TRUE, )[V3=="gene"]
  
  
  # Approach number 1 - downloading the GeneID name
  # This regex is a lookbehind (?<=) and grabs all digits (\\d) following 
  # the lookbehind match
  # gene_ids <- as.numeric(str_extract(tab_i$V9, pattern = "(?<=Dbxref=GeneID:)\\d+"))
  
  # Approach number 2 - using the gene name technique
  gene_names <- str_extract(tab_i$V9, pattern="(?<=gene=)[^;]*(?=;)")
  
  # Remove all but 1 NA
  gene_names <- unique(gene_names, incomparables=FALSE, fromLast=FALSE, by=seq_along(gene_names))
  
  # Rename to avoid name collision on merging
  output_table <- data.table(1, gene_names)
  names(output_table) <- c(accession_number_i, "gene_names")
  return(output_table)
})

# Reduce using the merge function, keeping all rows, accumulating genes each time
genes_merged <- purrr::reduce(gene_presence_list, dplyr::full_join, by="gene_names")


# Alternative way leaving it as a dataframe
# genes_merged <- Reduce(function(x, y){merge(x, y, all=TRUE)}, gene_presence_list)

# Convert into matrix
genes_merged[is.na(genes_merged)] <- 0

# Convert into a matrix
genes_merged <- column_to_rownames(genes_merged,var="gene_names")

# Save genes_merged as a csv
write.csv(genes_merged, "./code/data/merged_genes.csv", row.names=TRUE)

genes_merged <- read.csv("./code/data/merged_genes.csv",row.names = 1, header= TRUE)

# From the list, choose a reference gene
ref_gene <- as.numeric(genes_merged["aat",])

# Calculate Jaccard Similarity
gene_matches <- colSums(t(genes_merged)+ref_gene==2)
gene_totalhits <- colSums(t(genes_merged)+ref_gene>0)
fave_genes <- sort(gene_matches/gene_totalhits, decreasing = TRUE)

# Check the confusion matrix of the top 10 best matched genes
test <- lapply(names(head(fave_genes, 10)), function(gene_name){
  table(as.numeric(genes_merged["aat",]), 
        as.numeric(genes_merged[gene_name,]))
})

# plot the ordered data based on similarity score
x_steps = seq(1, length(fave_genes)-1, 1)
y_steps = fave_genes[-1]
df = data.frame(x=x_steps,y=y_steps)
ggplot(data=df,mapping = aes(x, y))+geom_point()
xlab("Rank Order")
ylab("Similarity Score")
