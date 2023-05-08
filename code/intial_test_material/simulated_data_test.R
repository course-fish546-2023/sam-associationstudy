# expand.grid produces permutations of each of the inputs and returns a data
# fram with each permutation in a column
perm_letters <- expand.grid(c(rep("A", 9), rep("B", 4), LETTERS), LETTERS, LETTERS)
# paste0 (and the more general paste) combines the three columns into a single
# character vector
list_of_genes <- paste0(perm_letters$Var1, perm_letters$Var2, perm_letters$Var3)

# Same method as above
perm_genome_names <- head(expand.grid(LETTERS, 1:40), 1000)
genome_names <- paste0(perm_genome_names$Var1, perm_genome_names$Var2)

# Loop over each genome to be created and produce a data frame
genome_list <- lapply(genome_names, function(genome_name){
  # work from the inside out on this one, it's deeply nested
  # sample(100:200) to get a random number of genes in a genome
  # sample(list_of_genes) to figure out exactly which genes those are
  # Fill in the "present" column with 1s since they're all found in this case
  v <- data.frame(genes=unique(sample(list_of_genes, size = sample(100:200, size = 1))), genome_name=1)
  names(v) <- c("genes", genome_name)
  v
})
# TODO: convert to library(data.table)s instead of data.frame and switch to merge() instead
# of full_join because it should be much faster?
# https://tysonbarrett.com/jekyll/update/2019/10/11/speed_of_joins/

# Checking that the start letter distribution matches the skew we added above
sapply(LETTERS, function(i)sum(grepl(paste0("^", i), genome_list[[1]]$genes)))


# Timing for full join reduce since it took forever first try (1+ minute)
dplyr::full_join(genome_list[[1]], genome_list[[2]], by="genes")
system.time(reduce(genome_list[1:100], dplyr::full_join, by="genes"))
nrow_timepoints <- c(seq(20, 100, by=10), seq(100, 300, by=50))
timing <- sapply(nrow_timepoints, function(n_genomes){
  system.time(reduce(head(genome_list, n_genomes), full_join, by="genes"))
})
plot(nrow_timepoints, timing["elapsed",])
lm_data <- data.frame(nrow_timepoints, time_elapsed=timing["elapsed",])
model_fit <- lm(formula = time_elapsed~poly(nrow_timepoints, 2), data = lm_data)
summary(model_fit)
pred_data <- data.frame(
  nrow_timepoints=1:1000,
  pred=predict(model_fit, newdata = data.frame(nrow_timepoints=1:1000))
)
plot(pred_data)
points(lm_data, col="red", pch=19)
plot(lm_data$nrow_timepoints, sqrt(lm_data$time_elapsed))

# Return to full join now that we know how long it'll take
# Reduce takes a list of arguments and repeatedly applies an aggregating function
# purrr's reduce is able to handle the additional by="genes" argument while
# the base R Reduce couldn't
# TODO: implement O(nlog(n)) version with binary tree joins??
gene_genome_matrix_raw <- purrr::reduce(genome_list, dplyr::full_join, by="genes")
gene_genome_matrix <- gene_genome_matrix_raw
gene_genome_matrix[is.na(gene_genome_matrix)] <- 0

gene_genome_matrix[1:10, 1:10]

gene_genome_matrix <- gene_genome_matrix[order(gene_genome_matrix$genes),]
rownames(gene_genome_matrix) <- gene_genome_matrix$genes
gene_genome_matrix$genes <- NULL

ref_gene <- as.numeric(gene_genome_matrix["AAA",])
# Calculate Manhattan distance (taxicab distance)
# sweep is a quirky function - applies FUN element-wise to each tuple of
# row[i] and ref_gene[i]
gene_distance <- rowSums(abs(sweep(gene_genome_matrix, 2, STATS = ref_gene, FUN = `-`)))
hist(gene_distance)
which(gene_distance==0)
head(sort(gene_distance))
table(as.numeric(gene_genome_matrix["AAA",]), 
      as.numeric(gene_genome_matrix["YVR",]))


# Jaccard on the full dataset takes forever (5+ minute)
# jacdist <- vegan::vegdist(gene_genome_matrix, method = "jaccard")
# Calculate jaccard for a single reference gene
# sweep output == 2 when both are 1, i.e. a good match to be rewarded
# sweep output >= 1 whenever either ref genome
gene_matches <- rowSums(sweep(gene_genome_matrix, 2, STATS = ref_gene, FUN = `+`)==2)
gene_totalhits <- rowSums(sweep(gene_genome_matrix, 2, STATS = ref_gene, FUN = `+`)>0)
fave_genes <- sort(gene_matches/gene_totalhits, decreasing = TRUE)

# ref_gene is automatically replicated until no more elements are needed
gene_matches <- colSums(t(gene_genome_matrix)+ref_gene==2)
gene_totalhits <- colSums(t(gene_genome_matrix)+ref_gene>0)

# Check on output distibution
test <- lapply(names(head(fave_genes, 10)), function(gene_name){
  table(as.numeric(gene_genome_matrix["narG",]), 
        as.numeric(gene_genome_matrix[gene_name,]))
})
