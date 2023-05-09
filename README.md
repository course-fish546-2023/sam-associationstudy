# Bacterial Reference Genome Association Study

##Database Generation
Once the genomes are downloaded, I will want to annotate them (or download the annotated version) and generate a list of all possible genes and then move that data into a sparse matrix.
 - download the ncbi files (ncbi-genome-download) (ideally download the annotations)
 - repeated full joins across the entire setup to get a large presence/absence matrix
      - most worried about this step and making sure that I can match up same proteins
 - convert matrix to a Sparsematrix.


##Database Search
From the sparse matrix, I need to choose a target gene, and separate out into two groups, those with and those without the target. I want to sum all the genes, and divide by total number of organisms in both groups, then look at the log fold or difference between the two groups. Sort the list by the largest difference/log fold and plot.

 - User input for target gene, take that column in the Sparsematrix, use that column as a filter on the matrix, and invert and use as a filter on the other matrix
 - count number of rows
 - sum columns in each matrix, then divide by the total number of rows
 - subtract or divide presence_matrix from absence_matrix
 - sort across the entire space, and see which genes are the most correlated

Current Progress:
Week 1+2: Not in Class

Week 3: Generated simulated data and created system to quickly generate a Jaccard score from a presence/absence matrix (under simulated_data_test.R)

Week 4: Worked on devloping a tool to download reference genomes and create a presence/absence matrix - used NCBI API to download one genome at a time (under download_script.R)

Week 5: Continued working on a download tool, eventually got NCBI datasets to work. Finalized the presence/absence matrix code. Prepared presentation. (under 01-ncbi_download.Rmd)

Week 6: Created a 17,022 genes x 17,056 genomes presence/absence matrix using premade code, ran some initial tests on genes of interest and saw expected correlations (genes within same operon were correlated, genes associated with central carbon metabolism were correlated with each other) (01-ncbi_download.Rmd into 02-full_matrix_test.Rmd)

Week 7: Generated a 17,022 x 17,022 correlation matrix with the Jaccard Similarity score of each gene-gene pair (currently running) (02-full_matrix_test.Rmd) - could also filter out all genes that appear in fewer than X genomes to speed up the calculation.

Week 8: Identify the strongest correlations and research those genes

Week 9: Compare the distribution of probabilities from each gene to each other to search for odd distributions

Week 10: Test other metrics instead of Jaccard and see if specific genes have similar correlations

