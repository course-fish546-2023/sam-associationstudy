Overalll Scope of the Project

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


## File Organization
### initial_test_materials 
shows the first test of simulated code, and also a theoretical download pipeline using the NCBI API to download genomes

### final_pipeline
01-ncbi_download - Rmd document showing the workflow to generate the full presence/absence matrix for all reference genomes. Saves to a CSV labelled merged_genes_full.csv in /data. Matrix approximately 17,000genes x 17,000genomes
02-full_matrix_test - Rmd document to take the matrix and find the Jaccard similarity score of every gene relative to a reference, and sort.
