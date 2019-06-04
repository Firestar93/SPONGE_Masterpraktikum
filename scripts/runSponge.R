library("data.table")
library("tidyr")
library("gridExtra")
library("ggrepel")
library("visNetwork")
library("ggplot2")
library("SPONGE")

#gene expression
genes <- fread("C:\\Users\\Elly\\Dropbox\\Uni\\Master\\2._Semester\\NEAP\\data\\short_genes_expr.tsv", header = T)
genes <- data.frame(genes[,-1], row.names=genes$V1)
gene_matrix <- as.matrix(genes)

#mirna expression
mirna <- fread("C:\\Users\\Elly\\Dropbox\\Uni\\Master\\2._Semester\\NEAP\\data\\short_miRNA_expr.tsv", header = T)
mirna <- data.frame(mirna[,-1], row.names=mirna$V1)
mirna_matrix <- as.matrix(mirna)

load("C:\\Users\\Elly\\Dropbox\\Uni\\Master\\2._Semester\\NEAP\\data\\targetscan.RData")

#######################################################################
#apply step A of the SPONGE workflow
#######################################################################
genes_miRNA_candidates <- sponge_gene_miRNA_interaction_filter(
         gene_expr = gene_matrix,
         mir_expr = mirna_matrix,
         mir_predicted_targets = targetscan_ensg,
         coefficient.threshold = 0)

#######################################################################
#apply step B of the SPONGE workflow
#######################################################################
ceRNA_interactions <- sponge(gene_expr = gene_expr,
                             mir_expr = mir_expr,
                             mir_interactions = genes_miRNA_candidates)

#######################################################################
#apply step C of the SPONGE workflow
#######################################################################


#number_of_datasets default: 1e+05 (sollte damit laufen) --> aus Geschw.gründen hier am Bsp geringer gemacht
mscor_null_model <- sponge_build_null_model(number_of_samples = nrow(gene_expr))

#plot the result of these simulations to see how dist is affected by gene-gene correlation and number of miRNAs
sponge_plot_simulation_results

#infer p-values
#method adds two additional columns to the ceRNA results, namely p.val and p.adj
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions, 
                                                   null_model = mscor_null_model)

#######################################################################
#apply step D of the SPONGE workflow
#######################################################################

#decide on an FDR cutoff for selecting significant ceRNA interactions, e.g. FDR < 0.01
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < 0.2),]

#plot resulting network
sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)
