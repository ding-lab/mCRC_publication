library(tidyverse)
message("useage: RemoveMicrobiomeGenes(obj)")

GetMicrobiomeGene <- function(obj){
	key_words = c('Escherichia', 'papillomavirus', 'Fusobacterium', 'Escherichia', 'Bacteroides')
	rownames(obj) %>% str_subset(paste(key_words, collapse = '|')) 
}

NonMicrobiomeGenes <- function(obj){
	message('Removing MicrobiomeGenes: ', str_c(GetMicrobiomeGene(obj), collapse = '\n'))
	setdiff(rownames(obj),GetMicrobiomeGene(obj))
}

source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Pipeline/SharedFunctions/Seurat/function_SubsetFeatures.R')
# SubsetFeatures <- function(object, features=NULL){
# 	# From DietSeurat
# 	if(is.null(features)){
# 		return(object)
# 	}
# 	for (assay in Assays(object = object)) {
#     if (!is.null(x = features)) {
#       features.assay <- intersect(
#         x = features,
#         y = rownames(x = object[[assay]])
#       )
#       if (!length(x = features.assay)) {
#         warning(message = paste0(
#           'No features found in assay ',
#           sQuote(x = assay),
#           ', removing...'
#         ))
#         object[[assay]] <- NULL
#         next
#       }
#       suppressWarnings(object[[assay]] <- subset(x = object[[assay]], features = features.assay))
#     }
#   }
#   return(object)
# }

RemoveMicrobiomeGenes <- function(obj){
	# Remove Microbiome genes
	genes_keep <- NonMicrobiomeGenes(obj)
	obj <- SubsetFeatures(obj, genes_keep)
	return(obj)
}

######################################################################################
# 2025-05-01 Mutation Probes
# Also remove mutation probes: pattern p.
GetMutationProbes <- function(obj){
  key_words = c('p[.]')
  rownames(obj) %>% str_subset(paste(key_words, collapse = '|')) 
}

# Remove mutation probes
NonMutationProbes <- function(obj){
  message('Removing MutationProbes: ', str_c(GetMutationProbes(obj), collapse = '\n'))
  setdiff(rownames(obj),GetMutationProbes(obj))
}

RemoveMutationProbes <- function(obj){
  # Remove Microbiome genes
  genes_keep <- NonMutationProbes(obj)
  obj <- SubsetFeatures(obj, genes_keep)
  return(obj)
}