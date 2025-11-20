# 2025-03-05. add min_transcripts option
# 2025-02-25. Updated using function_SlimSeurat_v2.R that allow using all genes instead of shared only
# Also add option to remove microbiome genes. Default is to remove microbiome genes

# conda activate seurat5_v2
# from: /diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/2_SlimMergedSeurat/script/1_slim_seurat_merge.R
# in v2, change use optparse to load the input sample sheet

library(optparse)

option_list = list(
	make_option(c("-i", "--sample_table"), type="character", default=NULL, help="Input sample sheet. Need to contain 3 columns: SectionID, CasetteID, FilePath"),
	make_option(c("-o", "--out_dir"), type="character", default=NULL, help="Output directory"),
	make_option(c("-n", "--analysis_name"), type="character", default=NULL, help="Output name"),
	make_option(c("-t", "--min_transcripts"), type="integer", default=5, help="Remove cells with less than min_transcripts. Default is 5"),
	make_option(c("-s", "--shared_genes_only"), type="logical", default=TRUE, help="Use shared genes only"),
	make_option(c("-m", "--remove_microbiome"), type="logical", default=TRUE, help="Remove microbiome genes"),
	make_option(c("-f", "--features_plot"), type="character", default="EPCAM,AR,KLK2,KLK3,FOLH1,PCA3,AMACR,ERG,TMPRSS2", help="Features to plot")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(dput(opt))

library(tidyverse)
library(Seurat)
library(qs)
library(googlesheets4)

# replace message with message with timestamp
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Pipeline/Timestamp/message_with_timestamp_v1.R')

# Assign parameters
out_dir = opt$out_dir
sheet_use = read_tsv(opt$sample_table) %>% setNames(c('SectionID','CasetteID','FilePath'))
name = opt$analysis_name

# out
#out_dir ="/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Prostate/Analysis/SlimXenium/1_createObject/out"


# Xenium mouse sheet  processing------------------------------------------
# gs4_deauth()
# sheet =read_sheet('https://docs.google.com/spreadsheets/d/1EPo3LxnrMsdoBwAFV3hSNku7vii-ONtpyx1oJQdbVEg/edit?gid=0#gid=0', sheet = 'FFPE processed')
# sheet_mouse = sheet %>% 
# 	filter(
# 		Species =='Human',
# 		SampleType == 'PCa',
# 		file.exists(`Output file path`),
# 		)
# sheet_use = sheet_mouse[,c('Section ID', 'Casette ID', 'Output file path')] %>%  # Casette ID as batch
# 	setNames(c('SectionID','CasetteID','FilePath'))
# dim(sheet_use)



## Create object ------------------------------------------
# functions: Timer, Slim seurat, Shared genes, getSize
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/2_SlimMergedSeurat/script/src/function_SlimSeurat_v3.R')
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/1_LeanSeurat/script/src/function_Timer.R')
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/1_LeanSeurat/script/src/function_sharedGenes.R')
source('/diskmnt/Projects/Users/simonmo/projects/MouseKidney/Revision4/5_Xenium/1_LeanSeurat/script/src/function_getSize.R')

# Load samples
Timer()
# Oneliner to create object 
sheet_use_formated = sheet_use %>% dplyr::rename(SampleID = SectionID, FolderPath = FilePath) 
obj_merged = MakeMergeSlimFromSampleTableLayer(sheet_use_formated, opt$shared_genes_only, min_ncount = opt$min_transcripts)
obj_merged = AddSampleLevelMeta(obj_merged, 
	meta_df = sheet_use_formated, 
	meta_sample_column = 'SampleID', 
	obj_sample_column = 'orig.ident')
Timer()


getSize(obj_merged)

# Remove microbiome genes
if(opt$remove_microbiome){
	source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Pipeline/XeniumProbe/RemoveMicrobiomeProbes_v1.R')
	obj_merged = RemoveMicrobiomeGenes(obj_merged)
}

# save 
qsave(obj_merged, str_glue('{out_dir}/1_xenium_slim_{name}.qs', nthreads = 30))

# Log Noramlize Data
obj_merged_data = NormalizeData(obj_merged)
getSize(obj_merged_data) # Obj Size incrased to 9.9GB. On Disk 2.5 GB

qsave(obj_merged_data, str_glue('{out_dir}/1_xenium_slim_{name}_normalized.qs', nthreads = 30))



# Plot with custom function
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/Tool/Seurat/Plot/ImagePlotFaster_v1.R')
library(ImageBinPlot)
#library(viridis)

pdf(str_glue('{out_dir}/2a_nCount_spatial.pdf'), width = 30, height = 30)
	print(ImagePlotFasterAllFOV(obj_merged, features = 'nCount_Xenium', max.cutoff="q85"))
dev.off()


pdf(str_glue('{out_dir}/2b_nFeature_spatial.pdf'), width = 30, height = 30)
	print(ImagePlotFasterAllFOV(obj_merged, features = 'nFeature_Xenium', max.cutoff="q85"))
dev.off()

# Plot features
features_plt = str_split(opt$features_plot, ',')[[1]] %>% intersect(rownames(obj_merged));print(features_plt)
for(feat in features_plt){
	pdf(str_glue('{out_dir}/3_features_{feat}_spatial.pdf'), width = 30, height = 30)
		print(ImagePlotFasterAllFOV(obj_merged, features = feat, max.cutoff="q85"))
		#print(ImageBinPlotFOVs(obj_merged, feature = feat, type = "expression", bin_size = 10, max_quantile = "q75",
  #same_scale = TRUE) & theme_black)
	dev.off()
}
