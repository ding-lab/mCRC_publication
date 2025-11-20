# written by Austin Southard-Smith
# some modification/feedback provided by Evan Chien-Wei Peng in v6
# 2023-10-11
# based on https://satijalab.org/seurat/articles/spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ
# addresses error mentioned https://github.com/satijalab/seurat/issues/7491
# update 2024-07-01: added in v5 of the script on 2024-07-01 to generate the `transcripts.csv.gz`file that is missing from the output of the Xenium onboard Analysis v3.0 software. I assume that Seurat will be made compatible with this in the future but right now they have not added any sort of translation functionality.
# update 2024-12-31: added in v6 of the script:
        # Output plots are now rasterized to (300 dpi for dotplots and 600 dpi for spatial image plots). This decreases the file size for plots of larger datasets
        # Added the Xenium.with.snvs.microbiome and Xenium.with.microbiome assays and their respective fovs. With v1 of the Xenium 5k panel some probes were generated to target transcripts expressed by the microbiome. Like the competing Variant/Reference probes targeting snvs and indels these are now removed from the Default Xenium assay present in the object and added as part of a separate assay so they are not normalized with the rest of the expression data during processing and analysis.
# update 2025-2-1: add flag to remove unneed cells  
# update 2025-4-17: add flag to keep a subset of cells         
library(Seurat)
library(tidyverse)
library(future)
library(optparse)
library(magrittr)
library(ggrastr)
set.seed(1234)
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL,
                help="path to data folder (e.g. xenium output folder: /diskmnt/primary/Xenium/data/20230803__195918__20230803_firstsamples/output-XETG00122__0010660__HT270P1-S1H1A1Fp1Us1_1__20230803__200008/)",
                metavar="character"),
    make_option(c('-s', "--sampleID"),
                type="character",
                default=NULL,
                help="SampleID that output file names will be based on",
                metavar="character"),
    make_option(c("-o", "--output"),
                type="character",
                default="./",
                help="output folder path",
                metavar="character"),
    make_option(c("--mols.qv.threshold"),
                type="integer",
                default=20,
                help="threshold for the qv score cutoff threshold. Transcripts above this will be included.",
                metavar="integer"),
    make_option(c("--pc_num"),
                type="integer",
                default=30,
                help="number of principal components to use",
                metavar="integer"),
    make_option(c("--counts_cutoff"),
                type="integer",
                default=0,
                help="retain only cells with nCounts greater than the value specified by the counts_cutoff",
                metavar="integer"),
    make_option(c("--features_cutoff"),
                type="integer",
                default=0,
                help="retain only cells with nCounts greater than the value specified by the counts_cutoff",
                metavar="integer"),
    make_option(c("--plot_features"),
                type='character',
                default=NA,
                help = "comma separated list of features to plot when the object is generated (for pancreas this looks like 'KRT7,AMY2A,TFF2,CFTR,EPCAM,MKI67,TOP2A,PPARG,PECAM1,VWF,MS4A1,CD3D,CD3E,CD4,MYC,INS,GCG,C7,PDGFRA,SRPX,MAMDC2,VCAN,MEDAG,MFAP5,GPRC5A,CAV1,PTPRC,IL7R,TMC5,MET,EHF,CAPN8,CD163,MS4A6A,PDGFRB,MYH11,KCNMA1'",
                metavar="character"),
    make_option(c("--with_snvs"),
                action="store_true",
                default=FALSE,
                help="Use this flag when the input xenium data also contains SNVs present. By default if this flag is used then the assay will only remove those Xenium probes with the pattern '_WT' or '_ALT' in the probe name. To specify other snv probes to remove use the `--snv_probes` flag. These snv probes need to be filtered out prior to any normalization or subsequent marker analysis. When used this flag results in the resulting object having two raw assays one is the 'Xenium' assay and the other is the 'Xenium.with.snvs' assay. The 'Xenium' assay has no 'snv' target probes in it and is used for subsequent SCT normalization and can be used for future non-mutation related analysis. There are then two whole slide FOVs generated. One is the normal 'fov' which is intended for use with the 'Xenium' and 'SCT' assay. The other is the 'fov.with.snvs' image which is intended for use with the 'Xenium.with.snvs' assay. In order to visualize the localization of the SNV probes on the section you need to use the 'fov.with.snvs' image. If this option is not specified then only the 'Xenium' assay and 'fov' image will be present in the resulting seurat object. If this flag is used with the --with_microbiome flag then an additional raw assay 'Xenium.with.snvs.microbiome' and full section FOV 'fov.with.snvs.microbiome' will be generated."),
    make_option(c("--with_microbiome"),
                action="store_true",
                default=FALSE,
                help="Use if the gene panel includes both SNVS and microbiome probes. By default if this flag is used then the assay will only remove those Xenium probes with the pattern 'Bacteroides', 'Escherichia', 'Fusobacterium', and 'papillomavirus' in the probe name. To specify other snv probes to remove use the `--microbiome_probes` flag. These snv probes need to be filtered out prior to any normalization or subsequent marker analysis. When used this flag results in the resulting object having two raw assays one is the 'Xenium' assay and the other is the 'Xenium.with.microbiome' assay. The 'Xenium' assay has no 'microbiome' target probes in it and is used for subsequent SCT normalization and can be used for future non-mutation related analysis. There are then two whole slide FOVs generated. One is the normal 'fov' which is intended for use with the 'Xenium' and 'SCT' assay. The other is the 'fov.with.microbiome' image which is intended for use with the 'Xenium.with.microbiome' assay. In order to visualize the localization of the SNV probes on the section you need to use the 'fov.with.microbiome' image. If this option is not specified then only the 'Xenium' assay and 'fov' image will be present in the resulting seurat object. If this flag is used with the --with_snvs flag then an additional raw assay 'Xenium.with.snvs.microbiome' and full section FOV 'fov.with.snvs.microbiome' will be generated."),
    make_option(c("--snv_probes"),
                type='character',
                default=NA,
                help="Only use this flag in combination with the `--with_snvs`. This flag does not do anything if it is used and `--with_snvs` is not also used. Use this flag followed by the path to a file containing a list of SNV probes that need to be removed (all probes to remove need to be present in this file) from the seurat object when the object is being generated. (e.g. '--with_snvs --snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_probe_names.tsv'). Any probe/gene specified in this file will be removed from the 'Xenium` assay and 'fov' image. These probes will still be retained and accessible from the 'Xenium.with.snvs' assay and 'fov.with.snvs' image."),
    make_option(c("--microbiome_probes"),
                type='character',
                default=NA,
                help="Only use this flag in combination with the `--with_snvs_microbiome`. This flag does not do anything if it is used and `--with_microbiome` is not also used. Use this flag followed by the path to a file containing a list of microbiome probes that need to be removed (all probes to remove need to be present in this file) from the seurat object when the object is being generated. (e.g. '--with_microbiome --microbiome_probes /diskmnt/Projects/Users/austins2/tools/Xenium_human_panel_5k_v1_DHYHMW_microbiome_probe_names.tsv'). Any probe/gene specified in this file will be removed from the 'Xenium` assay and 'fov' image. These probes will still be retained and accessible from the 'Xenium.with.microbiome' assay and 'fov.with.microbiome' image."),
    make_option(c("--cell_removal_csv"),
                type='character',
                default=NULL,
                help="input a csv file indicate the cells to remove"),
    make_option(c("--cell_keep_csv"),
                type='character',
                default=NULL,
                help="input a csv file indicate the cells to keep")    
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
    print_help(opt_parser)
    stop("Path to data is required (--input).n", call.=FALSE)
}
if (is.null(opt$sampleID)){
    print_help(opt_parser)
    stop("sampleID is required (--sampleID).n", call.=FALSE)
}

input = opt$input
sampleID = opt$sampleID
out_path = paste0(opt$output,"/",sampleID,"/")
dir.create(out_path)
num_pcs = opt$pc_num
cell_removal_csv = opt$cell_removal_csv
cell_keep_csv = opt$cell_keep_csv

#added in v5 of the script on 2024-07-01 to generate the `transcripts.csv.gz`file that is missing from the output of the Xenium onboard Analysis v3.0 software.
PATH = paste0(input,"/", "transcripts.parquet")
OUTPUT <- gsub('\\.parquet$', '.csv', PATH)
if (!(file.exists(paste0(OUTPUT,".gz")))) {
    library(arrow)
    CHUNK_SIZE <- 1e6
    parquet_file <- arrow::read_parquet(PATH, as_data_frame = FALSE)
    start <- 0
    while(start < parquet_file$num_rows) {
        end <- min(start + CHUNK_SIZE, parquet_file$num_rows)
        chunk <- as.data.frame(parquet_file$Slice(start, end - start))
        data.table::fwrite(chunk, OUTPUT, append = start != 0)
        start <- end
    }
    if(require('R.utils', quietly = TRUE)) {
        R.utils::gzip(OUTPUT)
    }
}

probe_hyphens <- function(probe_list) {
    probe_list <- gsub(":", "-", probe_list)
    probe_list <- gsub("\\.", "-", probe_list)
    probe_list <- gsub("_", "-", probe_list)
    probe_list <- gsub("\\+", "-plus-", probe_list)
    return(probe_list)
}

create_assay_and_fov <- function(counts, pixels, segmentations, list.to.remove, assay_name) {
    counts_working <- counts
    counts.no.remove <- counts_working[(!rownames(counts_working) %in% c(list.to.remove)),]
    rm(counts_working)
    assay.no.remove <- CreateAssayObject(counts = counts.no.remove)
    pixels_working <- pixels
    # print(dim(pixels))
    # print(head(pixels))
    pixels.indices.remove <- which(pixels_working[,"gene"] %in% list.to.remove, arr.ind = T)
    # print(head(pixels.indices.remove))
    # print(dim(pixels.indices.remove))
    pixels.no.remove <- pixels_working[(!rownames(pixels_working) %in% c(pixels.indices.remove)),]
    # print(dim(pixels.no.remove))
    # print(head(pixels.no.remove))
    rm(pixels_working)
    fov.no.remove <- CreateFOV(
        coords = segmentations,
        type = c("segmentation", "centroids"),
        molecules = pixels.no.remove,
        assay = assay_name
    )
    rm(counts.no.remove, pixels.no.remove)
    gc()
    return(list(assay.no.remove, fov.no.remove))
}

plan("multicore", workers = 30)
# The below errors out because it is reading in the Cell IDs incorrectly per https://github.com/satijalab/seurat/issues/7491
# xenium.obj <- LoadXenium(data.dir=input, fov = "fov", assay = "Xenium") #this might be fixed in Seurat5.
xenium.list <- ReadXenium(data.dir=input, outs = c("matrix", "microns"), type = c("centroids", "segmentations"), mols.qv.threshold = opt$mols.qv.threshold)
plan("sequential")

# fixing input matrix names
rownames(xenium.list$matrix$`Gene Expression`) <- gsub(":","-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub(":","-",xenium.list$microns$gene)
rownames(xenium.list$matrix$`Gene Expression`) <- gsub("\\.","-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub("\\.","-",xenium.list$microns$gene)
rownames(xenium.list$matrix$`Gene Expression`) <- gsub("_","-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub("_","-",xenium.list$microns$gene)
rownames(xenium.list$matrix$`Gene Expression`) <- gsub("\\+","-plus-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub("\\+","-plus-",xenium.list$microns$gene)
# extracting the list of features that variant/reference and microbiome probes will be filtered from
targets <- rownames(xenium.list$matrix[["Gene Expression"]])

if (opt$with_snvs && opt$with_microbiome) { 
    # compiling the list of snvs features that will be removed
    if (is.na(opt$snv_probes)) {
        wt_targets <- targets[grepl("-WT",targets)]
        alt_targets <- targets[grepl("-ALT",targets)]
        snv.row.names.remove <- c(wt_targets, alt_targets)
    } else {
        snv.row.names.remove <- readLines(opt$snv_probes)
    }
    snv.row.names.remove <- probe_hyphens(snv.row.names.remove)
    print("Removing the following SNV probes from the default Xenium matrix")
    print(snv.row.names.remove)
    
    if (is.na(opt$microbiome_probes)) {
        bacteroides_targets  <- targets[grepl("Bacteroides", targets)]
        escherichia_targets  <- targets[grepl("Escherichia", targets)]
        fusobacterium_targets  <- targets[grepl("Fusobacterium", targets)]
        hpv_targets  <- targets[grepl("papillomavirus", targets)]
        microbiome.row.names.remove <- c(bacteroides_targets, escherichia_targets, fusobacterium_targets, hpv_targets)
    } else {
        microbiome.row.names.remove <- readLines(opt$microbiome_probes)
    }
    microbiome.row.names.remove <- probe_hyphens(microbiome.row.names.remove)
    print("Removing the following Microbiome probes from the default Xenium matrix")
    print(microbiome.row.names.remove)
    
    all.row.names.remove <- c(snv.row.names.remove, microbiome.row.names.remove)
    # each assay and FOV will be associated with the same segmentation data since we aren't subsetting cells yet.
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    # generating the object with the initial assay with everything
    assay <- "Xenium.with.snvs.microbiome"
    fov <- "fov.with.snvs.microbiome" #If you are using Seurat v4.3.0 do not rename this as it will result in you being unable to plot molecules.
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay) #In Seurat 5 this creates a v5 assay.
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- all.coords
    # generating the assay and fov with genes and variant/reference probes
    assay <- "Xenium.with.snvs"
    fov <- "fov.with.snvs"
    with.snvs.list <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                           xenium.list$microns,
                                           segmentations.data,
                                           microbiome.row.names.remove,
                                           assay)
    xenium.obj[[assay]] <- with.snvs.list[[1]]
    xenium.obj[[fov]] <- with.snvs.list[[2]]
    # generating the assay and fov with genes and microbiome probes
    assay <- "Xenium.with.microbiome"
    fov <- "fov.with.microbiome"
    with.microbiome.list <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                           xenium.list$microns,
                                           segmentations.data,
                                           snv.row.names.remove,
                                           assay)
    xenium.obj[[assay]] <- with.microbiome.list[[1]]
    xenium.obj[[fov]] <- with.microbiome.list[[2]]
    # generating the assay and fov with only genes
    assay <- "Xenium"
    fov <- "fov"
    xenium.list.sub <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                                 xenium.list$microns,
                                                 segmentations.data,
                                                 all.row.names.remove,
                                                 assay)
    xenium.obj[[assay]] <- xenium.list.sub[[1]]
    xenium.obj[[fov]] <- xenium.list.sub[[2]]
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
} else if (opt$with_microbiome) {
    # compiling the list of snvs features that will be removed
    if (is.na(opt$microbiome_probes)) {
        bacteroides_targets  <- targets[grepl("Bacteroides", targets)]
        escherichia_targets  <- targets[grepl("Escherichia", targets)]
        fusobacterium_targets  <- targets[grepl("Fusobacterium", targets)]
        hpv_targets  <- targets[grepl("papillomavirus", targets)]
        microbiome.row.names.remove <- c(bacteroides_targets, escherichia_targets, fusobacterium_targets, hpv_targets)
    } else {
        microbiome.row.names.remove <- readLines(opt$microbiome_probes)
    }
    microbiome.row.names.remove <- probe_hyphens(microbiome.row.names.remove)
    print("Removing the following Microbiome probes from the default Xenium matrix")
    print(microbiome.row.names.remove)
    # each assay and FOV will be associated with the same segmentation data since we aren't subsetting cells yet.
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    # generating the object with the initial assay with everything
    assay <- "Xenium.with.microbiome"
    fov <- "fov.with.microbiome" #If you are using Seurat v4.3.0 do not rename this as it will result in you being unable to plot molecules.
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay) #In Seurat 5 this creates a v5 assay.
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- all.coords
    
    # generating the assay and fov with only genes
    assay <- "Xenium"
    fov <- "fov"
    xenium.list.sub <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                                 xenium.list$microns,
                                                 segmentations.data,
                                                 microbiome.row.names.remove,
                                                 assay)
    xenium.obj[[assay]] <- xenium.list.sub[[1]]
    xenium.obj[[fov]] <- xenium.list.sub[[2]]
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
} else if (opt$with_snvs) {
    # compiling the list of snvs features that will be removed
    if (is.na(opt$snv_probes)) {
        wt_targets <- targets[grepl("-WT",targets)]
        alt_targets <- targets[grepl("-ALT",targets)]
        snv.row.names.remove <- c(wt_targets, alt_targets)
    } else {
        snv.row.names.remove <- readLines(opt$snv_probes)
    }
    snv.row.names.remove <- probe_hyphens(snv.row.names.remove)
    print("Removing the following SNV probes from the default Xenium matrix")
    print(snv.row.names.remove)
    # each assay and FOV will be associated with the same segmentation data since we aren't subsetting cells yet.
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    # generating the object with the initial assay with everything
    assay <- "Xenium.with.snvs"
    fov <- "fov.with.snvs"
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay) #In Seurat 5 this creates a v5 assay.
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- all.coords
    # generating the assay and fov with only genes
    assay <- "Xenium"
    fov <- "fov"
    xenium.list.sub <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                                 xenium.list$microns,
                                                 segmentations.data,
                                                 snv.row.names.remove,
                                                 assay)
    xenium.obj[[assay]] <- xenium.list.sub[[1]]
    xenium.obj[[fov]] <- xenium.list.sub[[2]]
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
} else {
    assay <- "Xenium"
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay)
    xenium.obj[["BlankCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = xenium.list$matrix[["Negative Control Probe"]])
    fov <- "fov" #If you are using Seurat v4.3.0 do not rename this as it will result in you being unable to plot molecules.
    xenium.obj[[fov]] <- all.coords
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
}
#check to make sure that plotting with transcript/cell coordinates works
pdf(paste(out_path,"VlnPlot_QC_",sampleID,".pdf", sep=''),useDingbat=FALSE, width=8, height=6)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
dev.off()
# pdf(paste("ImageDimPlot_QC_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #This shows all of the dots on one plot
# ImageDimPlot(xenium.obj, fov = "fov", molecules = c("TFF2", "KRT7", "VCAN", "CFTR"), nmols = 20000)
# dev.off()
# pdf(paste("ImageDimPlot_QC_features_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# ImageFeaturePlot(xenium.obj, features = c("TFF2", "KRT7", "VCAN", "CFTR"), max.cutoff = c(20, 20, 15, 5), size = 0.25, cols = c("lightgrey", "blue", "yellow"), border.size= NA)
# dev.off()
# pdf(paste("ImageDimPlot_QC_features_2color_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# ImageFeaturePlot(xenium.obj, features = c("TFF2", "KRT7", "VCAN", "CFTR"), max.cutoff = c(20, 20, 15, 5), size = 0.25, cols = c("blue", "yellow"), border.size= NA)
# dev.off()
# save the raw
saveRDS(xenium.obj, paste(out_path,sampleID,"_raw.rds",sep=''))
# some cells will be empty so you need to remove them
# works fine Seruat v4.3.0.1 but not in v5.0.1
# counts <- xenium.obj@assays$Xenium@counts
# empty_cells <- colnames(counts)[(colSums(counts) == 0)]
source("/diskmnt/Projects/Users/austins2/tools/subset_obj_seurat.R")
if (sum(xenium.obj$nCount_Xenium <= opt$counts_cutoff) > 0) {
    xenium.empty <- subset_opt(xenium.obj, subset = nCount_Xenium <= opt$counts_cutoff)
    empty_cells = Cells(xenium.empty)
} else {
    empty_cells = c()
}
write.table(empty_cells,paste(out_path,"Empty_cells_",sampleID,'.tsv',sep=''),sep="\t",quote=FALSE)
rm(xenium.empty)
# cells_to_keep <- colnames(counts)[(colSums(counts) > opt$counts_cutoff)]
# xenium.sub <- subset(xenium.obj, cells = cells_to_keep)
#xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay)
xenium.sub <- subset_opt(xenium.obj, subset = nCount_Xenium > opt$counts_cutoff)
xenium.sub <- subset_opt(xenium.sub, subset = nFeature_Xenium > opt$features_cutoff)


if (!is.null(cell_removal_csv) && nzchar(cell_removal_csv)){
    message("Removing input cells")
    removal_data <- read.csv(cell_removal_csv, skip = 3, header = FALSE, stringsAsFactors = FALSE)
    cells_to_remove <- trimws(removal_data[[1]])
    head(cells_to_remove)

    cells <- setdiff(Cells(xenium.sub), cells_to_remove)
    xenium.sub <- subset_opt(xenium.sub, cells = cells)
    message("Input cells were removed!")
}

if (!is.null(cell_keep_csv) && nzchar(cell_keep_csv)) {
    message("Keep input cells")
    
    keep_data <- read.csv(cell_keep_csv, skip = 3, header = FALSE, stringsAsFactors = FALSE)
    cells_to_keep <- trimws(keep_data[[1]])
    head(cells_to_keep)
    
    cells <- intersect(Cells(xenium.sub), cells_to_keep)
    
    xenium.sub <- subset_opt(xenium.sub, cells = cells)
    message("Input cells were kept!")
}


print("Normalizing")
# normalize,umap,cluster, SCTv2 handles low number of genes very well compared to LogNormalize (see Ilya's results)
xenium.sub <- SCTransform(xenium.sub, assay = "Xenium", return.only.var.genes = F, vst.flavor="v2")
xenium.sub <- RunPCA(xenium.sub, npcs = num_pcs, features = rownames(xenium.sub))
xenium.sub <- RunUMAP(xenium.sub, dims = 1:num_pcs, reduction.name = paste("umap.",num_pcs,"PC", sep = ""), reduction.key = paste("UMAP",num_pcs,"PC_",sep=""))
xenium.sub <- FindNeighbors(xenium.sub, reduction = "pca", dims = 1:num_pcs)
xenium.sub <- FindClusters(xenium.sub, resolution = 0.3)
DefaultFOV(xenium.sub, assay='SCT') <- 'fov'
# the above line results in the following warning and may cause problems down the line. to do this properly you would copy the current 'fov' slot and save it to the FOV slot for the SCT assay I think. Not 100% sure how to do this yet.
# Warning: FOV 'fov' currently associated with assay 'Xenium', changing to 'SCT'

# save the processed object.
Idents(xenium.sub) <- "seurat_clusters"
# saveRDS(xenium.sub, file='HT270P1-S1H1A1Fp1Us1_1_processed.rds')
# adding in cell labels (here they are based on the xenium browser clusters)
# can't use the traditional method for seurat cluster adding and I think it has to do with the xenium browser clusters not being a factor.
xenium_browser_clusters <- read.table(paste(input,"/analysis/clustering/gene_expression_graphclust/clusters.csv",sep=''), sep = ",", header=T)
xenium_browser_clusters <- xenium_browser_clusters[(xenium_browser_clusters$Barcode %in% Cells(xenium.sub)),]
rownames(xenium_browser_clusters) <- xenium_browser_clusters$Barcode
xenium_browser_clusters$Barcode <- NULL
xenium.sub <- AddMetaData(object = xenium.sub, metadata = xenium_browser_clusters, col.name = "xenium_browser_graph_based_clusters")
# adding in cell labels (here they are based on the xenium browser clusters)
# can't use the traditional method for seurat cluster adding and I think it has to do with the xenium browser clusters not being a factor.
# could also be due to the xenium object not containing the right cell type labels.
#cell_types <- c("your","clusters","go","here")
#xenium.sub$cell_type_individual <- NA
#xenium_clusters = sort(unique(xenium.sub$xenium_browser_graph_based_clusters))
#for (i in 1:length(xenium_clusters)) {
#    value = xenium_clusters[i]
#    cell_type = cell_types[i]
#    xenium.sub$cell_type_individual[xenium.sub$xenium_browser_graph_based_clusters == value] <- cell_type
#}
# Setting the default values to be included in plots
DefaultAssay(xenium.sub) <- "SCT"
Idents(xenium.sub) <- "seurat_clusters"
DefaultFOV(xenium.sub, assay='SCT') <- 'fov'
DefaultBoundary(xenium.sub[["fov"]]) <- "segmentation"
# save the processed object.
saveRDS(xenium.sub, paste(out_path,sampleID,"_processed.rds",sep=''))
xenium.sub$original_Xenium_barcode <- rownames(xenium.sub@meta.data)
tmp_df <- xenium.sub@meta.data[,c("original_Xenium_barcode","seurat_clusters")]
colnames(tmp_df) <- c("cell_id","seurat_clusters")
write.table(tmp_df, paste0(out_path,sampleID,"_seurat_clusters.tsv"),row.names=FALSE, sep="\t", quote=FALSE)
colnames(tmp_df) <- c("cell_id","group")
write.table(tmp_df, paste0(out_path,sampleID,"_seurat_clusters.csv"),row.names=FALSE, sep=",", quote=FALSE)
print("generating_QC_plots")
pdf(paste(out_path,"DimPlots_QC_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6)
print(rasterize(DimPlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), group.by = "seurat_clusters", label=TRUE, label.size=4, raster=TRUE), layers='Point', dpi=300))
print(rasterize(DimPlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), group.by = "xenium_browser_graph_based_clusters", label=TRUE, label.size=4, raster=FALSE), layers='Point', dpi=300))
print(rasterize(FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = "nCount_Xenium", label=TRUE, label.size=4, raster=TRUE), layers='Point', dpi=300))
print(rasterize(FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = "nFeature_Xenium", label=TRUE, label.size=4, raster=TRUE), layers='Point', dpi=300))
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.color=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "segmentation"`
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.size=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "centroids"`
print(rasterize(ImageFeaturePlot(xenium.sub, features = "nCount_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("blue", "yellow"), border.color = NA), layers='Polygon', dpi=600)) #genes with high expression can have low misleading background so we set all zero value to lightgrey and max.cutoff is 90th percentile
print(rasterize(ImageFeaturePlot(xenium.sub, features = "nFeature_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("blue", "yellow"), border.color = NA), layers='Polygon', dpi=600)) #max.cutoff is 90th percentile
print(rasterize(ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", border.color=NA), layers='Polygon', dpi=600))
#ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", border.color = "white", border.size=0.001)
print(rasterize(ImageDimPlot(xenium.sub, fov = "fov", group.by="xenium_browser_graph_based_clusters", border.color=NA), layers='Polygon', dpi=600))
#ImageDimPlot(xenium.sub, fov = "fov", group.by="xenium_browser_graph_based_clusters", border.color = "white", border.size=0.001)
dev.off()
if (is.na(opt$plot_features)){
    print("skipping plotting of features as no features were provided")
} else {
    genes_filtered <-c()
    plot_features=str_split(opt$plot_features, pattern=',')[[1]]
    plot_features=unique(plot_features)
    sample_genes <- rownames(x = xenium.sub@assays$SCT@counts)
    
    for (gene in plot_features) {
        if (gene %in% sample_genes) {
            genes_filtered <- c(genes_filtered,gene)
        }
    }
    print("the following genes are unique and are present in the SCT assay for inclusion in plots")
    print(genes_filtered)
    pdf(paste(out_path,"FeaturePlots_markers_of_interest_SCT_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6)
    for (marker in genes_filtered) {
        print(rasterize(FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = marker, label=TRUE, label.size=4, raster=TRUE), layers='Point', dpi=600))
    }
    dev.off()
    pdf(paste(out_path,"ImageDimPlot_QC_features_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.color=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "segmentation"`
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.size=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "centroids"`
    for (marker in genes_filtered){
        print(rasterize(ImageFeaturePlot(xenium.sub, features = marker, max.cutoff = 'q95', size = 0.15, cols = c("lightgrey", "darkblue", "yellow"), border.color= NA), layers='Polygon', dpi=600)) #genes with high expression can have low misleading background so we set all zero value to lightgrey and max.cutoff is 90th percentile
        print(rasterize(ImageFeaturePlot(xenium.sub, features = marker, max.cutoff = 'q90', size = 0.15, cols = c("darkblue", "yellow"), border.color= NA), layers='Polygon', dpi=600)) #max.cutoff is 90th percentile
    }
    dev.off()
}

# Further examples of how to crop to new FOVs and 
# 
# options(future.globals.maxSize= +Inf)
# new_markers = c("CD68","CD4","FGG","KRT7","KRT19","ACTA2","PECAM1","CD8A")
# new_colors = c("magenta","cyan","yellow","pink","orange","red","green","white")
# Idents(nano.sub) <- "orig.ident"
# DefaultBoundary(xenium.sub[["fov"]]) <- "segmentation"
# pdf('test7_6.pdf',useDingbats = F, width = 8, height = 6)
# #ImageDimPlot(nano.obj, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", #size = 0.5, 
#              border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", #size = 0.5, 
#              #border.size = NA,
#              border.color = "white", border.size=0.001,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "fov", #axes = TRUE, 
#              #border.size = NA, 
#              border.color = "white", border.size = 0.001, cols = NA,
#              mols.size = 0.1,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes = F, stroke=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()

# #Zoom1
# # c(0, 3000), y = c(0, 3000)
# cropped.coords <- Crop(xenium.sub[["fov"]], x = c(0, 3000), y = c(0, 3000), coords = "plot")
# xenium.sub[["zoom1"]] <- cropped.coords
# DefaultBoundary(xenium.sub[["zoom1"]]) <- "segmentation"
# source("/diskmnt/Projects/Users/austins2/tools/subset_obj_seurat.R")
# pdf('test7_zoom1.pdf',useDingbats = F, width = 8, height = 6)
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "zoom1", group.by="seurat_clusters", #size = 0.5,
#              border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "zoom1", group.by="seurat_clusters", #size = 0.5,
#              #border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "zoom1", #axes = TRUE,
#              #border.size = NA,
#              #border.color = "white", border.size = 0.2, cols = "black",
#              mols.size = 0.5,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# xenium.sub = subset_opt(xenium.sub, subset = seurat_clusters %in% c(0,1,2))
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "zoom1", #axes = TRUE,
#              #border.size = NA,
#              #border.color = "white", border.size = 0.2, cols = "black",
#              mols.size = 0.5,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()

# # Zoom2 
# # x = c(1000, 4000), y = c(1000, 4000)
# cropped.coords <- Crop(xenium.sub[["fov"]], x = c(1000, 4000), y = c(1000, 4000), coords = "plot")
# xenium.sub[["zoom2"]] <- cropped.coords
# DefaultBoundary(xenium.sub[["zoom2"]]) <- "segmentation"
# pdf('test7_zoom2.pdf',useDingbats = F, width = 8, height = 6)
# #ImageDimPlot(nano.obj, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "zoom2", group.by="seurat_clusters", #size = 0.5, 
#              border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "zoom2", group.by="seurat_clusters", #size = 0.5, 
#              #border.size = NA,
#              border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot.FIX(xenium.sub, fov = "zoom2", axes = TRUE,
#              #border.size = NA, 
#              border.color = "white", border.size = 0.2, cols = "black",
#              mols.size = 0.3,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes=F, stroke=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()
