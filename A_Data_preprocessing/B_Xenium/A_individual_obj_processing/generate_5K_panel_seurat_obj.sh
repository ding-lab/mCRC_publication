#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate Seurat5.0.1

Rscript seurat_5.0.1_xenium_v6.1.R \
-i Xenium_OutPut_Path \
-s SampleID \
-o OutPut_Dir \
--plot_features 'EPCAM,MKI67,LGR5,PECAM1,CD4,CD8A,ACTA2,PDGFRA,CD68' \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_panel_5k_v1_DHYHMW_snv_probe_names.tsv \
--with_microbiome \
--microbiome_probes /diskmnt/Projects/Users/austins2/tools/Xenium_human_panel_5k_v1_DHYHMW_microbiome_probe_names.tsv \
--counts_cutoff 20 \
--features_cutoff 10 

conda deactivate