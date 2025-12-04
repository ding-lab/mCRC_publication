#!/bin/bash

Rscript seurat_5.0.1_xenium_v6.2.R \
-i Xenium_OutPut_path \
-s Sample_ID \
-o OutPut_dir \
--plot_features 'KRT7,TFF2,CFTR,EPCAM,MKI67,TOP2A,PECAM1,VWF,CD3E,ACTA2,C7,PDGFRA,PDGFRA,SPP1,SLC2A1,LAMC2,CD68,CD163' \
--with_snvs \
--snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_probe_names.tsv
