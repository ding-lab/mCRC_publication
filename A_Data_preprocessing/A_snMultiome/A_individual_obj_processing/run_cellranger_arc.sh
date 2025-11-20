#!/bin/bash

eval "$(conda shell.bash hook)"

export TMPDIR=/diskmnt/Projects/PECGS_analysis/Sn_Multiome/tmp/
export PATH=$PATH:/diskmnt/Projects/Users/rliu/Software/cellranger-arc-2.0.0/

sample=$1
projectfolder=/diskmnt/Projects/MetNet_analysis_2/Colorectal/epeng/sn_multiome/
workingfolder=${projectfolder}cellranger-arc/preprocessing/

echo "Processing sample: ${sample}"

# Create directories if they don't exist
mkdir -p "${projectfolder}/Scrublet/combo_merged/${sample}"
mkdir -p "${projectfolder}/rds_objects/withSoupX_emptyDrops/"

cd ${projectfolder}cellranger-arc/

if [ -d "$sample"/outs ]
then
	echo "Cellranger for $sample exists, moving to the next step"
else
	/diskmnt/Projects/Users/rliu/Software/cellranger-arc-2.0.0/cellranger-arc count \
	 --id ${sample} \
	 --reference /diskmnt/Datasets/Reference/Cellranger-ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
	 --libraries ${workingfolder}${sample}/library.csv \
	 --disable-ui --localcores 40 --localmem 100  
fi

if [ -f ${projectfolder}/Scrublet/combo_merged/${sample}/${sample}_combo_scrublet_output_table.csv ]
	then
		echo "Scrublet for $sample exists, moving to next step"
	else
		cd ${projectfolder}/Scrublet/RNA
		conda activate scrublet&&bash /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/Scrublet/scripts/scrublet-RNA-one-sample.sh ${projectfolder}/cellranger-arc/ $sample

		cd ${projectfolder}/Scrublet/ATAC
		bash /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/Scrublet/scripts/scrublet-ATAC-one-sample.sh ${projectfolder}/cellranger-arc/ $sample

		cd ${projectfolder}/Scrublet/combo_merged
		bash /diskmnt/Projects/Users/austins2/tools/scrublet-auto-combining.sh ${projectfolder}/cellranger-arc/ ${projectfolder}/Scrublet/RNA/ ${projectfolder}/Scrublet/ATAC/
		
		conda deactivate
fi

cd ${projectfolder}/rds_objects/
conda activate signac1.8&&Rscript /diskmnt/Projects/Users/allakarpova/scripts/multiome_sc_analysis/seurat_pipeline_v5.2.combo.R -s ${sample} \
 -d ${projectfolder}/cellranger-arc/${sample} \
 -m /diskmnt/Projects/Users/allakarpova/Tools/miniconda3/envs/signac1.8/bin/macs2 \
 -o ${projectfolder}/rds_objects/withSoupX_emptyDrops/ \
 --chrom_size /diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt \
 --prf_min 1000 --pct_min 15 --ns_max 5 --pc_first 2 --pc_num 50 \
 --scrublet ${projectfolder}/Scrublet/combo_merged/${sample}
conda deactivate

