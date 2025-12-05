# ffpefilter version 0.0.1 #
Use MicroSec for FFPE filtering. ffpefilter pipeline is a pipeline designed to do variant filtering for FFPE samples by using microsec tool. It uses the docker image scao/microsec:0.0.1 for efficiently launching jobs on compute1 cluster (LSF job system)


## Install ##

The only thing you need to do is to do "git clone https://github.com/ding-lab/ffpefiltering.git" 

## Usage ##

1. Create folder (rdir or use anyname you like) for running multiple samples 

2. Create subdirectory under rder for each sample that need to run ffpe filtering

For each sample direcotry, it should contain a maf file and a tumor bam file and associated index file for each sample named as:

samplename.formated.maf

samplename.T.bam

samplename.T.bam.bai

The maf file should contain the following columns, which can be extracted from the original maf file from somaticwrapper pipeline (remember to remove _T for the sample name in the Tumor_Sample_Barcode)

Hugo_Symbol     Chromosome      Start_Position  End_Position    Strand  Variant_Classification  Variant_Type    Reference_Allele        Tumor_Seq_Allele1       Tumor_Seq_Allele2       Tumor_Sample_Barcode

SAMD11  chr1    925818  925818  +       Intron  SNP     G       G       T       ALCH-ABBG-TTP1-A

3. Enter the folder which contains ffpefilter.pl script, type
 
perl ffpefilter.pl  --rdir --log --groupname --users --step 

rdir = full path of the folder holding files for this run (user must provide)

log = full path of the folder for saving log file; usually upper folder of rdir

groupname = job group name

users = user name for job group

step = run this pipeline step by step. (user must provide)

## Examples ##

See work_log_alchemist_test for how to submit the jobs.

## Contact ##

Yingduo Song and Song Cao
