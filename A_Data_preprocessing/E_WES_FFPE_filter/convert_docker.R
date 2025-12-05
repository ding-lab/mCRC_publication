#!/bin/bash

# set arguments
rm(list = ls())
args <- commandArgs(trailingOnly = T)
wd <- args[1]
maf <- args[2]
sample_input <- args[3]

## set up library and directory
setwd(wd)
library(openxlsx)
library(MicroSEC)
library(GenomeInfoDb)
library(dplyr)


## mutation file
df_CCLE <- utils::read.csv(maf,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = "\t",
    comment.char = "#"
)

df_CCLE <- df_CCLE %>% dplyr::filter(Tumor_Sample_Barcode == sample_input)
df_CCLE$Sample <- df_CCLE$Tumor_Sample_Barcode
df_CCLE$Chr <- df_CCLE$Chromosome
df_CCLE$Chr <- as.character(df_CCLE$Chr)
df_CCLE$Reference <- df_CCLE$Reference_Allele
df_CCLE$Tumor_Seq <- df_CCLE$Tumor_Seq_Allele2
seqlevelsStyle(df_CCLE$Chr) <- "UCSC"
df_CCLE <- df_CCLE %>% dplyr::filter(Chr != "M")
write.xlsx(df_CCLE, paste0(sample_input, "_Mutations.xlsx"))

## convert funcition
new <- function(mutation_file, organism) {
    Hugo_Symbol <- NULL
    Protein_Change <- NULL
    Start_position <- NULL
    End_position <- NULL
    Variant_Type <- NULL
    Reference <- NULL
    Tumor_Seq <- NULL
    Mut_type <- NULL
    Chr <- NULL
    Start <- NULL
    End <- NULL
    Ref <- NULL
    Alt <- NULL
    Alt_length_1 <- NULL
    Alt_length_2 <- NULL
    PRE_ins <- NULL
    PRE_del <- NULL
    Alt_ins <- NULL
    Alt_del <- NULL
    Alt_snv <- NULL
    Ref_ins <- NULL
    Ref_del <- NULL
    Ref_snv <- NULL
    Neighbor_start_1 <- NULL
    Neighbor_end_1 <- NULL
    Neighbor_start_2 <- NULL
    Neighbor_end_2 <- NULL
    Pre_Neighbor <- NULL
    Alt_indel <- NULL
    POST_ins <- NULL
    Post_Neighbor <- NULL
    Alt_length <- NULL
    Ref_indel <- NULL
    Pos <- NULL

    df_mutation <- read.xlsx(mutation_file, sheet = 1)
    fun_load_genome(organism)

    fun_genome <- function(x, y) {
        r <- NULL
        for (i in seq_len(length(x))) {
            r <- c(r, as.character(ref_genome[[x[i]]][y[i]]))
        }
        return(r)
    }

    fun_genome_2 <- function(x, y, z) {
        r <- NULL
        for (i in seq_len(length(x))) {
            r <- c(r, as.character(ref_genome[[x[i]]][y[i]:z[i]]))
        }
        return(r)
    }

    df_mutation <- df_mutation %>% mutate(
        Gene = Hugo_Symbol,
        HGVS.p = Protein_Change, Start = Start_position,
        End = End_position, Mut_type = Variant_Type, Ref = Reference,
        Alt = Tumor_Seq
    )
    df_mutation$Pos <- df_mutation$Start_Position

    df_mutation <- df_mutation %>% dplyr::mutate(Mut_type = tolower(Mut_type))
    df_mutation$Chr <- as.character(df_mutation$Chr)
    seqlevelsStyle(df_mutation$Chr) <- "UCSC"
    df_mutation$Pos <- as.integer(df_mutation$Pos)
    df_mutation$Start <- as.integer(df_mutation$Start)
    df_mutation$End <- as.integer(df_mutation$End)
    df_mutation <- df_mutation %>% dplyr::mutate(PRE_del = fun_genome(
        Chr,
        Start - 1
    ), PRE_ins = fun_genome(Chr, Start), POST_ins = fun_genome(
        Chr,
        End
    ), Alt_length_1 = nchar(Ref), Alt_length_2 = nchar(Alt))
    df_mutation <- df_mutation %>% dplyr::mutate(Mut_type = ifelse(Mut_type ==
        "snp", "snv", Mut_type))
    df_mutation <- df_mutation %>% dplyr::mutate(Mut_type = ifelse(Mut_type ==
        "dnp", "snv", Mut_type))
    df_mutation <- df_mutation %>% dplyr::mutate(Mut_type = ifelse(Mut_type ==
        "tnp", "snv", Mut_type))
    df_mutation <- df_mutation %>% dplyr::mutate(Mut_type = ifelse(Mut_type ==
        "onp", "snv", Mut_type))
    df_mutation <- df_mutation %>% dplyr::mutate(Alt_length = (((Alt_length_1 -
        Alt_length_2) + abs(Alt_length_1 - Alt_length_2)) / 2) +
        Alt_length_2, Ref_ins = ifelse(Mut_type == "ins", PRE_ins,
        ""
    ), Ref_del = ifelse(Mut_type == "del", paste(PRE_del,
        Ref,
        sep = ""
    ), ""), Ref_snv = ifelse(Mut_type == "snv",
        Ref, ""
    ), Alt_ins = ifelse(Mut_type == "ins", paste(PRE_ins,
        Alt,
        sep = ""
    ), ""), Alt_del = ifelse(Mut_type == "del",
        PRE_del, ""
    ), Alt_snv = ifelse(Mut_type == "snv", Alt,
        ""
    ))
    df_mutation <- df_mutation %>% dplyr::mutate(Alt_indel = paste(Alt_ins,
        Alt_del, Alt_snv,
        sep = ""
    ), Ref_indel = paste(Ref_ins,
        Ref_del, Ref_snv,
        sep = ""
    ))
    df_mutation <- df_mutation %>% dplyr::mutate(Neighbor_start_1 = Start -
        20, Neighbor_end_1 = Start - 1, Neighbor_start_2 = End +
        1, Neighbor_end_2 = End + 20)
    df_mutation <- df_mutation %>% dplyr::mutate(Pre_Neighbor = fun_genome_2(
        Chr,
        Neighbor_start_1, Neighbor_end_1
    ), Post_Neighbor = fun_genome_2(
        Chr,
        Neighbor_start_2, Neighbor_end_2
    ))
    df_mutation <- df_mutation %>% dplyr::mutate(Neighborhood_sequence = ifelse(Mut_type ==
        "ins", paste(Pre_Neighbor, Alt_indel, POST_ins, stringr::str_sub(
        Post_Neighbor,
        1, 19
    ), sep = ""), ifelse(Mut_type == "del", paste(Pre_Neighbor,
        Post_Neighbor,
        sep = ""
    ), paste(Pre_Neighbor, Alt, Post_Neighbor,
        sep = ""
    ))))
    df_mutation <- df_mutation %>% dplyr::mutate(
        Mut_type = paste(Alt_length,
            "-", Mut_type,
            sep = ""
        ), Ref = Ref_indel, Alt = Alt_indel,
        Pos = ifelse(stringr::str_detect(Mut_type, pattern = "del"), Pos -
            1, Pos)
    )
    df_mutation <- df_mutation %>% dplyr::select(
        -PRE_del, -PRE_ins,
        -POST_ins, -Alt_length_1, -Alt_length_2, -Alt_length,
        -Ref_ins, -Ref_del, -Ref_snv, -Alt_ins, -Alt_del, -Alt_snv,
        -Alt_indel, -Ref_indel, -Neighbor_start_1, -Neighbor_start_2,
        -Neighbor_end_1, -Neighbor_end_2, -Pre_Neighbor, -Post_Neighbor,
        -Hugo_Symbol, -Start_position, -End_position, -Variant_Type,
        -Reference, -Tumor_Seq, -Protein_Change, -Start, -End
    )
    return(df_mutation)
}

## run
df_mutation <- new(paste0(sample_input, "_Mutations.xlsx"), "hg38")
write.xlsx(df_mutation, paste0(sample_input, "_Mutations_converted.xlsx"))
