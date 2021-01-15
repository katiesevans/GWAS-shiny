library(fst)
library(tidyverse)


# make one condition file for all traits for 96 strain gwas rerun by gaotian 12/2020

load("~/Dropbox/AndersenLab/RCode/GWAS/Data/GWAS_96strainset/2170602_gwas6-8_easysorter_assayControl_regressed.Rda")
unique(allg_regressed$condition)

wd <- "~/Dropbox/AndersenLab/LabFolders/Katie/git/GWAS-shiny/data/"
for(d in unique(allg_regressed$condition)) {
    files <- grep(d, list.files(glue::glue("{wd}/mappings/")), value = T)
    dmap <- NULL
    for(f in files) {
        map <- readr::read_tsv(glue::glue("{wd}/mappings/{f}"))
        dmap <- rbind(dmap, map)
    }
    fst::write_fst(dmap, glue::glue("{wd}/{d}_pmd.fst"))
}


for(d in unique(allg_regressed$condition)) {
    map <- fst::read_fst(glue::glue("{wd}/mappings/96strain/{d}_pmd.fst")) %>%
        dplyr::select(marker:aboveBF, peak_id) %>%
        dplyr::left_join(qtlpeaks) %>%
        dplyr::distinct() %>%
        dplyr::mutate(condition = d,
                      trait = stringr::str_split_fixed(trait, "_", 2)[,1]) %>%
        dplyr::select(marker:log10p, condition, everything())
    
    
    fst::write_fst(map, glue::glue("{wd}/mappings/96strain/{d}_pmd.fst"))
    
}

paste(unique(allg_regressed$condition), collapse = "','")

t <- allg_regressed %>%
    dplyr::filter(!grepl("red|green|yellow", trait)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(trait) %>%
    dplyr::pull(trait)

paste(t, collapse = "','")

# trait data
pheno <- allg_regressed %>%
    dplyr::select(condition, strain, trait, phenotype)
fst::write_fst(pheno, "~/Dropbox/AndersenLab/LabFolders/Katie/git/GWAS-shiny/data/96strain_pheno.fst")


# add genotype information
geno_matrix <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/git/GWAS-shiny/data/Genotype_Matrix.tsv") %>%
    tidyr::gather(strain, allele, -(CHROM:ALT))


for(d in unique(allg_regressed$condition)) {
    map <- fst::read_fst(glue::glue("{wd}/mappings/96strain/{d}_pmd.fst")) %>%
        dplyr::mutate(EIGEN = -log10(0.05/489.044),
                      aboveEIGEN = log10p > EIGEN) %>%
        dplyr::select(marker:aboveBF, EIGEN, aboveEIGEN, everything())
    
    
    fst::write_fst(map, glue::glue("{wd}/mappings/96strain/{d}_pmd.fst"))
    
}
