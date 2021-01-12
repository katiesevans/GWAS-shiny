# library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(linkagemapping)
library(shiny)
library(DT)
library(plotly)
library(fst)
library(curl)
library('R.utils')
library(RSQLite)

# load RIAIL regressed phenotype data to get possible condition/traits
data("eQTLpeaks")
data("probe_info")
allRIAILsregressed <- fst::read_fst("data/allRIAILsregressed.fst")

# pxgplot with fake cross object
pxgplot_fake <- function (cross, map, parent = "N2xCB4856", tsize = 20) {
    peaks <- map %>% 
        na.omit() %>%
        dplyr::group_by(iteration) %>% 
        dplyr::filter(!is.na(var_exp)) %>% 
        dplyr::do(head(., n = 1))
    if (nrow(peaks) == 0) {
        stop("No QTL identified")
    }
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    pheno <- cross$pheno %>% dplyr::select_(map$trait[1])
    
    # get geno
    uniquechr <- unique(stringr::str_split_fixed(uniquemarkers, "_", 2)[,1])
    
    if(length(uniquechr) > 1) {
        geno <- data.frame(cross$geno[[uniquechr[1]]])
        for(i in 2:length(uniquechr)) {
            geno <- cbind(geno, cross$geno[[uniquechr[i]]])
        }
    } else {
        geno <- data.frame(cross$geno[[uniquechr]])
    }
    
    geno <- geno %>%
        dplyr::select(which(colnames(.) %in% uniquemarkers)) %>% 
        data.frame(., pheno)
    colnames(geno)[1:(ncol(geno) - 1)] <- sapply(colnames(geno)[1:(ncol(geno) -  1)],
                                                 function(marker) {paste(unlist(peaks[peaks$marker == gsub("\\.", "-", marker), c("chr", "pos")]), collapse = ":") })
    colnames(geno)[ncol(geno)] <- "pheno"
    split <- tidyr::gather(geno, marker, genotype, -pheno)
    split$genotype <- sapply(split$genotype, function(x) {
        if (is.na(x)) {
            return(NA)
        }
        if (parent == "N2xCB4856") {
            if (x == 1) {
                "N2"
            }
            else {
                "CB4856"
            }
        }
    })
    split$genotype <- factor(split$genotype, levels = c("N2", "CB4856", "LSJ2", "AF16", "HK104"))
    split <- split %>% 
        tidyr::drop_na(genotype) %>% 
        dplyr::mutate(chr = (as.character(stringr::str_split_fixed(marker, ":", 2)[, 1])), 
                      pos = as.numeric(as.character(stringr::str_split_fixed(marker, ":", 2)[, 2]))) %>% 
        dplyr::arrange(chr, pos)
    split$marker <- factor(split$marker, levels = unique(split$marker))
    ggplot2::ggplot(split) + ggplot2::scale_fill_manual(values = c(N2 = "orange", CB4856 = "blue", LSJ2 = "green", AF16 = "indianred", HK104 = "gold")) + 
        ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), alpha = 0.8, size = 0.5, width = 0.1) + 
        ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.shape = NA, alpha = 0.8) + 
        ggplot2::facet_wrap(~marker, ncol = 5) + ggplot2::theme_bw(tsize) + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", color = "black"), 
                       axis.text.y = ggplot2::element_text(face = "bold", color = "black"), 
                       axis.title.x = ggplot2::element_text(face = "bold", color = "black", vjust = -0.3), 
                       axis.title.y = ggplot2::element_text(face = "bold", color = "black"), 
                       strip.text.x = ggplot2::element_text(face = "bold", color = "black"), 
                       strip.text.y = ggplot2::element_text(face = "bold", color = "black"), 
                       plot.title = ggplot2::element_text(face = "bold", vjust = 1), 
                       legend.position = "none", 
                       panel.background = ggplot2::element_rect(color = "black", size = 1.2)) + 
        ggplot2::ggtitle(peaks$trait[1]) + 
        ggplot2::labs(x = "Genotype", y = "Phenotype")
}


# Define UI for application
ui <- fluidPage(
   
   # Application title
   titlePanel("Linkage Mapping Analysis - Andersen Lab"),
   
   # no sidebar layout
   shiny::wellPanel(
       # input well
       shiny::fluidRow(
           # conditions
           column(3, shiny::selectInput(inputId = "drug_input", label = "Condition:", choices = sort(unique(allRIAILsregressed$condition)), selected = NULL)),
           # set
           column(2, shiny::radioButtons(inputId = "set_input", label = "RIAIL set:", choices = c(1, 2), selected = 2, inline = TRUE)),
           # trait
           column(3, shiny::selectInput(inputId = "trait_input", label = "Trait:", choices = c("", unique(allRIAILsregressed$trait)), selected = NULL)),
           # qtl marker
           column(3, shiny::uiOutput("qtl"))
       ),
       # code download button
       shiny::downloadButton(outputId = "code_download", label = "Get code")
   ),
      
  # Main panel for plots
  mainPanel(width = 12,
            
     shiny::tabsetPanel(type = "tabs", id = "test",
                        # shiny::tabPanel("QTL Analysis: Condition", shiny::plotOutput("genplot", height = "800px")),
                        shiny::tabPanel("QTL Analysis: Condition", 
                                        shiny::uiOutput("cond_plot")),
                        shiny::tabPanel("QTL Analysis: Trait",   
                                        shiny::uiOutput("allplots")),
                        shiny::tabPanel("eQTL Overlap",
                                        plotly::plotlyOutput("eqtlplot", height = "500px"),
                                        DT::dataTableOutput("eqtl_data")),
                        shiny::tabPanel("Candidate Genes",
                                        # shiny::uiOutput("qtl"),
                                        shiny::uiOutput("candidate_genes")),
                        shiny::tabPanel("Help",
                                        shiny::uiOutput("help_md")))
  )
)


#########################################################################################
###########                          SERVER                                   ###########  
#########################################################################################

server <- function(input, output) {
    
    # show/hide tabs based on input information
    shiny::observeEvent(input$trait_input, {
        if(input$trait_input == ""){ 
            shiny::hideTab("test",target = "QTL Analysis: Trait") } else {
            shiny::showTab("test", target = "QTL Analysis: Trait", select = TRUE)
        }
    })
    
    shiny::observeEvent(input$whichqtl, {
        if(input$whichqtl == ""){ 
            shiny::hideTab("test",target = "eQTL Overlap")
            shiny::hideTab("test",target = "Candidate Genes")} else {
                shiny::showTab("test",target = "eQTL Overlap", select = TRUE)
                shiny::showTab("test",target = "Candidate Genes")
        }
    })
    
    #############################
    ####    LOAD DATA        ####
    #############################
    loaddata <- shiny::reactive({

        cond <- input$drug_input
        
        # load data
        # annotatedmap <- fst::read_fst(glue::glue("data/drug_data/{cond}-GWER.chromosomal.annotated.fst"))
        
        # I guess fst can't read directly from internet... this is my workaround for now.
        tmp_file <- tempfile()
        fst_url <- glue::glue("https://raw.githubusercontent.com/katiesevans/LM-shiny/main/drug_data/{cond}-GWER.chromosomal.annotated.fst")
        curl::curl_download(fst_url, tmp_file, mode="wb")
        
        # files
        annotatedmap <- fst::read_fst(tmp_file)
        peaks <- annotatedmap %>%
            na.omit()
        list(annotatedmap, peaks)
        
    })
    
    #############################
    ####    ALL LOD PLOT     ####
    #############################
    # plot alllodplots for each condition (or multiple conditions)
    output$cond_plot <- shiny::renderUI({
        showModal(modalDialog(footer=NULL, size = "l", tags$div(style = "text-align: center;", "Loading...")))
        
        # get drug
        cond <- input$drug_input
        strainset <- input$set_input
        
        # load data
        peaks <- loaddata()[[2]]
        
        # plot all traits
        newmap <- peaks %>%
            dplyr::filter(set == strainset)%>%
            arrange(desc(trait)) %>%
            dplyr::mutate(n2res = ifelse(eff_size < 0, "yes", "no"),
                          ci_l_pos = as.numeric(ci_l_pos),
                          ci_r_pos = as.numeric(ci_r_pos),
                          pos = as.numeric(pos),
                          lod = as.numeric(lod),
                          condition = cond) %>%
            dplyr::mutate(trait = stringr::str_split_fixed(trait, paste0(cond, "."), 2)[,2])
        
        newmap$trait <- factor(newmap$trait, levels = unique(newmap$trait))
        
        
        # get chromosome lengths
        chr_lens <- data.frame(chr = c("I", "II", "III", "IV", "V", "X"), 
                               start = rep(1,6), 
                               end = c(14972282,15173999,13829314,17450860,20914693,17748731),
                               condition = newmap$condition[1],
                               trait = newmap$trait[1])
        
        removeModal()
        
        ########### PLOT ############
        output$genplot <- shiny::renderPlot({
            
            # plot
            ggplot(newmap)+
                aes(x=pos/1E6, y=trait)+
                theme_bw() +
                viridis::scale_fill_viridis(name = "LOD") + 
                viridis::scale_color_viridis(name = "LOD") +
                geom_segment(aes(x = ci_l_pos/1e6, y = trait, xend = ci_r_pos/1e6, yend = trait, color = lod), size = 2, alpha = 1) +
                geom_segment(data=chr_lens, aes(x =  start/1e6, xend = end/1e6, y = trait, yend=trait), color='transparent', size =0.1) +
                geom_point(aes(fill = lod, shape = n2res), color = "black",size = 3, alpha = 1)+
                scale_shape_manual(values = c("yes" = 24, "no" = 25)) +
                xlab("Genomic position (Mb)") + ylab("") +
                guides(shape = FALSE) +
                theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                      axis.ticks.y = element_blank(),
                      legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10),
                      legend.key.size = unit(.75, "cm"),
                      panel.grid.major.x = element_line(),
                      panel.grid.major.y = element_line(),
                      panel.grid.minor.y = element_blank(),
                      axis.text.y = element_text(size = 10, face = "bold", color = "black"),
                      axis.title.x = element_text(size=12, face="bold", color= "black"),
                      axis.title.y = element_blank(),
                      strip.text.x = element_text(size=12, face="bold", color="black"),
                      strip.text.y = element_text(size=12, face="bold", color="black", angle = 0),
                      strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
                      plot.title = element_text(size=12, face="bold")) +
                facet_grid(condition ~ chr, scales = "free_x", space = "free")
        })
        
        ############## PEAKS DF ###############
        output$peaksdf <- DT::renderDataTable({
            newmap %>%
                dplyr::select(chr, pos, trait, lod, var_exp, eff_size, ci_l_pos, ci_r_pos) %>%
                dplyr::arrange(desc(trait)) %>%
                dplyr::mutate(lod = round(lod, digits = 4),
                              var_exp = round(var_exp, digits = 4),
                              eff_size = round(eff_size, digits = 4))
        })
        
        removeModal()
        
        # what to show
        tagList(
            shiny::plotOutput("genplot", height = "800px"),
            h3("QTL peaks:"),
            DT::dataTableOutput("peaksdf")
        )
        
        
    })

    #############################
    ####    LOD/PXG PLOT     ####
    #############################
    # plot linkage defaults (LOD, pxg, riail pheno)
    output$allplots <- shiny::renderUI({
        
        showModal(modalDialog(footer=NULL, size = "l", tags$div(style = "text-align: center;", "Loading...")))
        
        # define variables
        trt <- input$trait_input
        cond <- input$drug_input
        strainset <- input$set_input

        # load data
        annotatedmap <- loaddata()[[1]]

        # filter data
        traitmap <- annotatedmap %>%
            dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                          set == strainset)
        
        pheno <- allRIAILsregressed %>%
            dplyr::filter(condition == cond,
                          trait == trt)
        
        # drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, pheno, strainset)
        load("data/newcross.Rda")
        drugcross <- linkagemapping::mergepheno(newcross, pheno, strainset)
        rm(newcross)

        #########
        # Plots #
        #########
        
        ################### BAR PLOT RIAIL ###################
        output$barplot <- plotly::renderPlotly({
            # riail pheno plot
            riailset <- drugcross$pheno %>%
                dplyr::filter(set == strainset)
            phenodf <- pheno %>%
                dplyr::filter(strain %in% c(riailset$strain, "N2", "CB4856")) %>%
                dplyr::mutate(strain_fill = dplyr::case_when(strain == "N2" ~ "n2",
                                                             strain == "CB4856" ~ "cb",
                                                             TRUE ~ "RIL")) %>%
                dplyr::group_by(strain, trait) %>%
                dplyr::mutate(avg_phen = mean(phenotype, na.rm = T)) %>%
                dplyr::distinct(strain, trait, .keep_all = T) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(norm_pheno = ((avg_phen - min(avg_phen, na.rm = T)) / (max(avg_phen, na.rm = T) - min(avg_phen, na.rm = T)))) %>%
                dplyr::arrange(norm_pheno)
            
            phenodf$strain <- factor(phenodf$strain, levels = unique(phenodf$strain))
            
            # plot
            barplot <- phenodf %>%
                ggplot2::ggplot(.) +
                ggplot2::aes(x = strain, y = norm_pheno, fill = strain_fill, color = strain_fill) +
                ggplot2::geom_bar(stat = "identity") +
                ggplot2::scale_fill_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
                ggplot2::scale_color_manual(values = c("n2" = "orange", "cb" = "blue", "RIL" = "grey")) +
                theme_bw(15) +
                theme(axis.ticks.x = element_blank(),
                      axis.text.x = element_blank(),
                      panel.grid = element_blank(),
                      axis.title = element_text(face = "bold", color = "black"),
                      axis.text.y = element_text(face = "bold", color = "black"),
                      legend.position = "none") +
                labs(x = "Strain", y = "Normalized phenotype")
            
            # plotly
            plotly::ggplotly(barplot + aes(text = glue::glue("Strain: {strain}")), tooltip = "text")
        })
        
        ################# LOD & PXG PLOTS ###################
        output$lodplot <- shiny::renderPlot({
            lodplot <- linkagemapping::maxlodplot(traitmap) +
                theme(plot.title = element_blank())
            
            if(nrow(traitmap %>% na.omit()) > 0) {
                pxgplot <- pxgplot_fake(drugcross, traitmap) +
                    theme(plot.title = element_blank(),
                          panel.grid = element_blank())
            } else {
                pxgplot <- ggplot2::ggplot(traitmap) +
                    geom_blank() +
                    geom_text(x = 0.5, y = 0.5, label = "No QTL")
            }
            
            cowplot::plot_grid(lodplot, pxgplot, nrow = 2, align = "v", axis = "lr")
        })
        
        
        ################# PEAKS DF ################# 
        output$peaks <- DT::renderDataTable({
            traitmap %>%
                na.omit() %>%
                dplyr::arrange(chr, pos) %>%
                dplyr::select(marker, lod, var_exp, eff_size, ci_l_pos, ci_r_pos) %>%
                dplyr::mutate(lod = round(lod, digits = 4),
                              var_exp = round(var_exp, digits = 4),
                              eff_size = round(eff_size, digits = 4))
        })
        
        removeModal()
        
        # what to show
        tagList(
            h3("QTL plots:"),
            plotly::plotlyOutput("barplot", height = "300px"),
            shiny::plotOutput("lodplot", height = "600px"),
            h3("QTL peaks:"),
            DT::dataTableOutput("peaks")
        )
        
    })
    
    #############################
    ####    GET EQTL DATA    ####
    #############################
    get_eQTL <- shiny::reactive({
        
        # first, get QTL
        qtl_marker <- input$whichqtl
        
        if(qtl_marker != "") {
            # define variables
            trt <- input$trait_input
            cond <- input$drug_input
            strainset <- input$set_input
            
            # load and filter data
            peaks <- loaddata()[[2]] %>%
                dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                              set == strainset) %>%
                na.omit() %>%
                dplyr::filter(marker == qtl_marker)
            
            # made for multiple peaks but I think just one at a time is good.
            all_eQTL <- NULL
            for(i in 1:nrow(peaks)) {
                test <- eQTLpeaks %>%
                    dplyr::filter(chr == peaks$chr[i],
                                  ci_l_pos < peaks$ci_r_pos[i],
                                  ci_r_pos > peaks$ci_l_pos[i]) %>%
                    dplyr::mutate(QTL = peaks$pos[i])
                all_eQTL <- rbind(all_eQTL, test)
            }
            
            # combine to get gene names if available
            newprobes <- probe_info %>% 
                dplyr::select(probe:wbgene) %>% 
                dplyr::distinct() %>%
                dplyr::group_by(probe) %>%
                dplyr::mutate(all_genes = paste(unique(gene), collapse = ", "),
                              all_genes2 = paste(unique(wbgene), collapse = ", "),
                              all_genes3 = paste(all_genes, all_genes2, sep = ", ")) %>%
                dplyr::select(probe, gene = all_genes3) %>%
                dplyr::mutate(gene = gsub(", NA", "", gene)) %>%
                dplyr::distinct()
            all_eQTL <- all_eQTL %>%
                dplyr::left_join(newprobes, 
                                 by = c("trait" = "probe")) %>%
                # color by distant or local (within 1 Mb of peak)
                dplyr::mutate(class = dplyr::case_when(probe_chr != chr ~ "diff_chr",
                                                       (probe_start + (probe_stop - probe_start)/2) - pos > 1e6 ~ "distant",
                                                       TRUE ~ "cis"))
        } else {
            all_eQTL <- NULL
        }
        
    })
    
    # text output for eQTL trait
    output$trait_name <- renderText({
        input$trait_input
    })
    
    #############################
    ####    EQTL PLOT        ####
    #############################
    output$eqtlplot <- plotly::renderPlotly({
        
        showModal(modalDialog(footer=NULL, size = "l", tags$div(style = "text-align: center;", "Loading...")))
        
        all_eQTL <- get_eQTL()
        
        if(is.null(all_eQTL)) {
            # plot blank if no eqtl for this trait (no qtl)
            plotly::ggplotly(ggplot2::ggplot(NULL))
        } else {
            # factor to order chr
            all_eQTL$chr <- factor(all_eQTL$chr, levels = c("I", "II", "III", "IV", "V", "X"))
            all_eQTL$probe_chr <- factor(all_eQTL$probe_chr, levels = c("X", "V", "IV", "III", "II", "I"))
            
            # plot eQTL peaks
            tsize <- 12
            eplot <- all_eQTL %>%
                # dplyr::mutate(n2_res = ifelse(var_exp > 0, "no", "yes")) %>%
                ggplot(.) +
                aes(x = pos / 1e6, y = lod, color = class, size = var_exp) +
                # geom_rect(data=df_chr_length, aes(xmin =  start/1e6, xmax = stop/1e6, ymin = 1, ymax=1.01), color='transparent', fill='transparent', size =0.1, inherit.aes =  F) +
                geom_point() +
                # scale_shape_manual(values = c("yes" = 24, "no" = 25), guide = FALSE) +
                facet_grid(~chr, scales = "free", space = "free") +
                theme_bw(tsize) +
                labs(x = "QTL position (Mb)", y = "LOD") +
                # scale_fill_manual(values = c("cis" = "grey", "distant" = "yellow", "diff_chr" = "red"), name = "Class") +
                scale_color_manual(values = c("cis" = "grey", "distant" = "yellow", "diff_chr" = "red"), name = "Class") +
                scale_size_continuous(range = c(0.1,2), guide = "none") +
                theme(panel.grid = element_blank(),
                      legend.position = "right",
                      axis.title = element_text(face = "bold"),
                      axis.text.y = element_text(face = "bold", color = "black"),
                      axis.text.x = element_text(face = "bold", color = "black", angle = 90),
                      strip.text = element_text(face = "bold")) +
                scale_alpha(guide = "none") +
                geom_vline(aes(xintercept = QTL/1e6), linetype = "dashed", color = "blue")
            
            removeModal()
            
            plotly::ggplotly(eplot +
                                 aes(text = glue::glue("Probe: {trait}\n Gene: {gene}\n Probe_pos: {probe_chr}:{round(probe_start/1e6, digits = 3)} Mb \n {ifelse(var_exp > 0, 'CB4856 resistant', 'N2 resistant')}")), 
                             tooltip = "text")
        }
    
    })
    
    #############################
    ####    EQTL DATA        ####
    #############################
    output$eqtl_data <- DT::renderDataTable({
        
        all_eQTL <- get_eQTL()
        
        if(is.null(all_eQTL)) {
            data.frame(probe = NA, gene = NA, eQTL_pos = NA, lod = NA, var_exp = NA, eff_size = NA, ci_l_marker = NA, ci_r_marker = NA, probe_pos = NA, class = NA)
        } else {
            # clean up and print
            all_eQTL %>%
                dplyr::mutate(probe_pos = paste0(probe_chr, "_", probe_start))%>%
                dplyr::select(probe = trait, gene, eQTL_pos = marker, lod, var_exp, eff_size, ci_l_marker, ci_r_marker, probe_pos, class)
        }
        
    })
    
    ########## QTL INPUT SHINY ###########
    output$qtl <- shiny::renderUI({
        
        # get inputs
        trt <- input$trait_input
        cond <- input$drug_input
        strainset <- input$set_input
        # interval <- input$intervals
        
        peaks <- loaddata()[[2]]
        
        # filter data
        traitmap <- peaks %>%
            dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                          set == strainset)
        
        # output
        shiny::selectInput(inputId = "whichqtl", label = "Select QTL:", choices = c("", unique(traitmap$marker)), selected = NULL)

    })
    
    #############################
    ####  GENES GENES GENES  ####
    #############################
    output$candidate_genes <- shiny::renderUI({
        
        source("query_genes.R")
        
        showModal(modalDialog(footer=NULL, size = "l", tags$div(style = "text-align: center;", "Loading...")))
        
        # first, get QTL
        qtl_marker <- input$whichqtl
        
        # get inputs
        trt <- input$trait_input
        cond <- input$drug_input
        strainset <- input$set_input
        # interval <- input$intervals
        
        peaks <- loaddata()[[2]]
        
        # filter data
        traitmap <- peaks %>%
            dplyr::filter(trait == glue::glue("{cond}.{trt}"),
                          set == strainset)
        
        # get the region from the marker
        markerdf <- traitmap %>%
            dplyr::filter(marker == qtl_marker) %>%
            dplyr::mutate(reg = paste0(chr, ":", ci_l_pos, "-", ci_r_pos))
        
        # next, run query genes
        # first element: dataframe,second element: text
        test <- query_genes(markerdf$reg)
        
        # organize data and add wormbase links
        df <- test[[1]] %>%
            # dplyr::group_by(wbgene) %>%
            # dplyr::mutate(go_term = paste(go_term, collapse = "; "),
            #               go_name = paste(go_name, collapse = "; "),
            #               go_description = paste(go_description, collapse = "; ")) %>%
            # dplyr::distinct() %>%
            dplyr::select(-go_term, -go_description, -gene_class_description, -gene_id, -go_annotation) %>%
            dplyr::rename(GO_term = go_name) %>%
            dplyr::ungroup() %>%
            # dplyr::mutate(WormBase = paste0("https://wormbase.org/species/c_elegans/gene/", wbgene)) %>%
            dplyr::mutate(wbgene = paste0("<a href='",paste0("https://wormbase.org/species/c_elegans/gene/", wbgene),"'>",wbgene,"</a>"))
        
        removeModal()

        ############# GENE DF #############
        output$dataframe <- DT::renderDataTable(escape = FALSE, {
            df
        })
        
        ############ GENE TEXT ############
        output$gene_text <- shiny::renderUI({

            # print as bullets
            tags$div(
                tags$ul(
                    tags$li(test[[2]][1]),
                    tags$li(test[[2]][2]),
                    tags$li(test[[2]][3]),
                    tags$li(test[[2]][4]),
                    tags$li(test[[2]][5]),
                    tags$li(test[[2]][6]),
                    tags$li(test[[2]][7])
                )
            )
        })
        
        # what to print to tab
        tagList(
            shiny::uiOutput("gene_text"),
            shiny::h3("List of genes:"),
            DT::dataTableOutput("dataframe")
        )

    })
    
    #############################
    ####  MARKDOWN DOWNLOAD  ####
    #############################
    output$code_download <- shiny::downloadHandler(
        
        # Dynamic file name -- WHY IS THIS ONLY PARTIALLY WORKING?
        filename = function() {
            # load inputs
            trt <- input$trait_input
            cond <- input$drug_input
            strainset <- as.numeric(input$set_input)
            glue::glue("{gsub('-', '', Sys.Date())}_{input$drug_input}_{input$trait_input}_set{input$set_input}_report.html")
            },
        
        content = function(file) {
            # Copy the report file to a temporary directory before processing it, in
            # case we don't have write permissions to the current working dir (which
            # can happen when deployed).
            tempReport <- file.path(tempdir(), "report.Rmd")
            file.copy("report.Rmd", tempReport, overwrite = TRUE)

            # load inputs
            trt <- input$trait_input
            cond <- input$drug_input
            strainset <- as.numeric(input$set_input)
            qtl_marker <- input$whichqtl
            
            # load data
            annotatedmap <- loaddata()[[1]] %>%
                dplyr::filter(set == strainset)
            allRIAILsregressed <- fst::read_fst("data/allRIAILsregressed.fst")
            data("eQTLpeaks")
            data("probe_info")
            linkagemapping::load_cross_obj("N2xCB4856cross_full")
            # load("data/newcross.Rda")
            gene_annotations <- fst::read_fst("data/gene_annotations.fst")
            
            # Knit the document - use local environment to keep all the ^ above variables ^
            rmarkdown::render(tempReport, output_file = file)
            
        }
    )
    
    # render help markdown
    output$help_md <- shiny::renderUI({
        shiny::HTML(markdown::markdownToHTML(knitr::knit('README.md', quiet = TRUE)))
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)

