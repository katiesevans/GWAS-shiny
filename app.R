# library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
# library(linkagemapping)
library(shiny)
library(DT)
library(plotly)
library(fst)
library(curl)
library('R.utils')
library(RSQLite)

drugs <- c('2percDMSO','3percEthanol','cisplatin250','cisplatin500','vincristine',
           'bleomycin','G418','daunorubicin','docetaxel','etoposide','irinotecan',
           'dactinomycin','K12','OP50TRP','OP50B12','OP50','DA837','DA1877','HT115',
           'JUb71','JUb85','JUb87','abamectin','albendazole','amsacrine','bicuculline',
           'chlorpyrifos','fenbendazole-15','fenbendazole-30','mebendazole','MPTP',
           'arsenictrioxide','copper','deiquat','FUdR','nickel','silver','zinc',
           'chlorpromazine','monepantel','praziquantel','thymol','chlorothanil',
           'fluoxetine-125','thiabendazole-625','thiabendazole-125')

trts <- c('cv.EXT','cv.TOF','f.ad','f.L1','f.L2L3','f.L4','iqr.EXT','iqr.TOF',
          'mean.EXT','mean.norm.EXT','mean.TOF','median.EXT','median.norm.EXT',
          'median.TOF','n','norm.n','q10.EXT','q10.norm.EXT','q10.TOF','q25.EXT',
          'q25.norm.EXT','q25.TOF','q75.EXT','q75.norm.EXT','q75.TOF','q90.EXT',
          'q90.norm.EXT','q90.TOF','var.EXT','var.TOF')

phenodf <- fst::read_fst("data/96strain_pheno.fst")
genodf <- readr::read_tsv("~/Dropbox/AndersenLab/LabFolders/Katie/git/GWAS-shiny/data/Genotype_Matrix.tsv") %>%
    tidyr::gather(strain, allele, -(CHROM:ALT))


# Define UI for application
ui <- fluidPage(
   
   # Application title
   titlePanel("GWA Analysis - Andersen Lab"),
   
   # no sidebar layout
   shiny::wellPanel(
       # input well
       shiny::fluidRow(
           # conditions
           column(3, shiny::selectInput(inputId = "drug_input", label = "Condition:", choices = sort(drugs), selected = NULL)),
           # trait
           # columm(3, shiny::uiOutput("trait")),
           # this is not accurate, doesn't have traits for no QTL 
           column(3, shiny::selectInput(inputId = "trait_input", label = "Trait:", choices = c("", sort(trts)), selected = NULL)),
           # qtl marker
           column(3, shiny::uiOutput("qtl"))
       ),
       # code download button
       shiny::downloadButton(outputId = "code_download", label = "Get code")
   ),
      
  # Main panel for plots
  mainPanel(width = 12,
            
     shiny::tabsetPanel(type = "tabs", id = "test",
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
        fst_url <- glue::glue("https://raw.githubusercontent.com/katiesevans/GWAS-shiny/main/data/mappings/96strain/{cond}_pmd.fst")
        curl::curl_download(fst_url, tmp_file, mode="wb")
        
        # files
        annotatedmap <- fst::read_fst(tmp_file)
        
        ## right now it can show eigen on map but not choose it for region pxg
        peaks <- annotatedmap %>%
            dplyr::filter(peak_id)
        list(annotatedmap, peaks)
        
    })
    
    #############################
    ####    ALL PEAK PLOT    ####
    #############################
    # plot alllodplots for each condition (or multiple conditions)
    output$cond_plot <- shiny::renderUI({
        showModal(modalDialog(footer=NULL, size = "l", tags$div(style = "text-align: center;", "Loading...")))
        
        # get drug
        cond <- input$drug_input

        # load data
        peaks <- loaddata()[[2]]
        
        # plot all traits
        newmap <- peaks %>%
            arrange(desc(trait))
        
        newmap$trait <- factor(newmap$trait, levels = unique(newmap$trait))
        
        
        # get chromosome lengths
        chr_lens <- data.frame(CHROM = c("I", "II", "III", "IV", "V", "X"), 
                               start = rep(1,6), 
                               end = c(14972282,15173999,13829314,17450860,20914693,17748731),
                               condition = newmap$condition[1],
                               trait = newmap$trait[1])
        
        removeModal()
        
        ########### PLOT ############
        output$genplot <- shiny::renderPlot({
            
            # plot
            ggplot(newmap)+
                aes(x=peakPOS/1E6, y=trait)+
                theme_bw() +
                viridis::scale_fill_viridis(name = "log10(p)") + 
                viridis::scale_color_viridis(name = "log10(p)") +
                geom_segment(aes(x = startPOS/1e6, y = trait, xend = endPOS/1e6, yend = trait, color = log10p), size = 2, alpha = 1) +
                geom_segment(data=chr_lens, aes(x =  start/1e6, xend = end/1e6, y = trait, yend=trait), color='transparent', size =0.1) +
                geom_point(aes(color = log10p),size = 3, alpha = 1)+
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
                facet_grid(condition ~ CHROM, scales = "free_x", space = "free")
        })
        
        ############## PEAKS DF ###############
        output$peaksdf <- DT::renderDataTable({
            newmap %>%
                dplyr::select(CHROM, POS, trait, log10p, var.exp, startPOS, endPOS) %>%
                dplyr::arrange(desc(trait)) %>%
                dplyr::mutate(log10p = round(log10p, digits = 4),
                              var.exp = round(var.exp, digits = 4))
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
    ####      MAN PLOT       ####
    #############################
    # plot linkage defaults (LOD, pxg, riail pheno)
    output$allplots <- shiny::renderUI({
        
        showModal(modalDialog(footer=NULL, size = "l", tags$div(style = "text-align: center;", "Loading...")))
        
        # define variables
        trt <- input$trait_input
        cond <- input$drug_input

        # load data
        annotatedmap <- loaddata()[[1]]

        # filter data
        traitmap <- annotatedmap %>%
            dplyr::filter(trait == trt)
        
        pheno <- phenodf %>%
            dplyr::filter(condition == cond,
                          trait == trt)

        #########
        # Plots #
        #########
        
        ################### BAR PLOT WI ###################
        output$barplot <- plotly::renderPlotly({
            # riail pheno plot
            wipheno <- pheno %>%
                dplyr::group_by(strain, trait) %>%
                dplyr::mutate(avg_phen = mean(phenotype, na.rm = T),
                              col = dplyr::case_when(strain == "N2" ~ "N2",
                                                     strain == "CB4856" ~ "CB",
                                                     TRUE ~ "WI")) %>%
                dplyr::distinct(strain, trait, .keep_all = T) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(norm_pheno = ((avg_phen - min(avg_phen, na.rm = T)) / (max(avg_phen, na.rm = T) - min(avg_phen, na.rm = T)))) %>%
                dplyr::arrange(norm_pheno)
            
            wipheno$strain <- factor(wipheno$strain, levels = unique(wipheno$strain))
            
            # plot
            barplot <- wipheno %>%
                ggplot2::ggplot(.) +
                ggplot2::aes(x = strain, y = norm_pheno, fill = col, color = col) +
                ggplot2::geom_bar(stat = "identity") +
                ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB" = "blue", "WI" = "grey")) +
                ggplot2::scale_color_manual(values = c("N2" = "orange", "CB" = "blue", "WI" = "grey")) +
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
        
        ################# MAN & PXG PLOTS ###################
        output$manhatplot <- shiny::renderPlot({
            manp <- traitmap %>%
                dplyr::filter(CHROM != "MtDNA") %>%
                dplyr::mutate(aboveEIGEN = dplyr::case_when(aboveBF == T ~ "2",
                                                            aboveEIGEN == T ~ "1",
                                                            TRUE ~ "0")) %>%
                dplyr::distinct(marker, .keep_all = T) %>%
                ggplot2::ggplot(.) +
                ggplot2::aes(x = POS/1e6, y = log10p) +
                ggplot2::scale_color_manual(values = c("0" = "black", "1" = "deeppink2", "2" = "red")) +
                ggplot2::geom_rect(ggplot2::aes(xmin = startPOS/1e6,    # this is the plot boundary for LD and gene plots
                                                xmax = endPOS/1e6,    # this is the plot boundary for LD and gene plots
                                                ymin = 0, 
                                                ymax = Inf, 
                                                fill = "palevioletred1"), 
                                   color = "palevioletred1",fill = "palevioletred1",linetype = 2, alpha=.3)+
                ggplot2::geom_hline(ggplot2::aes(yintercept = EIGEN), color = "gray", alpha = .75, size = 1, linetype = 2) +
                ggplot2::geom_hline(ggplot2::aes(yintercept = BF), color = "gray", alpha = .75, size = 1) +
                ggplot2::geom_point(ggplot2::aes(color= factor(aboveEIGEN)), size = 0.7) +
                ggplot2::facet_grid( . ~ CHROM, scales = "free_x" , space = "free_x") +
                ggplot2::theme_bw(12) +
                ggplot2::theme(
                    axis.text = ggplot2::element_text(color = "black", face = "bold"),
                    axis.title = ggplot2::element_text(face = "bold", color = "black"),
                    strip.text = ggplot2::element_text(face = "bold", color = "black"),
                    plot.title = ggplot2::element_blank(),
                    panel.grid = ggplot2::element_blank(),
                    legend.position = "none",
                    panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
                ggplot2::labs(x = "Genomic position (Mb)",
                              y = expression(-log[10](italic(p))))
            
            if(nrow(traitmap %>% na.omit()) > 0) {
                pxgplot <- peaks %>% 
                    dplyr::filter(trait == trt) %>%
                    dplyr::left_join(pheno) %>%
                    dplyr::left_join(genodf) %>%
                    tidyr::drop_na(allele) %>%
                    dplyr::mutate(allele = factor(allele, levels = c(-1,1), labels = c("REF","ALT"))) %>% 
                    dplyr::mutate(norm_pheno = ((phenotype - min(phenotype, na.rm = T)) / (max(phenotype, na.rm = T) - min(phenotype, na.rm = T))),
                                  phenotype = norm_pheno) %>%
                    dplyr::group_by(allele) %>%
                    ggplot()+
                    aes(x = allele, y = phenotype, fill = allele)+
                    geom_jitter(size = 0.5, width = 0.1) +
                    geom_boxplot(alpha = 0.8, outlier.color = NA) +
                    scale_fill_manual(values = c("REF" = "grey", "ALT" = "palevioletred1")) +
                    ggplot2::facet_grid(~factor(marker, unique(peaks$marker)), scales = "free") + 
                    ggplot2::theme_bw(12) + 
                    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "bold", color = "black"), 
                                   axis.text.y = ggplot2::element_text(face = "bold", color = "black"), 
                                   axis.title.x = ggplot2::element_text(face = "bold", color = "black", vjust = -0.3),
                                   axis.title.y = ggplot2::element_text(face = "bold", color = "black"), 
                                   strip.text.x = ggplot2::element_text(face = "bold", color = "black"),
                                   strip.text.y = ggplot2::element_text(face = "bold", color = "black"),
                                   plot.title = ggplot2::element_blank(), 
                                   legend.position = "none", 
                                   panel.grid = element_blank(),
                                   panel.background = ggplot2::element_rect(color = NA, size = 0.6)) +
                    ggplot2::labs(x = "Genotype at QTL", y = trt)
                    
            } else {
                pxgplot <- ggplot2::ggplot(peaks) +
                    geom_blank() +
                    geom_text(x = 0.5, y = 0.5, label = "No QTL")
            }
            
            cowplot::plot_grid(manp, pxgplot, nrow = 2, align = "v", axis = "lr")
        })
        
        
        ################# PEAKS DF ################# 
        output$peaks <- DT::renderDataTable({
            traitmap %>%
                dplyr::filter(aboveEIGEN) %>%
                dplyr::arrange(CHROM, POS) %>%
                dplyr::mutate(thresh = dplyr::case_when(aboveBF == T ~ "BF",
                                                        aboveEIGEN == T ~ "EIGEN",
                                                        TRUE ~ "NS")) %>%
                dplyr::select(marker, CHROM, POS, log10p, peak_id, thresh, var.exp, startPOS, endPOS) %>%
                dplyr::mutate(log10p = round(log10p, digits = 4),
                              var.exp = round(var.exp, digits = 4))
        })
        
        removeModal()
        
        # what to show
        tagList(
            h3("QTL plots:"),
            plotly::plotlyOutput("barplot", height = "300px"),
            shiny::plotOutput("manhatplot", height = "600px"),
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

