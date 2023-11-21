# BiocManager::install("DEP")
library(here);library(data.table);library(dplyr);library('bayesbio')
library("pRolocdata");library(stringr);library(tibble);library(DEP);library(ggplot2)
output_folder = here::here('Output')
Hyperlopit_2018 = readr::read_tsv(here::here("Datasets","Processed","annot_hyperlopit.tsv"))
metabolic = openxlsx::read.xlsx(here::here('Datasets','Raw','metabolism_gene_list_Sabatini.xlsx')) |> unlist()
MS_data = openxlsx::read.xlsx(here::here('Datasets','Raw','MS protein ungrouped mean values.xlsx')) |> 
    dplyr::select(Accession,matches('^Abundance:\\.')) |> tibble::column_to_rownames('Accession') |> as.matrix()
names_MS_data = MS_data |> 
    colnames() |>
    stringr::str_extract(',\\.([:print:]{3,5})$') |> 
    stringr::str_remove_all('^,\\.')
names_MS_data = paste('Abundance',names_MS_data,rep(1:3,times = 5), sep = '_')
colnames(MS_data) = names_MS_data
HUMAN_9606 <- readr::read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% 
    purrr::set_names(c("Uniprot","Type","ID"))

dataset_name = "Breast_cancer_Garcia_chromatome"
input_matrix =MS_data %>% as.data.frame()
    set.seed(2023)
    experimental_design_DDA <-  data.frame(
        label = colnames(input_matrix),
        condition =  str_remove_all(colnames(input_matrix),"_[:digit:]*$"),
        replicate = str_remove_all(colnames(input_matrix),"^[:graph:]*_") %>% as.numeric()
    )
    data_unique_Etop <- input_matrix %>% rownames_to_column("name") %>% 
        left_join(HUMAN_9606 %>% 
                      subset(Type == "Gene_Name") %>% 
                      dplyr::select(-Type) %>% 
                      subset(!duplicated(Uniprot)), 
                  by = c("name" = "Uniprot")) %>%
        subset(!duplicated(name))
    Quant_columns <- which(colnames(data_unique_Etop) %in%colnames(input_matrix))# get LFQ column numbers
    data_se <- make_se(data_unique_Etop, Quant_columns, experimental_design_DDA)
    plot_frequency(data_se)+ggtitle(glue::glue("Protein_overlap ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_overlap ",dataset_name,".pdf")))
    data_filt <- filter_missval(data_se, thr = 1)
    plot_numbers(data_filt)+ggtitle(glue::glue("Protein_numbers ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_numbers ",dataset_name,".pdf")))
    plot_coverage(data_filt)
    ggsave(here::here(output_folder,glue::glue("plot_coverage ",dataset_name,".pdf")))
    data_filt@assays@data@listData[[1]][is.nan(data_filt@assays@data@listData[[1]])] <- NA 
    dev.off()
    pdf(here::here(output_folder,glue::glue("Protein_Missingness ",dataset_name,".pdf")), width = 10, height = 10) 
    plot_missval(data_filt)
    dev.off()
    # data_norm=data_filt
    data_norm <- normalize_vsn(data_filt)
    data_norm@assays@data@listData[[1]] <-  proDA::median_normalization(data_norm@assays@data@listData[[1]])
    # data_norm@assays@data@listData[[1]] <-  DEqMS::equalMedianNormalization(data_norm@assays@data@listData[[1]])
    
   
    
    # data_norm@assays@data@listData[[1]] <- input_matrix
    
     DEP::meanSdPlot(data_norm)
     ggsave(here::here(output_folder,glue::glue("normalize_vsn ",dataset_name,".pdf")))
    plot_normalization(data_se, data_norm)+ggtitle(glue::glue("Protein_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_normalisation ",dataset_name,".pdf")))
    
    
    dev.off()
        png(here::here(output_folder,glue::glue("Protein_Missingness_Abundance ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
        plot_detect(data_norm)
        dev.off()
        #ggsave(here::here(output_folder,glue::glue("Protein_missingness ",dataset_name,".pdf")))
        data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.1)
        multiple_imputation <-  imp4p::impute.mi(tab = data_norm@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                              conditions = experimental_design_DDA$condition %>% as.factor(),
                                              repbio = experimental_design_DDA$replicate %>% as.factor())
            
            
        rownames(multiple_imputation) <- rownames(data_norm@assays@data@listData[[1]])
        colnames(multiple_imputation) <- colnames(data_norm@assays@data@listData[[1]])
        data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
        data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm@assays@data@listData[[1]]),]
        
    plot_imputation(data_norm, data_imp)
    ggsave(here::here(output_folder,glue::glue("Protein_imputted ",dataset_name,".pdf")))
    pca_res <- prcomp(data_imp@assays@data@listData[[1]]  %>% na.omit() %>% t(), scale=FALSE)
    var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
    
    pca_res$x %>% 
        
        as.data.frame %>%
        rownames_to_column("Sample") %>% 
        # dplyr::select(!matches('T47D')) |> 
        mutate(Condition = str_remove(Sample,"_.$")) %>% 
        ggplot(aes(x=PC1,y=PC2, label = Sample, colour = Condition )) + geom_point(size=4) +
        ggrepel::geom_label_repel()+
        theme_bw(base_size=32) + 
        labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
             y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
        theme(legend.position="top") +
        ggtitle(dataset_name)+ 
        theme(plot.title = element_text(size = 20))
    ggsave(here::here(output_folder,glue::glue(dataset_name," PCA.pdf")))
    data_diff_all_contrasts <- DEP::test_diff(data_imp, "control", "Abundance_231")
    dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 0.5)
nucleotide_pathways = fread(here::here('Datasets','Raw','hsa_pathways.txt')) |> 
    subset(str_detect(Path_description, 'Purine')) |> pull(Path_id)
nucleotide_genes = fread(here::here('Datasets','Raw','KEGG_genes.csv')) |> 
    subset(pathway %in%nucleotide_pathways ) |> pull(ID)

    Comparisons_list <- list()
    for(i in (dep@elementMetadata %>% names() %>% str_subset("diff") )){
        n_accepted_NA_per_condition = 2
        contrast <- str_remove_all(i,"_diff")
        print(contrast)
        conditions <- contrast %>% str_match("([:graph:]*)_vs_([:graph:]*)$") %>% .[2:3]
        non_missing_in_all_comparison <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
                                               subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames(),
                                           data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
                                               subset(., rowSums(is.na(.))<=n_accepted_NA_per_condition) %>% rownames()) %>% unique()
        Imputted <- c(data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[1])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames(),
                      data_norm@assays@data@listData[[1]] %>% as.data.frame() %>% dplyr::select(contains(conditions[2])) %>% 
                          subset(., rowSums(is.na(.))>1) %>% rownames()) %>% unique()
        volcano_df <-  data.frame(log2_FC = dep@elementMetadata %>%  .[(glue::glue(contrast,"_diff"))] %>% unlist(),
                                  Uniprot = dep@elementMetadata$name,
                                  significant = dep@elementMetadata %>%  .[(glue::glue(contrast,"_significant"))] %>% unlist(),
                                  p.adj = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.adj"))] %>% unlist() ,
                                  p.val = dep@elementMetadata%>%  .[(glue::glue(contrast,"_p.val"))] %>% unlist()) %>% 
            subset(Uniprot %in% non_missing_in_all_comparison) %>% 
            
            mutate(Imputted_comparison = Uniprot %in% Imputted,
                
                   Single_Uniprot = Uniprot %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
            # left_join(Metabolic_proteins, by  = c("Single_Uniprot" = "Uniprot") ) %>% 
            # dplyr::rename(Metabolic_library = Behaviour) %>% 
            # left_join(Interesting_proteins, by = c("Single_Uniprot" = "Uniprot")) %>% 
            left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type), by  = c("Single_Uniprot" = "Uniprot")) #%>% 
        # mutate(Metabolic_library = if_else(is.na(Metabolic_library),"Non-Metabolic", Metabolic_library))
       qval =  qvalue::qvalue(volcano_df$p.val)
       volcano_df$qvalue = qval$qvalues
       GTP_enzymes = c('IMPDH1','GMPS','GUK1','^NME')
       volcano_df$significant = if_else(volcano_df$qvalue<0.05 & abs(volcano_df$log2_FC)>1,T,F)
       volcano_df= volcano_df |> mutate(
           category = case_when(
              ID =='IMPDH2' ~'IMPDH2',
             ID %in% nucleotide_genes ~ 'purine',
             ID %in% metabolic ~ 'metabolic',
             TRUE~'other'
           )
       )
       quantile_log2 = volcano_df |> 
           subset(Imputted_comparison==F) |> pull(log2_FC) |> abs() |> quantile() 
        volcano_df |> 
            subset(Imputted_comparison==F | ID == 'IMPDH2') %>% ggplot(aes(x = -log2_FC, y = -log10(p.val), 
                                                                           label = ID, fill = category))+
            geom_point(data = . %>% subset(significant == T & category == 'other'), alpha = 0.7, size = 3.9, shape =21, stroke=0, colour = 'grey50')+
            geom_point(data = . %>% subset(significant == F & category == 'other'), alpha = 0.2, size = 2.4, shape =21, stroke = 0, colour = 'grey50')+
            geom_point(data = . %>% subset(significant == T & category != 'other'), alpha = 0.7, size = 4, stroke =1, shape =21, colour = 'black')+
            geom_point(data = . %>% subset(significant == F & category != 'other'), alpha = 0.2, size = 2.5, stroke =1, shape =21, colour = 'black')+
            
            # geom_point(data = . %>% subset(category != 'other'), alpha = 0.2, size = 2.5, shape =21, colour = 'black')+
            geom_point(data = . %>% subset(significant == T & Uniprot == 'P12268'), shape =21,  stroke =1.5,size = 4, colour = 'black')+
            geom_point(data = . %>% subset(significant == F & Uniprot == 'P12268'), shape =21, stroke =1.5, size = 2.5, colour = 'black')+
            ggrepel::geom_text_repel(data = . %>% subset(significant == T & abs(log2_FC)>quantile_log2[4]))+
            ggrepel::geom_text_repel(data = . %>% subset(Uniprot == 'P12268'))+
            # ggrepel::geom_label_repel(data = . %>% subset(str_detect(ID,paste0(GTP_enzymes,collapse = '|'))),alpha = 0.8)+
            theme_bw()+scale_fill_manual(values = c('IMPDH2' = '#B8223A',
                                                      'purine' = '#EF633B',
                                                      'metabolic' = '#F7C736',
                                                      'other' = 'grey50'))+
            annotate("text", x = c(-2,2), y=0, label = str_remove_all(conditions,'Abundance_'))+
            
            ggtitle(glue::glue("Diff Present", contrast),
                    subtitle = dataset_name)
        ggsave(here::here(output_folder,glue::glue("Protein_volcano_significant",dataset_name," ",contrast,".pdf")), width = 7, height = 6)
        
        Comparisons_list[[i]] <- volcano_df
        
        
        # ego2 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T) ,
        #               OrgDb        = org.Hs.eg.db,
        #               ont          = "ALL",
        #               keyType = "UNIPROT",
        #               #nPerm        = 1000,
        #               minGSSize    = 50,
        #               maxGSSize    = 500,
        #               pvalueCutoff = 0.05,
        #               verbose      = FALSE)
        # if((ego2@result %>% nrow)>0){
        #     enrichplot::ridgeplot(ego2, showCategory = 68)+
        #         ggtitle(glue::glue(dataset_name," ALL-GSEA",i))
        #     ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," ALL-GSEA.pdf")), height = 20, width  = 15)}
        # # ego3 <- gseGO(geneList     = dep@elementMetadata %>% .[i] %>% unlist %>% set_names(dep@elementMetadata$name) %>% sort(decreasing = T),
        #               OrgDb        = org.Hs.eg.db,
        #               ont          = "BP",
        #               keyType = "UNIPROT",
        #               #nPerm        = 1000,
        #               minGSSize    = 50,
        #               maxGSSize    = 500,
        #               pvalueCutoff = 0.01,
        #               verbose      = FALSE)
        # if((ego3@result %>% nrow)>0){
        # 
        #     enrichplot::ridgeplot(ego3 , showCategory = 68)+
        #         ggtitle(glue::glue(dataset_name,"BP-GSEA",i))
        #     ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," BP-GSEA.pdf")), height = 20, width  = 15)}
        
        # gene_list <- volcano_df %>%
        #   dplyr::select(Single_Uniprot,log2_FC) %>%
        #   left_join(Human_hsa, by = c("Single_Uniprot" ="Uniprot")) %>%
        #   na.omit()%>%arrange(-log2_FC) %>%   pull(log2_FC,ID)
        # 
        # kk2 <- gseKEGG(geneList     = gene_list,
        #                organism     = 'hsa',
        #                minGSSize    = 10,
        #                pvalueCutoff = 0.05,
        #                verbose      = FALSE)
        # if((kk2@result %>% nrow)>0){
        # 
        #   enrichplot::ridgeplot(kk2, showCategory = 100)+
        #     ggtitle(glue::glue(dataset_name," KEGG ",i))
        #   ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," KEGG.pdf")), height = 20, width  = 15)}
        # 
        # 
        # mkk2 <- gseMKEGG(geneList = gene_list,
        #                  organism = 'hsa',
        #                  minGSSize = 10,
        #                  pvalueCutoff = 0.05)
        # if((mkk2@result %>% nrow)>0){
        # 
        #   enrichplot::ridgeplot(mkk2, showCategory = 68)+
        #     ggtitle(glue::glue(dataset_name," MKEGG ",i))
        #   ggsave(here::here(output_folder, glue::glue(i," ",dataset_name," MKEGG.pdf")), height = 20, width  = 15)}
        # 
        # 
    }
    names(Comparisons_list) <- names(Comparisons_list) |> str_remove_all('Abundance')
    openxlsx::write.xlsx(Comparisons_list,here::here('Datasets','Processed','Cell_line_comparison.xlsx'))
    
    listInput = list()
    for(i in names(Comparisons_list)){
        listInput[[i]]= Comparisons_list[[i]] |>
            subset(significant==T & Imputted_comparison==F) |> pull(ID)

    }
    number_significant = listInput |> unlist() |> table()
    upset(fromList(listInput), order.by = "freq")

regression_dt <-  data.table()
 for(i in rownames( data_imp@assays@data@listData[[1]])){
   regression_proteins = data_imp@assays@data@listData[[1]][i,] |> enframe() |> 
        mutate(
            cell_line = str_remove_all(name,'Abundance_') |> str_remove_all('_.$'),
            level = case_when(
                cell_line=='231'~5,
                cell_line=='BT474'~3,
                cell_line=='SKBR3'~4,
                cell_line=='MCF7'~1,
                cell_line=='T47D'~2,
                TRUE~20
                
        ))
    
    # create linear model
    linear_model <- lm( value ~level, regression_proteins)
    
    # print summary of linear model
    summary_dt = summary(linear_model)
    regression_dt_tmp = data.table(Uniprot =i ,
                                   pvalue = summary_dt$coefficients[2,4],
                                       Estimate =linear_model$coefficients[2] )
    regression_dt = rbind(regression_dt,regression_dt_tmp)
 }
fwrite(regression_dt,here::here('Datasets','Processed','agressiveness_regression.csv'))
regression_dt |> 
    inner_join(HUMAN_9606 |> subset(Type == 'Gene_Name')) |> 
    mutate(
        category = case_when(
            ID =='IMPDH2' ~'IMPDH2',
            ID %in% nucleotide_genes ~ 'purine',
            ID %in% metabolic ~ 'metabolic',
            TRUE~'other'
        ),
        padj = p.adjust(pvalue,method = 'BH')) |> 
    mutate(significant = padj<0.05 & ID %in% names(number_significant[number_significant>1])) |> 
    ggplot(aes(x = Estimate, y = -log10(pvalue),label = ID,fill = category))+
    geom_point(data = . %>% subset(significant == T & category == 'other'), alpha = 0.7, size = 3.9, shape =21, stroke=0, colour = 'grey50')+
    geom_point(data = . %>% subset(significant == F & category == 'other'), alpha = 0.2, size = 2.4, shape =21, stroke = 0, colour = 'grey50')+
    geom_point(data = . %>% subset(significant == T & category != 'other'), alpha = 0.7, size = 4, shape =21, stroke =1.5, colour = 'black')+
    geom_point(data = . %>% subset(significant == F & category != 'other'), alpha = 0.2, size = 2.5, shape =21,stroke =1.5, colour = 'black')+
    
    # geom_point(data = . %>% subset(category != 'other'), alpha = 0.2, size = 2.5, shape =21, colour = 'black')+
    geom_point(data = . %>% subset(significant == T & Uniprot == 'P12268'), shape =21, size = 4, colour = 'black')+
    geom_point(data = . %>% subset(significant == F & Uniprot == 'P12268'), shape =21, size = 2.5, colour = 'black')+
    theme_bw()+scale_fill_manual(values = c('IMPDH2' = '#B8223A',
                                                           'purine' = '#EF633B',
                                                           'metabolic' = '#F7C736',
                                                           'other' = 'grey50'))+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T ))+
    # ggrepel::geom_text_repel(data = . %>% subset(Uniprot == 'P12268'))+
    # ggrepel::geom_text_repel(data = regression_dt[Uniprot =='P12268',])+
    ggtitle('Correlation with Aggressiveness across all samples',
            subtitle = 'showing significant proteins (in correlation) which were also significant agaisnt 231 in at least 2 cell lines')
ggsave(here::here(output_folder,glue::glue("Linear regression aggressiveness ",dataset_name,".pdf")), width = 7, height = 6)
IMPDH2_etop = readxl::read_xls(here::here('Datasets','Raw','msb202211267-sup-0005-datasetev3.xls'), sheet = 4)
IMPDH2_etop |> 
    mutate(
        category = case_when(
            ID =='IMPDH2' ~'IMPDH2',
            ID %in% nucleotide_genes ~ 'purine',
            ID %in% metabolic ~ 'metabolic',
            TRUE~'other'
        )) |> 
     ggplot(aes(x = -log2_FC, y = -log10(p.val),label = ID,fill = category))+
    geom_point(data = . %>% subset(significant == T & category == 'other'), alpha = 0.7, size = 3.9, shape =21, stroke=0, colour = 'grey50')+
    geom_point(data = . %>% subset(significant == F & category == 'other'), alpha = 0.2, size = 2.4, shape =21, stroke = 0, colour = 'grey50')+
    geom_point(data = . %>% subset(significant == T & category != 'other'), alpha = 0.7, size = 4, shape =21, stroke =1.5, colour = 'black')+
    geom_point(data = . %>% subset(significant == F & category != 'other'), alpha = 0.2, size = 2.5, shape =21,stroke =1.5, colour = 'black')+
    
    # geom_point(data = . %>% subset(category != 'other'), alpha = 0.2, size = 2.5, shape =21, colour = 'black')+
    geom_point(data = . %>% subset(significant == T & Uniprot == 'P12268'), shape =21, size = 4, colour = 'black')+
    geom_point(data = . %>% subset(significant == F & Uniprot == 'P12268'), shape =21, size = 2.5, colour = 'black')+
    theme_bw()+scale_fill_manual(values = c('IMPDH2' = '#B8223A',
                                            'purine' = '#EF633B',
                                            'metabolic' = '#F7C736',
                                            'other' = 'grey50'))+
    ggrepel::geom_text_repel(data = . %>% subset(significant == T ))+
    # ggrepel::geom_text_repel(data = . %>% subset(Uniprot == 'P12268'))+
    # ggrepel::geom_text_repel(data = regression_dt[Uniprot =='P12268',])+
    ggtitle('Etop no release (left) against Etop 24hr release U2OS DNA damage')
ggsave(here::here(output_folder,glue::glue("Etop_DNA_damage_IMPDH2.pdf")), width = 7, height = 6)

    Significant_proteins <- data_imp@assays@data@listData[[1]] %>%
        as.data.frame() %>%
        rownames_to_column("ProteinGroup") |> 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>%
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>%
        subset(ID %in% c('IMPDH2','MTHFD1','CAD','UMPS','IMPDH1','DHODH','ATIC','APRT','GMPS','PAICS','HPRT1','PRPS2') ) |> 
        tidyr::pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>%
        mutate(Condition = str_remove_all(Condition,"_.$")) %>%
        ungroup() %>%
          group_by(ProteinGroup,Condition,ID) %>%
        dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = T)) %>%
        tidyr::pivot_wider(names_from = "Condition",values_from = "Mean_Abundance") %>%

        #mutate(duplicated = BiocGenerics::duplicated(Uniprot))
        ungroup %>%
        subset(!is.na(ID)) %>%
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>%
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>%
        mutate(ID = if_else(duplicated(ID),paste0(ID,"_1"),ID)) %>%


        column_to_rownames("ID") %>%
        dplyr::select(where(is.numeric)) %>%
        mutate(Rowmean = rowMeans(.),
               across(where(is.numeric),~.x- Rowmean)) %>%
        dplyr::select(-Rowmean) %>%
        as.matrix()
    paletteLength <- 30
    myColor <- colorRampPalette(c("grey90",'white', "#642669"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1),
                  seq(max(Significant_proteins, na.rm = T)/paletteLength, max(Significant_proteins, na.rm = T), length.out=floor(paletteLength/2)))

    pdf(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".pdf")), width = 7, height = 14)
    pheatmap::pheatmap(Significant_proteins[c('MTHFD1','IMPDH2','PRPS2','HPRT1','PAICS','GMPS','ATIC','APRT','IMPDH1','DHODH','CAD','UMPS'),
                                            c('Abundance_MCF7','Abundance_T47D','Abundance_BT474','Abundance_SKBR3','Abundance_231')],
                       cluster_cols = F,fontsize_row = 6, clustering_distance_rows = "euclidean", cluster_rows = F,
                       # scale = "row",
                       
                       main = glue::glue(dataset_name, " \nSignificant Proteins Normalised - imputted"),color=myColor, breaks=myBreaks)
    dev.off()
    # if(names(dev.cur()) != "null device"){dev.off()}
    #clusters <- NbClust::NbClust(Significant_proteins, method = "kmeans")$Best.partition
    
    
    
    purrr::imap(.x = Comparisons_list, ~.x  %>% select(Uniprot, log2_FC, significant,ID) %>% rename(
        "{.y}_FC" := log2_FC,
        "{.y}_significant" := significant)) %>% 
        purrr::reduce(full_join) %>% 
        ggplot(aes(x = c_c_vs_c_uc_diff_FC,
                   y = n_c_vs_n_uc_diff_FC, label = ID))+geom_point()+
        ggrepel::geom_label_repel(data = . %>% subset((c_c_vs_c_uc_diff_FC>1.5 &n_c_vs_n_uc_diff_FC<(-1.5)) |
                                                          (c_c_vs_c_uc_diff_FC<(-1.5) &n_c_vs_n_uc_diff_FC>1.5)))+
        ggtitle("Nuclear and Cytoplasmic changes upon Confinement")
    ggsave(here::here(output_folder, glue::glue(i," ",dataset_name,"Combined_Cyto_Nucle_FC.pdf")), height = 20, width  = 15)
    
    data_imp@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
        #glimpse() %>% 
        pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
        group_by(Condition,Uniprot) %>% 
        summarise(Mean_Abundance = mean(Abundance, na.rm = T),
                  Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
        ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
        geom_point()+
        facet_wrap("Condition")+
        ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
                subtitle = dataset_name)
    ggsave(here::here(output_folder,glue::glue("CV_conditions_imputted",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
    data_norm@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
        #glimpse() %>% 
        pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
        mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
        group_by(Condition,Uniprot) %>% 
        summarise(Mean_Abundance = mean(Abundance, na.rm = T),
                  Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
        ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
        geom_point()+
        facet_wrap("Condition")+
        ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
                subtitle =  dataset_name)
    ggsave(here::here(output_folder,glue::glue("CV_conditions_unimputted",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
 
    data_matrices <- list(Imputted = data_imp@assays@data@listData[[1]],
                          Unimputted = data_norm@assays@data@listData[[1]], 
                          DEPs = set_names(Comparisons_list, dep@elementMetadata %>% names() %>% str_subset("diff")))
    write.xlsx(data_matrices[1:2], here::here("Datasets","Processed",glue::glue(dataset_name, "matrices.xlsx")), overwrite = T,rowNames = T)
    write.xlsx(data_matrices[[3]], here::here("Datasets","Processed",glue::glue(dataset_name, "volcano_DFs.xlsx")), overwrite = T)
# Making proteomic ruler
Uniprot_length_Mass <- here::here("Datasets","Raw", "Uniprot_Molecular_sizes.tab") %>% 
    readr::read_tsv() %>% 
    dplyr::select(-`Entry name`) %>% 
    purrr::set_names(c("Uniprot","Length","Mass")) %>% 
    mutate(Mass = Mass/1000)
chromatome_Abudances = openxlsx::read.xlsx(here::here('Datasets','Raw','MS protein ungrouped mean values.xlsx')) |> 
    dplyr::select(Accession,matches('^Abundance:\\.')) |> as.data.frame() |> column_to_rownames('Accession') |> 
    dplyr::select(matches('47|MCF7|231')) |> 
    rowMeans() |> enframe(name = 'Accession',value = 'Chromat_breast')
HUMAN_9606 <- readr::read_tsv(here::here("Datasets","Raw", "HUMAN_9606_idmapping.dat"),
                       col_names = FALSE) %>% purrr::set_names(c("Uniprot","Type","ID"))
Proteomic_Ruler <- here::here("Datasets","Processed","CCLE_prot_Ruler.txt") %>% read.delim() %>% .[-1,] %>% 
    dplyr::select(matches("Copy|Uniprot_Acc|accuracy"))%>% 
    remove_rownames() %>% 
    column_to_rownames("Uniprot_Acc") %>% 
    #mutate(across(where(is.numeric), as.numeric)) %>% 
    purrr::set_names(.,str_remove_all(names(.),"Copy\\.number\\.")) %>% 
    mutate(across(contains("_"),~log2(as.numeric(.x))),
           across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
    # subset(!is.nan(U2OS_BONE)) %>% 
    # subset(.,rowSums(is.na(.))<(ncol(.)/3)) %>%
    subset(Absolute.quantification.accuracy != "low") %>%
    dplyr::select(-Absolute.quantification.accuracy) %>%
    janitor::clean_names()

    MS_data_counts = openxlsx::read.xlsx(here::here('Datasets','Raw','MS protein ungrouped mean values.xlsx')) |> 
        dplyr::select(Accession,matches('Peptides')) |> dplyr::select(-matches('Mascot')) 
    
    # MS_data_counts$condition = 'breast_chrom'
    MS_data_counts$Uniprot =MS_data_counts$Accession |> str_remove_all(';[:print:]*$')
    MS_data_counts = inner_join(MS_data_counts,Uniprot_length_Mass)
    MS_data_counts = inner_join(MS_data_counts,chromatome_Abudances)
    colnames(MS_data_counts) = c('protein_group','Total_pept','Unique_razor_pept','Uniprot','Length','Mass','Chromatome')
    MS_data_counts = MS_data_counts[,c(1,7,3,2,4,5,6)]
    readr::write_tsv(MS_data_counts,here::here("Datasets","Processed",glue::glue("For_Proteomic_ruler_",dataset_name,".txt")))
    

    input_matrix <- here::here("Datasets","Processed",'From_proteomic_ruler.txt') %>% 
        read.delim() %>% .[-1,] %>% 
        janitor::clean_names() %>% 
        subset(absolute_quantification_accuracy != "low") %>% 
        remove_rownames() %>% 
        dplyr::select(protein_group, contains("copy")) %>% as.data.frame() %>% 
        column_to_rownames("protein_group") %>%
        mutate(across(everything(), ~as.numeric(.x))) %>% 
        purrr::set_names(.,str_remove_all(colnames(.),"copy_number_")) %>% 
        mutate(across(where(is.numeric),~log2(as.numeric(.x))),
               across(where(is.numeric), ~if_else(is.infinite(.x), NaN,.x))) %>% 
        as.matrix() %>% 
        proDA::median_normalization()
    Breast_WCE = Proteomic_Ruler |> dplyr::select(matches('mcf7|t47d|skbr3|bt474|mdamb231')) |> 
        as.matrix() |> rowMeans() |> enframe(name = 'Uniprot',value = 'WCE_breast')
comparison_ruler = input_matrix |> as.data.frame()|> rownames_to_column('Uniprot') |> 
    inner_join(Breast_WCE)
comparison_ruler$chromatome= comparison_ruler$chromatome |> scale() %>% .[,1]
comparison_ruler$WCE_breast= comparison_ruler$WCE_breast |> scale()%>% .[,1]
comparison_ruler$Enrichment = comparison_ruler$chromatome-comparison_ruler$WCE_breast
Hyperlopit_2018 |> group_by(svm)
left_join(comparison_ruler,Hyperlopit_2018) |> 
    subset(`final.assignment`!='unknown') |> 
    # subset(!is.na(markers) & markers!='unknown') |> 
    ggplot(aes(x= `final.assignment`, y = Enrichment))+
    geom_boxplot()+
    theme_bw()+
    ggtitle('Relative enrichment of proteins in T47D MCF7 and MDAMB231 Samples')
ggsave(here::here(output_folder,glue::glue("Chromatin enrichment ",dataset_name,".pdf")))
library(clusterProfiler)
ego3 <- gseGO(geneList     = comparison_ruler |> pull(Enrichment,Uniprot) |> sort(decreasing = T),
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              keyType = 'UNIPROT',
              verbose      = FALSE)
ridgeplot(ego3)+ggtitle('GO CC Enrichment of proteomic-ruler chrom-vs-WCE')
ggsave(here::here(output_folder,glue::glue("Chromatin enrichment GOCC ",dataset_name,".pdf")))

# essentiality
sample_info <-  read.csv(fs::path_rel(here::here("Datasets",'Raw',"sample_info.csv")), stringsAsFactors = FALSE)

Achilles_file <- fs::path_rel(here::here("Datasets",'Raw',"Achilles_gene_effect.csv")) #Raw from https://depmap.org/portal/download/ 2/11/2020
Achilles <- inner_join(sample_info[,1:2],read.csv(Achilles_file), by = "DepMap_ID")[,-1] %>%
    column_to_rownames(var = "stripped_cell_line_name") %>% t() %>%
    magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2])
breast = sample_info |> subset(lineage == 'breast') |> dplyr::select(stripped_cell_line_name,lineage_molecular_subtype,lineage_sub_subtype,lineage_subtype )
GTP_enzymes = c('IMPDH1','IMPDH2','GMPS','GUK1','^NME')
Achilles_GTP = Achilles[str_detect(rownames(Achilles),paste(GTP_enzymes,collapse = '|')),] |> as.data.frame() |> 
    dplyr::select(any_of(breast$stripped_cell_line_name)) |> rownames_to_column('Gene') |> 
    tidyr::pivot_longer(-Gene,names_to = 'stripped_cell_line_name', values_to = 'Essentiality')
Achilles_GTP=Achilles_GTP |> inner_join(breast) |> 
    mutate(TBNC  = if_else(lineage_sub_subtype == 'ERneg_HER2neg','TNBC','non-TNBC')) 
    ggplot(Achilles_GTP,aes(colour = TBNC, y= Essentiality, x = TBNC))+
    
    geom_boxplot()+
    geom_jitter(alpha = 0.5, height = 0)+
    facet_wrap('Gene', scales = 'free_x')+
    ggtitle('Pattern of Essentiality of GTP-synthesis Genes ERneg_HER2neg vs nonERneg_HER2neg',
            subtitle = 'none of them are significant')+
    theme_bw()
ggsave(here::here(output_folder,glue::glue("GTP synthesis essentiality ",dataset_name,".pdf")))

# opencell 
open_cell_interactors = fread(here::here('Datasets','Raw','opencell-protein-interactions.csv')) |> 
    subset(target_gene_name %in% c('TOP2A','PARP1' )) |> 
    # inner_join(HUMAN_9606_idmapping |> subset(Type == 'Gene_Name'), by = c('interactor_gene_name' = 'ID') ) |> 
    mutate(
        category = case_when(
            interactor_gene_name =='IMPDH2' ~'IMPDH2',
            interactor_gene_name %in% nucleotide_genes ~ 'purine',
            interactor_gene_name %in% metabolic ~ 'metabolic',
            TRUE~'other'
        ))
fwrite(open_cell_interactors,here::here('Datasets','Processed','TOP2_PARP1_opencell_interactors.csv'))
overlap = open_cell_interactors$interactor_gene_name |> table()
overlap = overlap[overlap>1] |> names()
    ggplot(open_cell_interactors,aes(x=enrichment, y = pval,fill = category, label = interactor_gene_name))+
    geom_point(data = open_cell_interactors |> subset(category=='other'),alpha = 0.4, 
               shape = 21, colour = 'white', size = 2.5)+
        geom_point(data = open_cell_interactors |> subset(category!='other'),alpha = 1,
                   shape = 21, size = 3, colour = 'black')+
        geom_point(data = open_cell_interactors |> subset(category=='IMPDH2'),alpha = 1,
                   shape = 21, size = 3, colour = 'black')+
        # geom_point
        ggrepel::geom_text_repel(data = open_cell_interactors |> 
                                     subset(category != 'other'), size = 5,max.overlaps = 100, alpha = 0.3)+
        ggrepel::geom_text_repel(data = open_cell_interactors |> 
                                     subset((interactor_gene_name %in% overlap & category != 'other' )|
                                     enrichment>7.5), size = 5,max.overlaps = 100)+
      
    facet_wrap('target_gene_name',scales = 'free_y')+
    theme_bw()+scale_fill_manual(values = c('IMPDH2' = '#B8223A',
                                              'purine' = '#EF633B',
                                              'metabolic' = '#F7C736',
                                              'other' = 'grey50'))
    ggsave(here::here(output_folder,glue::glue("opencell_TOP2A_PARP1_other.pdf")), width = 9, height = 8)
    

