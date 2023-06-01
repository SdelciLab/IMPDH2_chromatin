# BiocManager::install("DEP")
library(here);library(data.table);library(dplyr);library('bayesbio')
library("pRolocdata");library(stringr);library(tibble);library(DEP);library(ggplot2)
output_folder = here::here('Output')
Hyperlopit_2018 = readr::read_tsv(here::here("Datasets","Processed","annot_hyperlopit.tsv"))
MS_data = openxlsx::read.xlsx(here::here('Datasets','Raw','MS protein ungrouped mean values.xlsx')) |> 
    dplyr::select(Accession,matches('Normalized')) |> tibble::column_to_rownames('Accession') |> as.matrix()
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
    pdf(here::here(output_folder,glue::glue("Protein_Missingness ",dataset_name,".pdf")), width = 10, height = 10) 
    plot_missval(data_filt)
    dev.off()
    data_norm=data_filt
    # data_norm <- normalize_vsn(data_filt)
    # data_norm@assays@data@listData[[1]] <-  proDA::median_normalization(data_norm@assays@data@listData[[1]])
    data_norm@assays@data@listData[[1]] <-  DEqMS::equalMedianNormalization(data_norm@assays@data@listData[[1]])
    
   
    
    # data_norm@assays@data@listData[[1]] <- input_matrix
    
     DEP::meanSdPlot(data_norm)
     ggsave(here::here(output_folder,glue::glue("normalize_vsn ",dataset_name,".pdf")))
    plot_normalization(data_se, data_norm)+ggtitle(glue::glue("Protein_norm ",dataset_name))
    ggsave(here::here(output_folder,glue::glue("Protein_normalisation ",dataset_name,".pdf")))
    
    
    
    if(data_norm@assays@data@listData[[1]] %>% is.na() %>% any()){
        png(here::here(output_folder,glue::glue("Protein_Missingness_Abundance ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
        plot_detect(data_norm)
        dev.off()
        #ggsave(here::here(output_folder,glue::glue("Protein_missingness ",dataset_name,".pdf")))
        data_imp <- DEP::impute(data_norm, fun = "MinProb", q = 0.01)
        multiple_imputation <-  imp4p::impute.mi(tab = data_norm@assays@data@listData[[1]],#methodMNAR = "impute.pa",
                                              conditions = experimental_design_DDA$condition %>% as.factor(),
                                              repbio = experimental_design_DDA$replicate %>% as.factor())
            
            
        rownames(multiple_imputation) <- rownames(data_norm@assays@data@listData[[1]])
        colnames(multiple_imputation) <- colnames(data_norm@assays@data@listData[[1]])
        data_imp@assays@data@listData[[1]] <- multiple_imputation#  rbind(to_not_impute,multiple_imputation)
        data_imp@assays@data@listData[[1]] <- data_imp@assays@data@listData[[1]][rownames(data_norm@assays@data@listData[[1]]),]}
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
    dep <- add_rejections(data_diff_all_contrasts, alpha = 0.05, lfc = 1)
    
 
    # data_imp@assays@data@listData[[1]] %>%
    #     as.data.frame() %>%
    #     rownames_to_column("ProteinGroup") %>%
    #     mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>%
    #     left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>%
    #     pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>%
    #     mutate(Condition= factor(Condition)) %>%
    #     group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>%
    #     mutate(DMSO = mean(c_across(contains("dmso")), na.rm = T)) %>%
    #     mutate(across(where(is.numeric), ~.x-DMSO)) %>%
    #     dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>%
    #     left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>%
    #     ggplot(aes(x = Condition, y  = Abundance, colour = ProteinGroup,group= ProteinGroup, label = ID, alpha= Significant))+
    #     geom_line()+
    #     geom_point()+
    #     scale_alpha_manual(values = c(0.3,1))+
    #     ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DDA$label,1) & Significant == T))+
    #     ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DDA$label,1) & Significant == F))+
    #     theme(legend.position = "none") +
    #     facet_wrap("Behaviour")+
    #     scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    #     ggtitle("Interesting DDR proteins Detected",
    #             "Significant - opaque, non-significant Transparent")
    # ggsave(here::here(output_folder,glue::glue("Known_Behaviour ",dataset_name,".pdf")), width = 20, height = 20)
    # data_norm@assays@data@listData[[1]] %>%
    #     as.data.frame() %>%
    #     rownames_to_column("ProteinGroup") %>%
    #     mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>%
    #     left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>%
    #     pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>%
    #     mutate(Condition= factor(Condition)) %>%
    #     group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>%
    #     mutate(Samples_median = median(c_across(where(is.numeric)), na.rm = T)) %>%
    #     mutate(across(where(is.numeric), ~.x-Samples_median)) %>%
    #     dplyr::select(-Samples_median) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>%
    #     left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>%
    #     subset(Significant == T) %>%
    #     ggplot(aes(x = Condition, y  = Abundance, colour = ProteinGroup,group= ProteinGroup, label = ID))+
    #     #geom_line()+
    #     geom_point()+
    #     #scale_alpha_manual(values = c(0.3,1))+
    #     #ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DDA$label,1) & Significant == T))+
    #     #ggrepel::geom_label_repel(data = . %>% subset(Condition == tail(experimental_design_DDA$label,1) & Significant == F))+
    #     theme(legend.position = "none") +
    #     facet_wrap("ID")+
    #     scale_x_discrete(guide = guide_axis(n.dodge = 2))+
    #     ggtitle(glue::glue("Significant - Proteins",dataset_name))
    # ggsave(here::here(output_folder,glue::glue("Significant_proteins ",dataset_name,".pdf")), width = 20, height = 20)
    Significant_protein_list <- data.frame(ProteinGroup = dep@elementMetadata$name,
                                           MCF7_231 = dep@elementMetadata$Abundance_MCF7_vs_Abundance_231_significant,
                                           T47D_231 = dep@elementMetadata$Abundance_T47D_vs_Abundance_231_significant,
                                           SKBR3_231 = dep@elementMetadata$Abundance_SKBR3_vs_Abundance_231_significant,
                                           BT474_231 = dep@elementMetadata$Abundance_BT474_vs_Abundance_231_significant) %>%
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>%
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot)))

    
    Significant_proteins <- data_imp@assays@data@listData[[1]] %>% 
        as.data.frame() %>% 
        rownames_to_column("ProteinGroup") %>% 
        #mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        #left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
        #left_join(Interesting_proteins) %>% 
        #subset(!is.na(Behaviour)) %>% 
        #pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #mutate(Condition= factor(Condition, levels= paste(rep(c("DMSO","T0","T24"), each= 3), rep(1:3,3),sep="_"))) %>% 
        #group_by(ProteinGroup) %>% pivot_wider(names_from = "Condition",values_from = Abundance) %>% 
        #mutate(DMSO = mean(c(DMSO_1,DMSO_2,DMSO_3))) %>% mutate(across(where(is.numeric), ~.x-DMSO)) %>%
        #dplyr::select(-DMSO) %>% pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        left_join(dep@elementMetadata$significant %>% set_names(dep@elementMetadata$name) %>% enframe(name = "ProteinGroup", "Significant")) %>% 
        subset(Significant == T) %>% dplyr::select(-Significant) %>% 
        pivot_longer(contains("_"), names_to = "Condition", values_to = "Abundance") %>% 
        #   mutate(Condition = str_remove_all(Condition,"_.")) %>% 
        # ungroup() %>% 
        #   group_by(ProteinGroup,Condition) %>% 
        # dplyr::summarise(Mean_Abundance = mean(Abundance, na.rm = T)) %>% 
        pivot_wider(names_from = "Condition",values_from = "Abundance") %>% 
        mutate(Uniprot= ProteinGroup %>% str_remove_all(";[:graph:]*$") %>% str_remove_all("-[:graph:]*$")) %>% 
        left_join(HUMAN_9606 %>% subset(Type == "Gene_Name") %>% dplyr::select(-Type) %>% subset(!duplicated(Uniprot))) %>% 
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
    paletteLength <- 20
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    # length(breaks) == length(paletteLength) + 1
    # use floor and ceiling to deal with even/odd length pallettelengths
    myBreaks <- c(seq(min(Significant_proteins, na.rm = T), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(max(Significant_proteins, na.rm = T)/paletteLength, max(Significant_proteins, na.rm = T), length.out=floor(paletteLength/2)))
    
    png(here::here(output_folder,glue::glue("Heatmap_Significant ",dataset_name,".pdf")), width = 2500, height = 3800,res  =300) 
    pheatmap::pheatmap(Significant_proteins[,experimental_design_DDA$label %>% sort()],cluster_cols = F,fontsize_row = 6, clustering_distance_rows = "euclidean", cluster_rows = T,
                       # scale = "row",
                       main = glue::glue(dataset_name, " \nSignificant Proteins Normalised - imputted"),color=myColor, breaks=myBreaks)
    dev.off()
    # if(names(dev.cur()) != "null device"){dev.off()}
    #clusters <- NbClust::NbClust(Significant_proteins, method = "kmeans")$Best.partition
    
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
        volcano_df %>% ggplot(aes(x = log2_FC, y = -log10(p.val), label = ID, colour = log2_FC, alpha = significant))+
            geom_point()+
            geom_point(data = . %>% subset(Imputted_comparison == T), colour = "grey50")+
            ggrepel::geom_label_repel(data = . %>% subset(significant == T))+
            annotate("text", x = c(-0.5,0.5), y=0, label = rev(conditions))+
            
            ggtitle(glue::glue("Diff Present", contrast),
                    subtitle = dataset_name)
        ggsave(here::here(output_folder,glue::glue("Protein_volcano_significant",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
        
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
    # replicates <- experimental_design_DDA %>% group_by(condition) %>% sample_n(2) %>% ungroup() %>% pull(label,condition)
    # input_df <- data_norm@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot")
    # named_vector = replicates[1:2]
    # rep_1 = named_vector[1]
    # rep_2 = named_vector[2]
    # MA_replicates(rep_1,rep_2,input_df)
    # 
    # MA_replicates <-  function(named_vector,named_vector_2,input_df){
    #     print(named_vector)
    #     print(named_vector_2)
    #     
    #     input_df %>% mutate(
    #        "{unique(names(named_vector))}" :=  diff({{named_vector}},{{named_vector_2}}))
    # }
    # 
    # data_norm@assays@data@listData[[1]] %>% as.data.frame() %>%  rownames_to_column("Uniprot") %>% 
    #     na.omit() %>% 
    #     mutate()
    #     #glimpse() %>% 
    #     pivot_longer(-Uniprot,names_to = "Condition", values_to = "Abundance") %>% 
    #     mutate(Condition = str_remove_all(Condition, "_.$")) %>% 
    #     group_by(Condition,Uniprot) %>% 
    #     summarise(Mean_Abundance = mean(Abundance, na.rm = T),
    #               Protein_CV = sd(Abundance, na.rm = F)/Mean_Abundance) %>% 
    #     ggplot(aes(x = Mean_Abundance, y = Protein_CV))+
    #     geom_point()+
    #     facet_wrap("Condition")+
    #     ggtitle("Coefficient of Variation in Conditions NAs in Sd retained",
    #             subtitle =  dataset_name)
    # ggsave(here::here(output_folder,glue::glue("CV_conditions_imputted",dataset_name," ",contrast,".pdf")), width = 10, height = 15)
    # 
    data_matrices <- list(Imputted = data_imp@assays@data@listData[[1]],
                          Unimputted = data_norm@assays@data@listData[[1]], 
                          DEPs = set_names(Comparisons_list, dep@elementMetadata %>% names() %>% str_subset("diff")))
    write.xlsx(data_matrices[1:2], here::here("Datasets","Processed",glue::glue(dataset_name, "matrices.xlsx")), overwrite = T,rowNames = T)
    write.xlsx(data_matrices[[3]], here::here("Datasets","Processed",glue::glue(dataset_name, "volcano_DFs.xlsx")), overwrite = T)
    
    return(data_matrices)
}