# CCLE IMPDH2 Sensitivity
library("tidyverse")
library("openxlsx")
### CPTAC Metastasis
#https://cptac-data-portal.georgetown.edu/study-summary/S039

#panel 1F, BRD4 and IMPDH2 between breast subtypes ##
S039_BRCA_propective_clinical_data <- read.xlsx(here::here('Datasets','Raw',"S039_BRCA_propective_clinical_data_r1.xlsx")) %>% 
    .[,c("Participant_ID","Estrogen_Receptor_Status_by_IHC","Progesterone_Receptor_Status_by_IHC","HER2_ERBB2_Status_by_IHC")] %>% na.omit()

S039_BRCA_propective_clinical_data<- S039_BRCA_propective_clinical_data %>%
    mutate(IHC = case_when(
        Estrogen_Receptor_Status_by_IHC == "Positive" & Progesterone_Receptor_Status_by_IHC == "Positive" & HER2_ERBB2_Status_by_IHC == "Positive" ~ "TP",
        Estrogen_Receptor_Status_by_IHC == "Positive" & Progesterone_Receptor_Status_by_IHC == "Negative" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "HR+",
        Estrogen_Receptor_Status_by_IHC == "Positive" & Progesterone_Receptor_Status_by_IHC == "Positive" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "HR+",
        Estrogen_Receptor_Status_by_IHC == "Negative" & Progesterone_Receptor_Status_by_IHC == "Positive" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "HR+",
        Estrogen_Receptor_Status_by_IHC == "Negative" & Progesterone_Receptor_Status_by_IHC == "Negative" & HER2_ERBB2_Status_by_IHC == "Positive" ~ "HER2+",
        Estrogen_Receptor_Status_by_IHC == "Negative" & Progesterone_Receptor_Status_by_IHC == "Negative" & HER2_ERBB2_Status_by_IHC == "Negative" ~ "TN",
        TRUE ~ "Other")) %>% 
    subset(!(IHC == "Other")) 

CPTAC_S039 <- read_tsv(here::here('Datasets','Raw',"CPTAC2_Breast_Prospective_Collection_BI_Proteome.tmt10.tsv")) %>%
    dplyr::select(matches("Gene|Ratio")) %>%
    pivot_longer(-Gene, names_to = "Sample", values_to = "Abundance" ) %>%
    subset(!str_detect(Sample,"Unshared")) %>%
    mutate(Sample = str_remove(Sample," Log Ratio"))
sample_ids <- read.xlsx(here::here('Datasets','Raw',"S039_Breast_Cancer_Prospective_Collection_Specimens.xlsx"))[,c(1:4)]
sample_ids <- left_join(sample_ids, S039_BRCA_propective_clinical_data[,c("Participant_ID","IHC")], 
                        by = c(`Participant.Protocol.Identifier.:.Collection.Protocol.Registration` = "Participant_ID")) %>% na.omit() %>%
    mutate(IHC = if_else(Sample.Type == "Adjacent_Normal", "Adjacent", IHC))%>%
    mutate(IHC = factor(IHC, levels = c("Adjacent","HR+","HER2+", "TP","TN")))



CPTAC_S039 <- CPTAC_S039 %>% left_join(sample_ids[,c("Specimen.Label","IHC")] %>% set_names(c("Sample","IHC")))

Genes_of_interest <- c("BRD4" ,"IMPDH2" , "EPCAM" ,"ERBB2" ,"ESR1" , "VIM","PGR")
CPTAC_S039_GOI<-CPTAC_S039 %>% subset(Gene %in% Genes_of_interest) %>% na.omit() %>% 
    mutate(Gene = factor(Gene,levels =Genes_of_interest ))

CPTAC_S039_GOI%>%
    ggplot((aes(x = IHC,y = Abundance, colour  =IHC)))+
    geom_boxplot()+geom_point(alpha = 0.5) + facet_wrap("Gene", ncol = 2)  +
    #      geom_rect(data = subset(tp,Gene %in% c('BRD4',"IMPDH2")),fill = "#D54FA4",xmin = -Inf,xmax = Inf, 
    #             ymin = -Inf,ymax = Inf,alpha = 0.1) + 
    theme_bw() +
    scale_color_manual(values = c(Adjacent = "#F7A8D7",
                                  `HR+` = "#FF9EE2",
                                  `HER2+` = "#FF67B0",
                                  TP = "#F62681",
                                  TN = "#BD0055"))+
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())+
    ggtitle("CPTAC proteomics data", subtitle =  "Breast cancer Patient samples")
ggsave(here::here('Output','CPTAC_IMPDH2_proteome.pdf'))

##### Cell Sensitivity to BRD4 inh ####
sample_info <-  read.csv(fs::path_rel(here::here('Datasets','Raw',"sample_info.csv")), stringsAsFactors = FALSE)
HUMAN_9606_idmapping <-
    readr::read_tsv(fs::path_rel(here::here('Datasets','Raw',"HUMAN_9606_idmapping.dat")),
                    col_names = FALSE) %>%
    setNames(c("Uniprot", "Type", "ID"))
##panel 1g##
loading_CCLE_TPM <- function(x,y){
    RNA_seq <- read_tsv(x) %>% dplyr::select(-transcript_ids) %>% 
        mutate(gene_id = str_remove_all(gene_id,"\\.[:graph:]*$")) %>%
        column_to_rownames("gene_id") %>% 
        set_names(str_remove_all(colnames(.),"_[:graph:]*$")) 
    RNA_seq <- RNA_seq[rowSums(RNA_seq == 0) <= (ncol(RNA_seq)*0.75), ] %>%
        .[rowSums(.)>100,] %>% 
        rownames_to_column("gene_id")%>%
        left_join(inner_join(y %>% 
                                 subset(Type == "Ensembl"), 
                             y %>% 
                                 subset(Type == "Gene_Name"), 
                             by = "Uniprot")  %>% 
                      dplyr::select(contains("ID")) %>% 
                      distinct() %>% 
                      set_names(c("gene_id","Gene_Name")), by = "gene_id")
    top_id <- data.frame(Gene_Name = RNA_seq$Gene_Name,
                         gene_id = RNA_seq$gene_id,
                         Abundance = RNA_seq %>% dplyr::select(-c(gene_id,Gene_Name)) %>%  
                             rowSums())%>%  group_split(Gene_Name) %>% 
        map_df(.x = .,~.x %>% arrange(desc(Abundance)) %>% head(1)) %>% pull(gene_id)
    RNA_seq <- RNA_seq %>% subset((gene_id %in% top_id) & !duplicated(Gene_Name)) %>% 
        na.omit() %>% dplyr::select(-gene_id)%>% remove_rownames()%>%
        column_to_rownames("Gene_Name")  %>% `+`(.,0.01)%>%
        `/`(.,matrixStats::rowMedians(as.matrix(.),na.rm = T)) %>% log2() %>% as.matrix
    RNA_seq
    
}
CCLE_RNA_seq_file <- fs::path_rel(here::here('Datasets','Raw',"CCLE_RNAseq_rsem_genes_tpm_20180929.txt"))#Raw from https://portals.broadinstitute.org/ccle/data 3/18/2020
CCLE_RNA_seq<-loading_CCLE_TPM(CCLE_RNA_seq_file,HUMAN_9606_idmapping)
extract_top_btm_expressing_cell_lines<-function(Gene_Name){
    #Gene_Name = "IMPDH2"
    ordered <- CCLE_RNA_seq[Gene_Name,]%>% unlist()%>%#.[names(.) %in% (sample_info %>% subset(lineage == "breast") %>% pull(stripped_cell_line_name))] %>% 
        sort(decreasing = T)
    list(
        top_cell_lines = head(ordered,round(length(ordered)/10)) %>% names(),
        btm_cell_lines = tail(ordered,round(length(ordered)/10)) %>% names()) 
}
extreme_cell_lines<- map(c("IMPDH2"),extract_top_btm_expressing_cell_lines) %>%
    set_names(c("IMPDH2"))
sec_screen_info <- read_csv(here::here('Datasets','Raw',"secondary-screen-replicate-collapsed-treatment-info.csv"))
DNA_inh <- sec_screen_info %>% subset(str_detect(target, "PARP|TOP")) %>% pull("broad_id") %>% paste(., collapse = "|")
sec_screen <- data.table::fread(here::here('Datasets','Raw',"secondary-screen-replicate-collapsed-logfold-change.csv")) %>% 
    inner_join(.,sample_info[,1:2], by = c("V1" = "DepMap_ID"))%>%
    #column_to_rownames(var = "stripped_cell_line_name")  %>% 
    dplyr::select(matches(DNA_inh)|contains("stripped"))  %>% 
    pivot_longer(cols = -stripped_cell_line_name, names_to = "column_name", values_to = "Sensitivity") %>%
    left_join(sec_screen_info)
extreme_cell_lines_IMPDH2<-rbind(data.frame(stripped_cell_line_name = extreme_cell_lines$IMPDH2$top_cell_lines,
                                           Expression = "High"),
                                data.frame(stripped_cell_line_name = extreme_cell_lines$IMPDH2$btm_cell_lines,
                                           Expression = "Low"))

sec_screen_DNAi <- sec_screen %>% inner_join(extreme_cell_lines_IMPDH2) 
for(comp_name in unique(sec_screen_DNAi$name)){
    sec_screen_DNAi_tmp = sec_screen_DNAi |> 
        subset(name ==comp_name ) |> 
        subset(dose == max(dose,na.rm = T)) 
    pvalue_test= t.test(sec_screen_DNAi_tmp |> subset(Expression=='Low') |> pull(Sensitivity),
                        sec_screen_DNAi_tmp |> subset(Expression=='High')|> pull(Sensitivity))
    sec_screen_DNAi_tmp |> ggplot(aes(x = Expression,y = Sensitivity,colour = Expression))+
        geom_boxplot()+theme_bw()+ theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
        geom_jitter(width = 0.1)+
        ggtitle(comp_name,
                subtitle = glue::glue('pvalue {pvalue_test$p.value}'))+
        scale_colour_manual(values = c('grey0','grey50'))
    ggsave(here::here('Output',glue::glue('{comp_name} IMPDH2 Expression.pdf')))
}