require(ggplot2)
require(ggrepel)
require(ComplexHeatmap)
require(ggpubr)
require(poolr)
require(reshape2)
require(rstatix)
require(gginnards)
require(ggthemes)
require(ggeasy)
require(readxl)
require(colorRamp2)
require(circlize)
require(tidyverse)
require(gridExtra)
require(ggplotify)
require(EnhancedVolcano)
require(cowplot)
require(edgeR)
require(MAGeCKFlute)
require(ggExtra)

## -----------------------------------
## Process DepMap data
## -----------------------------------

depmap <- readRDS("./datasets/depmap/depmap_21Q3_screen_data.rds")
metadata = depmap$meta_data
syn_cell_lines <-
  rownames(depmap$meta_data[grep("Synovial", depmap$meta_data$subtype), ])
non_syn_cell_lines <-
  rownames(depmap$meta_data[-grep("Synovial", depmap$meta_data$subtype), ])

## -----------------------------------
## Genes to map 
## -----------------------------------

Genes_Syn = readxl::read_excel("./src/aml-2019/resources/Synovial Sarcoma List for Heatmaps Pramod.xlsx", sheet =1) %>%
  as.data.frame()
Genes_Syn = Genes_Syn$`Gene/Compound`
Genes_to_map = readxl::read_excel("./src/aml-2019/resources/Complexes_for_HeatMap.xlsx", sheet =2) %>%
  as.data.frame()
Genes_to_map = Genes_to_map$`Gene Names`
int = intersect(Genes_Syn, Genes_to_map)
Genes_added = c("SSX6P", "PCGF3", "SUMO2", "SS18", "BRD9", "SSX8P",
                "PIAS1", "TGFBRAP1", "SSX1", "RING1", "UBE2I", "AMT",
                "UBA2", "ATRX", "GEMIN2", "CCND2", "KAT6A", "MDM4",
                "USP7", "UBE2E1", "SSX3", "SRSF2", "SUMO1P3", "SMARCD1",
                "WDR5", "CTNNBIP1", "CTNNB1", "AKT1", "TADA2B", "SYNJ1", 
                "SAE1", "CBX4", "TNFRSF6B", "UTY", "CD9")
map_genes = c(Genes_to_map, Genes_added)
length(map_genes)
int_map = intersect(Genes_Syn,map_genes)

## -----------------------------------
## Process DepMap rnai gene fitness scores for genes of interest and Synovial sarcoma cell lines
## -----------------------------------

dependency_score = depmap$gene_fitness_rnai_depmap %>%
  as.data.frame()
ds_df = dependency_score %>%
  dplyr::filter(rownames(depmap$gene_fitness_rnai_depmap) %in% Genes)
ds_df_syn = ds_df %>%
  as.data.frame() %>%
  dplyr::select(all_of(syn_cell_lines))
ds_df_syn[is.na(ds_df_syn)] <- 0   
ds_df_syn = ds_df_syn %>% 
  dplyr::filter(!if_all(everything(), ~ . == 0))
ds_df_syn = ds_df_syn %>%
  select(where(~ any(. != 0)))

## -----------------------------------
## Figure 1C Heatmap
## -----------------------------------
col_fun = colorRamp2(c(-2, 0, 2), rev(c("red", "white", "#4682B4")))
ht2 = Heatmap(as.matrix(ds_df_syn),  cluster_rows = TRUE,
              cluster_columns = TRUE,
              row_labels = rownames(as.matrix(ds_df_syn)),
              show_column_dend = FALSE,
              show_row_dend = FALSE,
              column_title = "Synovial_cell_lines",
              row_names_gp = gpar(fontsize = 6),
              col = col_fun)

## reorder genes in the heatmap
clusters= draw(ht2)
gene_cluster = row_order(clusters)                       
clu_df <- lapply(1:length(gene_cluster), function(i){
  print(i)
  out <- data.frame(coordinates = rownames(ds_df_syn)[gene_cluster[[i]]],
                    Cluster = paste0("cluster", i), stringsAsFactors = FALSE)
  return(out)
}) %>%
  do.call(rbind, .)

clu_df_names = rev(clu_df$coordinates)
ds_df_syn_t_df = ds_df_syn[match(rev(clu_df$coordinates), rownames(ds_df_syn)),]
mark_at_1 = which(rownames(ds_df_syn_t_df) %in% map_genes)
mark_1 = rownames(ds_df_syn_t_df)[c(2,  3,  4, 5,  6,  7,  8,  9,
                                        10, 12, 13, 14, 16, 20,
                                        22, 23, 28, 29, 31, 32, 36,
                                        37, 38, 39, 42, 43, 47, 51,
                                        53, 56, 57, 64, 66, 71, 
                                        72, 73, 78, 90)]
ha_1 = rowAnnotation(foo = anno_mark(at = mark_at_1, labels = mark_1))
pdf(file = "heatmap_synovial_rnai_legends_genes_marked.pdf", height = 10, width = 10)
Heatmap(as.matrix(ds_df_syn_t_df),  cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        column_title = "Synovial_cell_lines",
        row_names_gp = gpar(fontsize = 6),
        col = col_fun,
        right_annotation = ha_1,
        row_labels = rownames(as.matrix(ds_df_syn_t_df)),
        row_names_side = "left")
dev.off()

## heatmap for nonsynovial sarcoma cell lines
ds_df_non_syn = ds_df %>%
  as.data.frame() %>%
  dplyr::select(all_of(non_syn_cell_lines))
ds_df_non_syn[is.na(ds_df_non_syn)] <- 0   
ds_df_non_syn = ds_df_non_syn %>%
  select(where(~ any(. != 0)))
ds_df_non_syn_t =t(apply(ds_df_non_syn,1,sort))
ds_df_non_syn_t_df = ds_df_non_syn[match(rownames(ds_df_syn_t_df), rownames(ds_df_non_syn)),]
ds_df_non_syn_t_S =t(apply(ds_df_non_syn_t_df,1,sort))
dat1 = ds_df_non_syn_t_S
dat1 = t(dat1)
dat1m <- melt(cbind(dat1, ind = rownames(dat1)), id.vars = c('ind'))
dd =dat1m %>% group_by(Var2) %>% 
  arrange(desc(Var1))%>% 
  mutate(pct = value- mean(value)/sd(value), right = cumsum(pct), left=lag(right, default=-1471), max = 10704)
g = ggplot(dd) + 
  geom_rect(aes(xmin=max, xmax=left, ymin=as.numeric(Var2)-.4, ymax=as.numeric(Var2)+.4, fill= ((value)))) + 
   scale_fill_distiller(palette = "RdBu")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme_void()
grob_1 = grid.grabExpr(draw(Heatmap(as.matrix(ds_df_syn_t_df),  cluster_rows = FALSE,
                                  cluster_columns = TRUE,
                                  show_row_names=F ,
                                  show_column_names = F,
                                  show_column_dend = FALSE,
                                  show_row_dend = FALSE,
                                  row_names_gp = gpar(fontsize = 6),
                                  col = col_fun,
                                  right_annotation = ha_1,
                                  row_labels = rownames(as.matrix(ds_df_syn_t_df)),
                                  row_names_side = "left"))) 

## combined heatmap of synovial and non synovial cell lines
p2 = plot_grid(g, grob_1)
file <- paste0("Synovial_rnai_color_heatmap", fileext = ".pdf")
save_plot(file, p2, ncol = 2, base_asp = 0.5, base_height = 12)

## -----------------------------------
## Figure 1E Square Plot
## -----------------------------------

## read square plot input data
df_squareview = read.csv("./dataframe_squareview_rnai.csv", header = T) 
df_squareview$X = NULL
non_int_map = setdiff(Genes_Syn,int_map)
df_squareview_1 = df_squareview %>%
  mutate(map = ifelse((Gene %in% int_map), "mapped",
                         ifelse((Gene %in% non_int_map), "unmapped", "genes")))
data = df_squareview_1
data = as.data.frame(data, stringsAsFactors = FALSE)
data = data[!(is.na(data[,"logfc"])|is.na(data[,"CERES"])), ]
data$Label = data$Gene
g = ggplot(data=df_squareview_1, aes_string(x="logfc", y="CERES"))+
 geom_point(size = 1, alpha = 0.8, colour="gray") +
  theme_classic()
data1 = df_squareview_1[df_squareview_1$map == "mapped",]

## ggplot square view 
g <- g + geom_point(
  data = data1,
  aes_string(x = "logfc", y = "CERES"),
  colour = "#1034A6",
  fill = "#1034A6",
  size = 2,
  alpha = 1
) 
g1 = g+
  geom_label_repel(
    data = data1,
    label = data1$Gene,
    check_overlap = F,
    nudge_x = 0.05,
    nudge_y = 0.05
   
  )
p = g1 + labs(y = "Difference in Average CERES score in Synovial vs non-Synovial cell lines", 
             x = "log2FC diff expression in Synovial vs non-Synovial cell lines") 
p = ggMarginal(p, type="histogram", binwidth = 0.2, size=10, fill = "#bdbdbd")

pdf(file = "Synovial_rnai_square_plot_mmaped_genes.pdf", height = 20, width = 15)
ggMarginal(p, type="histogram", binwidth = 0.2, size=10, fill = "#bdbdbd")
dev.off()


## -----------------------------------
## Figure 2B, 2H, 2I and 2J
## -----------------------------------
Sweave("Day14_vs_Day0_sample_rra_summary.Rnw");


## -----------------------------------
## Figure 2C
## -----------------------------------
Sweave("Mouse_vs_Day0_sample_rra_summary.Rnw");


## -----------------------------------
## Figure 2F, 2H
## -----------------------------------
## Read MAGeCKMLE output table
count_table = read.table("./SySa-DepMap_co_ctrl_sample.count_normalized.txt", header = T)
file3 = "./SySa.mle.sample.gene_summary.txt"
gdata = ReadBeta(file3)
gdata$HumanGene = TransGeneID(gdata$Gene, fromType = "symbol", toType = "symbol",
                              fromOrg = "mmu", toOrg = "hsa")
gdata_cc = NormalizeBeta(gdata, id = "HumanGene", samples=c("Day14", "Mouse"), 
                         method="cell_cycle")

## Run FluteMLE 
FluteMLE(gdata_cc, treatname="Mouse", ctrlname="Depmap", proj="MLE_Mouse_Depmap_gene_Summary_hsa_cc", organism="hsa", incorporateDepmap = TRUE)
FluteMLE(gdata, treatname="Day14", ctrlname="Depmap", proj="MLE_Day14_Depmap_gene_Summary_hsa", organism="hsa", incorporateDepmap = TRUE)


## -----------------------------------
## Figure S2 
## -----------------------------------
Sweave("SySa-DepMap_co_ctrl_sample_countsummary.Rnw");





