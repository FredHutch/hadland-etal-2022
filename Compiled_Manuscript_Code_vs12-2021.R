### Project: Compiled R Scripts for Nature Communications manuscript, Hadland et al 2021; updated 2021-12-23
rm(list = ls())
options(warn=-1) 
###load libraries
library(monocle)
library(VGAM)  
library(viridis)                   
library(dplyr)
library(mygene)  
library(ggpubr)
library(scales)
sessionInfo()

###define file path
RES_DIR <- file.path("~/...")     	### change the file path here to match the folder where RDS and related files are saved locally
									### Monocle cds, 10X output files, and fastq files available at NCBI GEO, accession #GSE145886


### AGM-EC scRNAseq data analysis###
###LOAD UNPROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "SAMPLE1-4.RDS"))

agm <- updateCDS(agm)    
agm <- estimateSizeFactors(agm)
agm <- estimateDispersions(agm)
disp_table = dispersionTable(agm)
disp_table = disp_table %>% mutate(excess_disp =
  (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 2500)$gene_id)
agm = setOrderingFilter(agm, top_subset_genes)
agm <- preprocessCDS(agm, num_dim = 11)                   
agm <- reduceDimension(agm, reduction_method = 'UMAP')
agm <- clusterCells(agm,method = 'louvain', res = 1e-6, verbose = T)

###LOAD PROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "SAMPLE1-4_P.RDS.RDS"))    

### simple plotting theme removes axes elements, text, and legend for UMAP plots
simple_theme <-  theme(text = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank())

### UMI and genes plotted = Sup Fig 1a
agm <- detectGenes(agm, min_expr=0.1)
pData(agm)$UMI <- Matrix::colSums(exprs(agm))
summary(pData(agm)$UMI)
summary(pData(agm)$num_genes_expressed)
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
        panel.grid = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size=15))
sample_color <- c("EC4" = "navy", "EC2" = "darksalmon", "EC1" = "tomato", "EC3" = "firebrick4")
ggplot(pData(agm), aes(x=sample, y=UMI, fill=sample)) +geom_boxplot() + scale_y_continuous(trans='log10', limits=c(100,100000), labels=number) + simple_theme2 +
                   scale_fill_manual(values = sample_color)
ggplot(pData(agm), aes(x=sample, y=num_genes_expressed, fill=sample)) +geom_boxplot() +
                  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 +
                   scale_fill_manual(values = sample_color) 


#### plot UMAP by sample = Fig 1c (left panel)
sample_color2 <- c("4_EC4" = "navy", "2_EC2" = "darksalmon", "1_EC1" = "tomato", "3_EC3" = "firebrick4")
plot_cell_clusters(agm,
                   color_by = 'sample',
                   cell_size = 0.8,
                   show_group_id = F) +
                   scale_color_manual(values = sample_color2) +simple_theme

#### plot UMAP by clusters = Fig 1c (right panel)
plot_cell_clusters(agm,
                   color_by = 'Cluster',
                   cell_size = 0.8,
                   show_group_id = F) +simple_theme



#### Identify cluster-specific genes = Supplementary Data 1
start <- Sys.time()
spatial_res <- principalGraphTest(agm, relative_expr = TRUE, k = 25, cores = detectCores() - 2, verbose = FALSE)
end <- Sys.time()
end - start
cluster_marker_res <-
  find_cluster_markers(agm,
                       spatial_res,
                       group_by = 'Cluster',
                       morans_I_threshold = 0.25) 
genes <-
  (cluster_marker_res %>%
  dplyr::filter(mean > 0.2, percentage > 0.2, specificity >0.75) %>%  
  dplyr::group_by(Group) %>%
  dplyr::top_n(1000, wt = specificity))
S_top_genes <- genes %>%
               dplyr::filter(Group == 1) %>%
               dplyr::top_n(n = 500, wt  = specificity)
write.csv (S_top_genes, file.path(RES_DIR, "P_AGMEC_S_top_genes.csv"))  ### Save csv file of top genes in supportive AGM-EC
NS_top_genes <- genes %>%
               dplyr::filter(Group == 2) %>%
               dplyr::top_n(n = 500, wt  = specificity)
write.csv (NS_top_genes, file.path(RES_DIR, "P_AGMEC_NS_top_genes.csv")) ### Save csv file of top genes in non-supportive AGM-EC


### Gene set-scores for GO gene sets identified in supportive vs non-supportive AGM-EC = Fig 1d
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  return(cds)
}
go1_genes <- c("Cyr61", "Dab2", "Cd44", "Thbs1", "Cxcl12", "Gata3", "Icam1", "Itga3", "Vcam1", "Fbln2", "Selp", "Fn1", "Bcl10", "Csf1", "Col8a1") #positive cell adh
go2_genes <- c("Gstp1", "Nes", "Birc5", "Hmgb2", "Cyr61", "Dab2", "Cd44", "Thbs1", "Cxcl12", "Gata3", "Ptprz1", "Icam1", "Edn1", "Slc7a11", "Mt1", "Timp1", "Fn1", "Bcl10", "Ctgf", "Ednrb", "Lgal3", "Csf1", "Plaur", "Sncg", "Lifr") #neg reg cell death
go3_genes <- c("Cyr61", "Dab2", "Thbs1", "Cxcl12", "Icam1", "Itga3", "Vcam1", "Fn1", "Ctgf", "Ptprz1") #integrin binding
go4_genes <- c("Nes", "Fst", "Hmgb2", "Cyr61", "Dab2", "Cd44", "Thbs1", "Cxcl12", "Vcam1", "Gata3", "Ptprz1", "Fn1", "Icam1", "Csf1", "Edn1", "Itga3", "Timp1", "Fabp4", "Rnd1", "Nrg1", "Ctg1", "Ednrb", "Lgals3", "Plaur", "Lifr", "Srm") #signaling receptor binding
agm <- estimate_score(agm, markers = go1_genes)   ### repeat for each gene set above = Fig 1d
plot_cell_clusters(agm,
           color_by = 'score',
           cell_size = 0.8,
           show_group_id = F) +
    scale_color_gradient2(low="gray80", mid = "gray80", high="red3", midpoint = (max(pData(agm)$score)/3), space = "Lab") +simple_theme


### Plot gene x expression heatmap on UMAP = Fig 1e (Genes= Cdh5, Sox17, Tbx20, Fbln2, Icam1 ,Tgm2, Cola5, Fn1, Cxcl12, Mki67), Sup Fig 1b (Gene = Mki67)
plot_cell_clusters(agm,                 
                   markers = c('x'),
                   color_by=marker, 
                   show_group_id = T, cell_size = 0.8)+
            scale_color_viridis(option = "magma") +simple_theme


###### Identify expressed genes in supportive AGM-EC for ligand analysis = Supplementary Data 5, Tab C
agm1 <- agm[,pData(agm)$sample != "EC4"]   ### remove non-supportive EC
table(pData(agm1)$sample)
agm1 <- detectGenes(agm1, min_expr = 0.1)                                              
fData(agm1)$use_for_rec_list <- fData(agm1)$num_cells_expressed > 0.1 * ncol(agm1)
agm1_genes_common <- row.names(subset(fData(agm1), fData(agm1)$use_for_rec_list == TRUE))
gene.list = agm1_genes_common
agm1_genes_common <- getGenes(gene.list, fields='symbol', species = 'mouse')
RecLig = read.csv(file.path(RES_DIR, "RecLigPairs.csv"))            # read csv file RecLig of published database of Receptor Ligand pairs
AgmEC_Ligands <- intersect(RecLig$Ligand, agm1_genes_common$symbol)  # create list of Ligands
write.csv(AgmEC_Ligands, file.path(RES_DIR, "AgmEC_Ligands.csv"))   # save csv file Supportive AGM-EC ligands










### AGM E10 to E11 V+61+E+ scRNAseq data analysis ###
###LOAD UNPROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample5-8.RDS"))
### define embryonic stage by sample source
pData(agm)$sample <- plyr::revalue(as.character(pData(agm)$sample),
                                        c("agm3_SI-GA-C11" = 'E10',
                                          "agm4_SI-GA-C12" = 'E10',
                                          "agm1_SI-GA-B9" = 'E11',
                                          "agm2_SI-GA-B10" = 'E11'))
agm <- updateCDS(agm)   
agm <- estimateSizeFactors(agm)
agm <- estimateDispersions(agm)
disp_table = dispersionTable(agm)
disp_table = disp_table %>% mutate(excess_disp =
  (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 1000)$gene_id)
agm = setOrderingFilter(agm, top_subset_genes)
agm <- preprocessCDS(agm, num_dim = 6)                   
agm <- reduceDimension(agm, reduction_method = 'UMAP')   
agm <- clusterCells(agm, method = 'louvain', res = 4e-4, verbose = T)

###LOAD PROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample5-8_P.RDS"))

### UMI and genes plotted = Sup Fig 5d
agm <- detectGenes(agm, min_expr=0.1)
pData(agm)$UMI <- Matrix::colSums(exprs(agm))
summary(pData(agm)$UMI)
summary(pData(agm)$num_genes_expressed)
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
        panel.grid = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size=15))
ggplot(pData(agm), aes(x=sample, y=UMI, fill=sample)) +geom_boxplot() +  scale_y_continuous(trans='log10', limits=c(100,200000), labels=number) + simple_theme2 
ggplot(pData(agm), aes(x=sample, y=num_genes_expressed, fill=sample)) +geom_boxplot()  +
                  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 


### plot by Cluster in UMAP = Fig 3a
plot_cell_clusters(agm,
                   color_by = 'Cluster',
                   cell_size = 0.8,
                   show_group_id = F) +simple_theme

### Plot gene x expression in UMAP = Fig 3b (genes = Cdh5, Dll4, Runx1, Gfi1)
### Sup Fig 6 (HSC-primed HE genes), Sup Fig 7 (AGM HSC signature genes), Sup Fig 9c (Supportive AGM-EC genes)
plot_cell_clusters(agm,                 
                   markers = c('x'),                            
                   color_by=marker, show_group_id = F,
                   cell_size = 0.8) + scale_color_viridis(option = "magma") +simple_theme 

### plot by embryo stage = Sup Fig 5e
plot_cell_clusters(agm,
                   color_by = 'sample',
                   cell_size = 1,
                   show_group_id = F)  +
                  simple_theme


### Classify cells as arterial EC, non-arterial EC, hematopoietic, somite = Figure 3c
Cdh5_id <- row.names(subset(fData(agm), gene_short_name == "Cdh5"))
Efnb2_id <- row.names(subset(fData(agm), gene_short_name == "Efnb2"))
Dll4_id <- row.names(subset(fData(agm), gene_short_name == "Dll4"))
Hey1_id <- row.names(subset(fData(agm), gene_short_name == "Hey1"))
Runx1_id <- row.names(subset(fData(agm), gene_short_name == "Runx1"))
Gfi1_id <- row.names(subset(fData(agm), gene_short_name == "Gfi1"))
Nr2f2_id <- row.names(subset(fData(agm), gene_short_name == "Nr2f2"))
Nrp2_id <- row.names(subset(fData(agm), gene_short_name == "Nrp2"))
Meox1_id <- row.names(subset(fData(agm), gene_short_name == "Meox1"))
Meox1_id <- row.names(subset(fData(agm), gene_short_name == "Meox1"))
Pax1_id <- row.names(subset(fData(agm), gene_short_name == "Pax1"))
Myf5_id <- row.names(subset(fData(agm), gene_short_name == "Myf5"))
Pax3_id <- row.names(subset(fData(agm), gene_short_name == "Pax3"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "EC2Art", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Dll4_id,] > 0 & x[Efnb2_id,] > 0 & x[Hey1_id,] > 0 & x[Nr2f2_id,] ==0 & x[Nrp2_id,] == 0 & x[Runx1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "EC1Ven", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Nr2f2_id,] > 0 & x[Nrp2_id,] > 0 & x[Meox1_id,] == 0 & x[Runx1_id,] == 0 & x[Efnb2_id,] == 0 & x[Hey1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "Hem", classify_func=function(x) {x[Gfi1_id,] > 0 & x[Runx1_id,] > 0})
cth <- addCellType(cth, "Som", classify_func=function(x) {(x[Meox1_id,] > 0 & x[Pax1_id,] >0) | (x[Pax3_id,] > 0 & x[Myf5_id,] > 0)})
agm <- classifyCells(agm, cth, method="markers-only")
cell_type_color <- c("1_EC1Ven" = "skyblue",
                     "2_EC2Art" = "royalblue3",
                     "3_Hem" = "red3",
                    "4_Som" = "burlywood4", 
                    "5_Unknown" = tgray)
### plot by CellType in UMAP = Figure 3c
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme


#### Pseudotime trajectory analysis = Fig 3d
agm <- partitionCells(agm)
agm <- learnGraph(agm,  RGE_method = 'SimplePPT')
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
root_node_ids = get_correct_root_state(agm,
                                      cell_phenotype =
                                        'Cluster', "2")
agm <- orderCells(agm, root_pr_nodes = root_node_ids)

### remove clusters not included in EHT psuedotime
agmx <- agm[,pData(agm)$Cluster != "3"]
agmx <- agmx[,pData(agmx)$Cluster != "4"]
agmx <- agmx[,pData(agmx)$Cluster != "6"]
agmx <- agmx[,pData(agmx)$Cluster != "7"]
plot_cell_trajectory(agmx, cell_size = 0.8) + 
    scale_color_viridis(option = "viridis") +simple_theme


### heatmap of AEC and HSC-related marker genes by pseudotime = Fig 3e
markers <- c("Procr", "Flt3", "Efnb2", "Nrp1", "Myb", "Hlf", "Kit", "Pdzk1ip1", "Sox17", "Fgd5", "Adgrg1", "Mecom", "Prom1", "Gata2", "Foxc2", "Il7r", "Gfi1", "Runx1", "Spi1", "Hlf", "Csf3r", "Mpo", "Dll4", "Esam", "Mycn", "Itga2b", "Spn", "Ptprc", "Tek", "Pbx1", "Tal1", "Meis1", "Etv6", "Ikzf2", "Nfe2", "Lyl1", "Cd44", "Pbx1", "Bmx", "Vwf", "Hhex", "Ace", "Sox7", "Gja4", "Gja5", "H19", "Ctnnal1", "Mllt3", "Jam3", "Gata3", "Fcgr3", "Hoxb5", "Hif1a", "Foxo1", "Kdm6b", "Bmi1", "Zfp521", "Ptk7", "Sirt6", "Lin28b", "Erg", "Fubp1", "Ssbp2", "Smurf2", "Ash1l", "Cited2", "Ezh1", "Cul4a", "Hoxb4", "Id1", "Id2", "Kmt2a", "Hoxa7", "Hoxa9")
marker_genes <- row.names(subset(fData(agmx),
                                 gene_short_name %in% markers))
diff_test_res <- differentialGeneTest(agmx[marker_genes,],
                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
pal <- colorRampPalette(c("navy", "white", "red"))
plot_pseudotime_heatmap(agmx[sig_gene_names,], hmcols = pal(500), show_rownames = TRUE)


####### global analysis of genes DE over pseudotime = Supplementary Data 3
agmx <- detectGenes(agmx, min_expr = 1)  
diff_test_res <- differentialGeneTest(agmx, fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res <-subset(diff_test_res, num_cells_expressed > (0.05 * ncol(agmx)))   ## select genes expressed in >5% cells
diff_test_res <-subset(diff_test_res, qval<0.01)   
diff_test_res <-diff_test_res %>% arrange(qval)
write.csv(diff_test_res, file.path(RES_DIR, "EHT_DEG_PT.csv"))   ### Supplementary Data 3 


## Gene set scores in UMAP = Fig 3f, Sup Fig 9b
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  return(cds)
}
AEC_genes <- c("8430408G22Rik", "Clu", "Crip1", "Fbln2", "Gja4", "Hey1", "Mecom", "Sat1", "Sema3g", "Sox17", "Tm4sf1", "Tsc22d1")  ### AEC genes from Kalucka et al 2020
HE_sig_genes <- c("Neurl3", "Phgdh", "Sfrp2", "Nupr1", "Mycn", "Hlf", "Gfi1", "Gck", "Ift57", "Eya2", "Lmo1")   ### HE sig genes from Hou et al 2020
preHSC_sig_genes <- c("Angpt1", "Art4", "Calcrl", "Chad", "Clu", "Crebrf", "Ctnnal1", "Cyp26b1", "Dnmt3b", "Fam84b", "Fgfr3", "Gcnt2", "Gfi1", "Hif3a", "Hlf", "Mamdc2", "Meis1", "Mllt3", "Mpl", "Msi2", "Mycn", "Nkx2-3", "Ocln", "Pdzk1ip1", "Procr", "Rab38", "Serpinf1", "Setbp1", "Slc18a2", "Trim47", "Trim6", "Vdr", "Vipr2", "Vwf")  ### pre-HSC and BM HSC common signature genes from Zhou et al 2016
AGM_HSC_sig_genes <- c("Hes1", "Hey1", "Notch4", "Sox17", "Il6st", "Nos3", "Flt1", "Kdr", "Ly6a", "Cdkn1c", "Pdzk1ip1", "Procr", "Ramp2", "Trim47", "Vwf", "Meis2", "Cdh5", "Gfi1", "Pbx1", "Ptprb", "Bmp2k", "Nrp1", "Ptpru", "Dll4")  #### AGM HSC-specific genes from Vink et al 2020
HSC_sig_genes <- c("Procr", "Pdzk1ip1", "Ltb", "Mllt3", "Ifitm1", "Gimap1", "Gimap6", "Limd2", "Trim47", "Neil2", "Vwf", "Pde1b", "Neo1", "Sqrdl", "Sult1a1", "Cd82", "Ramp2", "Ubl3", "Ly6a", "Cdkn1c", "Fgfr3", "Cldn10", "Ptpn14", "Mettl7a1", "Smtnl1", "Ctsf", "Gstm1", "Sox18", "Fads3")  ### HSC sig genes from Wilson et al 2015
go1_genes <- c("Cyr61", "Dab2", "Cd44", "Thbs1", "Cxcl12", "Gata3", "Icam1", "Itga3", "Vcam1", "Fbln2", "Selp", "Fn1", "Bcl10", "Csf1", "Col8a1") #positive cell adh
go3_genes <- c("Cyr61", "Dab2", "Thbs1", "Cxcl12", "Icam1", "Itga3", "Vcam1", "Fn1", "Ctgf", "Ptprz1") #integrin binding
agm <- estimate_score(agm, markers = AEC_genes) ### repeat for each gene set score = Fig 3f, Sup Fig 9b
plot_cell_clusters(agm,
           color_by = 'score',
           cell_size = 0.8,
           show_group_id = F) + simple_theme +
    scale_color_gradient2(low="gray80", mid="gray80", high="red3", midpoint = (((max(pData(agm)$score)-min(pData(agm)$score))/2)+min(pData(agm)$score)), space = "Lab")


## Re-Classify cells as HSC Precursor and Progenitor cell types
Gfi1_id <- row.names(subset(fData(agm), gene_short_name == "Gfi1"))
Pbx1_id <- row.names(subset(fData(agm), gene_short_name == "Pbx1"))
Cdkn1c_id <- row.names(subset(fData(agm), gene_short_name == "Cdkn1c"))
Mycn_id <- row.names(subset(fData(agm), gene_short_name == "Mycn"))
Pdzk1ip1_id <- row.names(subset(fData(agm), gene_short_name == "Pdzk1ip1"))
Mllt3_id <- row.names(subset(fData(agm), gene_short_name == "Mllt3"))
Procr_id <- row.names(subset(fData(agm), gene_short_name == "Procr"))
Dll4_id <- row.names(subset(fData(agm), gene_short_name == "Dll4"))
Vwf_id <- row.names(subset(fData(agm), gene_short_name == "Vwf"))
Flt3_id <- row.names(subset(fData(agm), gene_short_name == "Flt3"))
Il7r_id <- row.names(subset(fData(agm), gene_short_name == "Il7r"))
Csf3r_id <- row.names(subset(fData(agm), gene_short_name == "Csf3r"))
Fcgr3_id <- row.names(subset(fData(agm), gene_short_name == "Fcgr3"))
Ptprc_id <- row.names(subset(fData(agm), gene_short_name == "Ptprc"))
Itga2b_id <- row.names(subset(fData(agm), gene_short_name == "Itga2b"))
Spn_id <- row.names(subset(fData(agm), gene_short_name == "Spn"))
Runx1_id <- row.names(subset(fData(agm), gene_short_name == "Runx1"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "EC2Art", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Dll4_id,] > 0 & x[Efnb2_id,] > 0 & x[Hey1_id,] > 0 & x[Nr2f2_id,] ==0 & x[Nrp2_id,] == 0 & x[Runx1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "EC1Ven", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Nr2f2_id,] > 0 & x[Nrp2_id,] > 0 & x[Meox1_id,] == 0 & x[Runx1_id,] == 0 & x[Efnb2_id,] == 0 & x[Hey1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "HSC", classify_func=function(x) {x[Dll4_id,] >0 & x[Pdzk1ip1_id,] >0 & x[Procr_id,] > 0 & x[Vwf_id,] > 0 & x[Cdkn1c_id,] > 0 & x[Mycn_id,] > 0 & x[Pbx1_id,] > 0 & x[Mllt3_id,] > 0 & x[Gfi1_id,] > 0 & x[Flt3_id,] == 0 & x[Il7r_id,] == 0 & x[Fcgr3_id,] == 0 & x[Csf3r_id,] == 0})
cth <- addCellType(cth, "Prog", classify_func=function(x) {x[Runx1_id,] >0 & (x[Flt3_id,] > 1 | x[Il7r_id,] > 1 | x[Fcgr3_id,] > 1 | x[Csf3r_id,] > 1)})
agm <- classifyCells(agm, cth, method="markers-only")

cell_type_color <- c("1_EC1Ven" = tgray,
                     "2_EC2Art" = "royalblue3",
                    "3_HSC" = tgray,
                    "4_Prog" = tgray, 
                    "5_Unknown" = tgray)
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme  ### Fig 3g AEC

cell_type_color <- c("1_EC1Ven" = tgray,
                     "2_EC2Art" = tgray,
                     "3_HSC" = "red3",
                    "4_Prog" = tgray, 
                    "5_Unknown" = tgray)
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme  ### Fig 3h HSC precursors

### clasify cells with HSC precursor type by E10 vs E11 stage
pData(agm)$sample_cluster <- as.character(paste(pData(agm)$sample, pData(agm)$Cluster, pData(agm)$CellType))
pData(agm)$typeE10 = ((pData(agm)$sample_cluster=="E10 1 HSC")|(pData(agm)$sample_cluster=="E10 2 HSC")|(pData(agm)$sample_cluster=="E10 5 HSC"))
pData(agm)$typeE11 = ((pData(agm)$sample_cluster=="E11 1 HSC")|(pData(agm)$sample_cluster=="E11 2 HSC")|(pData(agm)$sample_cluster=="E11 5 HSC"))
pData(agm)$HSC_type = (pData(agm)$typeE10==TRUE | pData(agm)$typeE11==TRUE)
pData(agm)$HSC_class <- as.character(paste(pData(agm)$sample, pData(agm)$HSC_type))
table(pData(agm)$HSC_class)

cell_type_color <- c("1_E10 FALSE" = tgray,
                    "2_E10 TRUE" = "red3",
                    "3_E11 FALSE" = tgray,
                    "4_E11 TRUE" = tgray)
plot_cell_clusters(agm,
                   color_by = 'HSC_class',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme   ### Fig 3h E10 subset
cell_type_color <- c("1_E10 FALSE" = tgray,
                    "2_E10 TRUE" = tgray,
                    "3_E11 FALSE" = tgray,
                    "4_E11 TRUE" = "red3")
plot_cell_clusters(agm,
                   color_by = 'HSC_class',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme   ### Fig 3h E11 subset

## Re-Classify cells as HE subtype of HSC precursors based on absence ot hematopoietic marker gene = Ptprc, Spn, and Itga2b
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "HE", classify_func=function(x) {x[Dll4_id,] >0 & x[Pdzk1ip1_id,] >0 & x[Procr_id,] > 0 & x[Vwf_id,] > 0 & x[Cdkn1c_id,] > 0 & x[Mycn_id,] > 0 & x[Pbx1_id,] > 0 & x[Mllt3_id,] > 0 & x[Gfi1_id,] > 0 & x[Flt3_id,] == 0 & x[Il7r_id,] == 0 & x[Fcgr3_id,] == 0 & x[Csf3r_id,] == 0 & x[Ptprc_id,] == 0 & x[Itga2b_id,] == 0 & x[Spn_id,] == 0})
cth <- addCellType(cth, "HSC", classify_func=function(x) {x[Dll4_id,] >0 & x[Pdzk1ip1_id,] >0 & x[Procr_id,] > 0 & x[Vwf_id,] > 0 & x[Cdkn1c_id,] > 0 & x[Mycn_id,] > 0 & x[Pbx1_id,] > 0 & x[Mllt3_id,] > 0 & x[Gfi1_id,] > 0 & x[Flt3_id,] == 0 & x[Il7r_id,] == 0 & x[Fcgr3_id,] == 0 & x[Csf3r_id,] == 0 & (x[Ptprc_id,] > 0 | x[Itga2b_id,] > 0 | x[Spn_id,] > 0)})
agm <- classifyCells(agm, cth, method="markers-only")
cell_type_color <- c("1_HE" = "red3",
					 "2_HSC" = tgray,
                    "3_Unknown" = tgray)
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme  ### Fig 3h HE subset


## Re-Classify cells to plot progenitor cell type
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "EC2Art", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Dll4_id,] > 0 & x[Efnb2_id,] > 0 & x[Hey1_id,] > 0 & x[Nr2f2_id,] ==0 & x[Nrp2_id,] == 0 & x[Runx1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "EC1Ven", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Nr2f2_id,] > 0 & x[Nrp2_id,] > 0 & x[Meox1_id,] == 0 & x[Runx1_id,] == 0 & x[Efnb2_id,] == 0 & x[Hey1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "HSC", classify_func=function(x) {x[Dll4_id,] >0 & x[Pdzk1ip1_id,] >0 & x[Procr_id,] > 0 & x[Vwf_id,] > 0 & x[Cdkn1c_id,] > 0 & x[Mycn_id,] > 0 & x[Pbx1_id,] > 0 & x[Mllt3_id,] > 0 & x[Gfi1_id,] > 0 & x[Flt3_id,] == 0 & x[Il7r_id,] == 0 & x[Fcgr3_id,] == 0 & x[Csf3r_id,] == 0})
cth <- addCellType(cth, "Prog", classify_func=function(x) {x[Runx1_id,] >0 & (x[Flt3_id,] > 1 | x[Il7r_id,] > 1 | x[Fcgr3_id,] > 1 | x[Csf3r_id,] > 1)})
agm <- classifyCells(agm, cth, method="markers-only")

cell_type_color <- c("1_EC1Ven" = tgray,
                     "2_EC2Art" = tgray,
                     "3_HSC" = tgray,
                    "4_Prog" = "orange3", 
                    "5_Unknown" = tgray)
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) #+simple_theme  ### Fig 3h Progenitors


###Gene set scores by cell type - violin plots = Fig 3j, repeat for each gene set score in Fig 3f/3j
agmz <- agm[,pData(agm)$CellType != "Unknown"]
df <- data.frame(pData(agmz)$CellType, pData(agmz)$score)
names(df) <- c("CellType", "Score")
df$CellType <- factor(df$CellType, levels = c("EC1Ven", "EC2Art", "HSC", "Prog"))
cell_type_color <- c("EC1Ven" = "skyblue",
                     "EC2Art" = "royalblue3",
                    "HSC" = "red3",
                    "Prog" = "orange3")
p <- ggplot(df, aes(x= CellType, y=Score, fill = CellType)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.color=NA)       
p <- p + theme_bw()
p <- p + scale_fill_manual(values = cell_type_color)   
p <- p + xlab(NULL) + ylab(NULL) 
p <- p + theme(plot.margin = unit(c(0,3.5,0,3.5), "cm")) + ylim(0,NA)
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank()) 
p

### heatmap of Notch genes by pseudotime = Fig 3k
markers <- c("Hey1","Hey2","Nrarp","Notch1","Gata2","Hes1","Cdca7","Notch2","Csf3r")
marker_genes <- row.names(subset(fData(agmx),
                                 gene_short_name %in% markers))
diff_test_res <- differentialGeneTest(agmx[marker_genes,],
                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
pal <- colorRampPalette(c("navy", "white", "red"))
plot_pseudotime_heatmap(agmx[sig_gene_names,], hmcols = pal(500), show_rownames = TRUE) ### Fig 3k


### CellType classification for cluster 6 = Sup Fig 5g
Gata1_id <- row.names(subset(fData(agm), gene_short_name == "Gata1"))
Klf1_id <- row.names(subset(fData(agm), gene_short_name == "Klf1"))
Gypa_id <- row.names(subset(fData(agm), gene_short_name == "Gypa"))
Pf4_id <- row.names(subset(fData(agm), gene_short_name == "Pf4"))
C1qb_id <- row.names(subset(fData(agm), gene_short_name == "C1qb"))
Mrc1_id <- row.names(subset(fData(agm), gene_short_name == "Mrc1"))   
Fcnb_id <- row.names(subset(fData(agm), gene_short_name == "Fcnb"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "Mac", classify_func=function(x) {x[C1qb_id,] > 0 & x[Mrc1_id,] > 0})
cth <- addCellType(cth, "Meg", classify_func=function(x) {x[Gata1_id,] > 1 & x[Pf4_id,]>1})
cth <- addCellType(cth, "Ery", classify_func=function(x) {x[Gypa_id,] > 1 & x[Klf1_id,] > 1})
cth <- addCellType(cth, "Gran", classify_func=function(x) {x[Fcnb_id,] > 0})
agm <- classifyCells(agm, cth, method="markers-only")
cell_type_color <- c("1_Ery" = "coral4",   
                     "2_Gran" = "blueviolet",
                    "3_Mac" = "darkblue",
                    "4_Meg" = "coral",
                    "5_Unknown" = tgray)
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 0.8,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) #+simple_theme  ### Sup Fig 5g



#### ANALYSIS TO Identify candidate receptors in HSC precursors (Supplementary Data 5, Tab E)
RES_DIR2 <- file.path("~/analysis/AGMLigRec") ### folder where receptor-ligand pairs database csv file is located

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "EC2Art", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Dll4_id,] > 0 & x[Efnb2_id,] > 0 & x[Hey1_id,] > 0 & x[Nr2f2_id,] ==0 & x[Nrp2_id,] == 0 & x[Runx1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "EC1Ven", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Nr2f2_id,] > 0 & x[Nrp2_id,] > 0 & x[Meox1_id,] == 0 & x[Runx1_id,] == 0 & x[Efnb2_id,] == 0 & x[Hey1_id,] == 0 & x[Gfi1_id,] == 0})
cth <- addCellType(cth, "HSC", classify_func=function(x) {x[Dll4_id,] >0 & x[Pdzk1ip1_id,] >0 & x[Procr_id,] > 0 & x[Vwf_id,] > 0 & x[Cdkn1c_id,] > 0 & x[Mycn_id,] > 0 & x[Pbx1_id,] > 0 & x[Mllt3_id,] > 0 & x[Gfi1_id,] > 0 & x[Flt3_id,] == 0 & x[Il7r_id,] == 0 & x[Fcgr3_id,] == 0 & x[Csf3r_id,] == 0})
cth <- addCellType(cth, "Prog", classify_func=function(x) {x[Runx1_id,] >0 & (x[Flt3_id,] > 1 | x[Il7r_id,] > 1 | x[Fcgr3_id,] > 1 | x[Csf3r_id,] > 1)})
agm <- classifyCells(agm, cth, method="markers-only")

#### Identify candidate receptors in HSC precursors (Supplementary Data 5, Tab E)
agm1 <- agm[,pData(agm)$CellType == "HSC"]  
agm1 <- detectGenes(agm1, min_expr = 0.1)                                                
fData(agm1)$use_for_rec_list <- fData(agm1)$num_cells_expressed > 0.05 * ncol(agm1)
agm1_genes_common <- row.names(subset(fData(agm1), fData(agm1)$use_for_rec_list == TRUE))
gene.list = agm1_genes_common
agm1_genes_common <- getGenes(gene.list, fields='symbol')
RecLig = read.csv(file.path(RES_DIR2, "RecLigPairs.csv"))  # read csv file of Receptor Ligand pairs
AgmHSC_receptors <- intersect(RecLig$Receptor, agm1_genes_common$symbol)       
write.csv(AgmHSC_receptors, file.path(RES_DIR2, "AgmHSC_receptors.csv"))          # save csv file AGM HSC precursors receptors

#### Identify candidate ligands in primary AGM arterial EC (Supplementary Data 5, Tab D)
agm2 <- agm[,pData(agm)$CellType == "EC2Art"]  
agm2 <- detectGenes(agm2, min_expr = 0.1)                                                
fData(agm2)$use_for_rec_list <- fData(agm2)$num_cells_expressed > 0.1 * ncol(agm2)
agm2_genes_common <- row.names(subset(fData(agm2), fData(agm2)$use_for_rec_list == TRUE))
gene.list = agm2_genes_common
agm2_genes_common <- getGenes(gene.list, fields='symbol')
RecLig = read.csv(file.path(RES_DIR2, "RecLigPairs.csv"))  # read csv file RecLig of Receptor Ligand pairs
AgmArtEC_ligands <- intersect(RecLig$Ligand, agm2_genes_common$symbol)       
write.csv(AgmArtEC_ligands, file.path(RES_DIR2, "AgmArtEC_ligands.csv"))          # save csv file AGM Arterial EC ligands













### scRNAseq analysis of HSC colony derived following AGM-EC co-culture (Fig 4, Sup Fig 8)
###LOAD UNPROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample9.RDS")) 

agm <- updateCDS(agm)     
agm <- estimateSizeFactors(agm)
agm <- estimateDispersions(agm)
disp_table = dispersionTable(agm)
disp_table = disp_table %>% mutate(excess_disp =
  (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 1000)$gene_id)
agm = setOrderingFilter(agm, top_subset_genes)
agm <- preprocessCDS(agm, num_dim = 7)                    
agm <- reduceDimension(agm, reduction_method = 'UMAP')   
agm <- clusterCells(agm,
                        method = 'louvain',
                        res = 1e-4,   
                        verbose = T)
### Redefine clusters as EC or Hematopoietic
pData(agm)$Cluster <- plyr::revalue(as.character(pData(agm)$Cluster),
                                        c("1" = 'EC',
                                        "2" = 'Hem',
                                        "3" = 'EC'))

###LOAD PROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample9_P.RDS"))

### Plot in UMAP by Ptprc expression = Fig 4c
plot_cell_clusters(agm,                 
                   markers = c('Ptprc'),                          
                   color_by=marker, 
                   show_group_id = T, cell_size = 1) +
                   scale_color_viridis(option = "magma") +simple_theme

### Plot in UMAP by cluster (EC vs Hem) = Sup Fig 8a
cluster_color <- c("EC" = "blue", "Hem" = "red")
plot_cell_clusters(agm,
                   color_by = 'Cluster',
                   cell_size = 1,
                   show_group_id = F) +scale_color_manual(values = cluster_color) +simple_theme

### UMI and genes plotted = Sup Fig 8c
agm <- detectGenes(agm, min_expr=0.1)
pData(agm)$UMI <- Matrix::colSums(exprs(agm))
summary(pData(agm)$UMI)
summary(pData(agm)$num_genes_expressed)
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
        panel.grid = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size=15))
sample_color <- c("EC" = "blue", "Hem" = "red")
ggplot(pData(agm), aes(x=Cluster, y=UMI, fill=Cluster)) +geom_boxplot() +  scale_y_continuous(trans='log10', limits=c(100,100000), labels=number) + simple_theme2 +
   scale_fill_manual(values = sample_color) 
ggplot(pData(agm), aes(x=Cluster, y=num_genes_expressed, fill=Cluster)) +geom_boxplot()  +
                  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 +   scale_fill_manual(values = sample_color) 


## Classify celltypes as HSC and HPC based on signature genes
Runx1_id <- row.names(subset(fData(agm), gene_short_name == "Runx1"))
Pdzk1ip1_id <- row.names(subset(fData(agm), gene_short_name == "Pdzk1ip1"))
Vwf_id <- row.names(subset(fData(agm), gene_short_name == "Vwf"))
Cd48_id <- row.names(subset(fData(agm), gene_short_name == "Cd48"))
Itgal_id <- row.names(subset(fData(agm), gene_short_name == "Itgal"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "HPC", classify_func=function(x) {x[Itgal_id,] > 0 | x[Cd48_id,] > 0})
cth <- addCellType(cth, "HSC", classify_func=function(x) {x[Runx1_id,] >0 & x[Pdzk1ip1_id,] >0 & x[Vwf_id,] >0 & x[Cd48_id,] == 0 & x[Itgal_id,] == 0})
agm <- classifyCells(agm, cth, method="markers-only")
tgray=alpha("gray", 0.1)
cell_type_color <- c("1_HPC" = "navy",
                    "2_HSC" = "red3",
                    "3_Unknown" = tgray)

### plot UMAP by CellType = Sup Fig 8f
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme

### Plot gene x in UMAP = Fig 4e (genes = Runx1, Pdzk1ip1, Vwf,Fgd5, Procr, Itgal, Cd48, Cdkn1c, Mki67)
plot_cell_clusters(agm,                 
                   markers = c('x'),                          
                   color_by=marker, 
                   show_group_id = T, cell_size = 1) +
                   scale_color_viridis(option = "magma") +simple_theme

### Compare HSC gene set scores (Wilson et al 2015) between HSC and HPC CellTypes = Sup Fig 8h (continues below for combined HSC & HPC from both colony types)
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  pData(cds)$index = log(aggregate_score + 1)
  return(cds)
}
HSC_sig_genes <- c("Procr", "Pdzk1ip1", "Ltb", "Mllt3", "Ifitm1", "Gimap1", "Gimap6", "Limd2", "Trim47", "Neil2", "Vwf", "Pde1b", "Neo1", "Sqrdl", "Sult1a1", "Cd82", "Ramp2", "Ubl3", "Ly6a", "Cdkn1c", "Fgfr3", "Cldn10", "Ptpn14", "Mettl7a1", "Smtnl1", "Ctsf", "Gstm1", "Sox18", "Fads3")  ### HSC sig genes from Wilson 2015
agm <- estimate_score(agm, markers = HSC_sig_genes) 
agmz <- agm[,pData(agm)$Cluster == "Hem"]
agmz <- agmz[,pData(agmz)$CellType != "Unknown"]
df1 <- data.frame(pData(agmz)$CellType, pData(agmz)$score)
names(df1) <- c("CellType", "Score")
df1$CellType <- factor(df1$CellType, levels = c("HSC", "HPC"))


###### Identify receptors expressed in HSC CellType 
agm1 <- agm[,pData(agm)$CellType == "HSC"]  
table(pData(agm1)$CellType)
agm1 <- detectGenes(agm1, min_expr = 0.1)                                               
fData(agm1)$use_for_rec_list <- fData(agm1)$num_cells_expressed > 0.2 * ncol(agm1)
agm1_genes_common <- row.names(subset(fData(agm1), fData(agm1)$use_for_rec_list == TRUE))
gene.list = agm1_genes_common
agm1_genes_common <- getGenes(gene.list, fields='symbol')
RecLig = read.csv(file.path(RES_DIR2, "RecLigPairs.csv"))  # read csv file RecLig of Receptor Ligand pairs
HSC1_receptors <- intersect(RecLig$Receptor, agm1_genes_common$symbol)        # create list of Receptors
write.csv(HSC1_receptors, file.path(RES_DIR2, "HSC1_receptors.csv"))          # save csv file receptors expressed in HSC generated in vitro in colony type 1









### scRNAseq analysis of HSC colony derived following AGM-EC co-culture (Fig 4, Sup Fig 8)
###LOAD UNPROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample10.RDS"))

agm <- updateCDS(agm)  
agm <- estimateSizeFactors(agm)
agm <- estimateDispersions(agm)
disp_table = dispersionTable(agm)
disp_table = disp_table %>% mutate(excess_disp =
  (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 1000)$gene_id)
agm = setOrderingFilter(agm, top_subset_genes)
agm <- preprocessCDS(agm, num_dim = 15)                    
agm <- reduceDimension(agm, reduction_method = 'UMAP')
agm <- clusterCells(agm,
                        method = 'louvain',
                        res = 1e-4,    ### adjust resolution of clusters
                        verbose = T)
### Redefine clusters as EC or Hematopoietic
pData(agm)$Cluster <- plyr::revalue(as.character(pData(agm)$Cluster),
                                        c("1" = 'Hem',
                                        "2" = 'EC',
                                        "3" = 'EC', "4"='EC'))

###LOAD PROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample10_P.RDS"))

### Plot in UMAP by Ptprc expression = Fig 4d
plot_cell_clusters(agm,                 
                   markers = c('Ptprc'),                          
                   color_by=marker, 
                   show_group_id = T, cell_size = 1) +
                   scale_color_viridis(option = "magma") +simple_theme

### Plot in UMAP by cluster (EC vs Hem) = Sup Fig 8b
cluster_color <- c("EC" = "blue", "Hem" = "red")
plot_cell_clusters(agm,
                   color_by = 'Cluster',
                   cell_size = 1,
                   show_group_id = F) +scale_color_manual(values = cluster_color) +simple_theme


### UMI and genes plotted = Sup Fig 8d
agm <- detectGenes(agm, min_expr=0.1)
pData(agm)$UMI <- Matrix::colSums(exprs(agm))
summary(pData(agm)$UMI)
summary(pData(agm)$num_genes_expressed)
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
        panel.grid = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size=15))
sample_color <- c("EC" = "blue", "Hem" = "red")
ggplot(pData(agm), aes(x=Cluster, y=UMI, fill=Cluster)) +geom_boxplot() +  scale_y_continuous(trans='log10', limits=c(100,100000), labels=number) + simple_theme2 +
   scale_fill_manual(values = sample_color) 
ggplot(pData(agm), aes(x=Cluster, y=num_genes_expressed, fill=Cluster)) +geom_boxplot()  +
                  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 +   scale_fill_manual(values = sample_color) 


## Classify celltypes as HSC and HPC based on signature genes
Runx1_id <- row.names(subset(fData(agm), gene_short_name == "Runx1"))
Pdzk1ip1_id <- row.names(subset(fData(agm), gene_short_name == "Pdzk1ip1"))
Vwf_id <- row.names(subset(fData(agm), gene_short_name == "Vwf"))
Cd48_id <- row.names(subset(fData(agm), gene_short_name == "Cd48"))
Itgal_id <- row.names(subset(fData(agm), gene_short_name == "Itgal"))
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "HPC", classify_func=function(x) {x[Itgal_id,] > 0 | x[Cd48_id,] > 0})
cth <- addCellType(cth, "HSC", classify_func=function(x) {x[Runx1_id,] >0 & x[Pdzk1ip1_id,] >0 & x[Vwf_id,] >0 & x[Cd48_id,] == 0 & x[Itgal_id,] == 0})
agm <- classifyCells(agm, cth, method="markers-only")

tgray=alpha("gray", 0.1)
cell_type_color <- c("1_HPC" = "navy",
                    "2_HSC" = "red3",
                    "3_Unknown" = tgray)

### plot UMAP by CellType = Fig 4g, Sup Fig 8g
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 1,
                   show_group_id = F) +
                   scale_color_manual(values = cell_type_color) +simple_theme

### Plot gene x in UMAP = Fig 4f (genes = Runx1, Pdzk1ip1, Itgal)
plot_cell_clusters(agm,                 
                   markers = c('x'),                          
                   color_by=marker, 
                   show_group_id = T, cell_size = 1) +
                   scale_color_viridis(option = "magma") +simple_theme


#### Psuedotime trajectory analysis = Fig 4h left panel
agm <- partitionCells(agm)
agm <- learnGraph(agm,  RGE_method = 'SimplePPT')
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
root_node_ids = get_correct_root_state(agm,
                                      cell_phenotype =
                                        'Cluster', "Hem")
agm <- orderCells(agm, root_pr_nodes = root_node_ids, reverse=TRUE)
plot_cell_trajectory(agm, cell_size = 0.5) +simple_theme +scale_color_viridis(option = "viridis")


### limit pseudotime analysis to Hem cells
agmx <- agm[,pData(agm)$Cluster == "Hem"]
agmx <- agmx[,pData(agmx)$Pseudotime >0]

##### Heatmap of gene expression over pseudotime for HSC and HPC marker genes = Fig 4h right panel
markers <- c("Cebpa", "Ebf1", "Cd48", "Itgal", "Hlf", "Cdkn1c", "Ly6a", "Mllt3", "Esam", "Procr", "Fgd5", "Pdzk1ip1", "Vwf", "Hes1", "Gata2", "Gata3", "Nrarp", "Mecom", "Pbx1")
marker_genes <- row.names(subset(fData(agmx),
                                 gene_short_name %in% markers))
diff_test_res <- differentialGeneTest(agmx[marker_genes,],
                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.05))
pal <- colorRampPalette(c("navy", "white", "red"))
plot_pseudotime_heatmap(agmx[sig_gene_names,], hmcols = pal(500),
                        cores = 1, show_rownames = T, num_clusters=2)


### global analysis of genes DE over pseudotime = Supplementary Data 4
diff_test_res <- differentialGeneTest(agmx,
                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res <-subset(diff_test_res, num_cells_expressed > (0.05 * ncol(agmx)))
diff_test_res <-subset(diff_test_res, qval<0.01) 
diff_test_res <-diff_test_res %>% arrange(qval)
write.csv(diff_test_res, file.path(RES_DIR, "HSCtoHPC_DEG_PT.csv"))     

### Compare HSC gene set scores (Wilson et al 2015) between HSC and HPC CellTypes
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  pData(cds)$index = log(aggregate_score + 1)
  return(cds)
}
HSC_sig_genes <- c("Procr", "Pdzk1ip1", "Ltb", "Mllt3", "Ifitm1", "Gimap1", "Gimap6", "Limd2", "Trim47", "Neil2", "Vwf", "Pde1b", "Neo1", "Sqrdl", "Sult1a1", "Cd82", "Ramp2", "Ubl3", "Ly6a", "Cdkn1c", "Fgfr3", "Cldn10", "Ptpn14", "Mettl7a1", "Smtnl1", "Ctsf", "Gstm1", "Sox18", "Fads3")  ### HSC sig genes from Wilson 2015
agm <- estimate_score(agm, markers = HSC_sig_genes) 
agmz <- agm[,pData(agm)$Cluster == "Hem"]
agmz <- agmz[,pData(agmz)$CellType != "Unknown"]
df2 <- data.frame(pData(agmz)$CellType, pData(agmz)$score)
names(df2) <- c("CellType", "Score")
df2$CellType <- factor(df2$CellType, levels = c("HSC", "HPC"))

##combine data from HSC and HPC from first colony type (above) and plot HSC gene set score (Wilson 2015) as violin plot = Sup Fig 8h
df3 <- rbind(df1, df2)
df3$CellType <- factor(df3$CellType, levels = c("HSC", "HPC"))
cell_type_color <- c("HSC" = "red3", "HPC" = "navy")
p <- ggplot(df3, aes(x= CellType, y=Score, fill = CellType)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + stat_compare_means()
p <- p+ stat_mean()
p <- p + theme_bw()
p <- p + scale_fill_manual(values = cell_type_color)  
p <- p + xlab(NULL) + ylab(NULL) ##renames axis label
p <- p + theme(plot.margin = unit(c(0,3.5,0,3.5), "cm")) + ylim(0,NA)
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p



### Identify expressed genes in HSC Cell Type
agm1 <- agm[,pData(agm)$CellType == "HSC"] 
table(pData(agm1)$CellType)
agm1 <- detectGenes(agm1, min_expr = 0.1)                                               
fData(agm1)$use_for_rec_list <- fData(agm1)$num_cells_expressed > 0.2 * ncol(agm1)
agm1_genes_common <- row.names(subset(fData(agm1), fData(agm1)$use_for_rec_list == TRUE))
gene.list = agm1_genes_common
agm1_genes_common <- getGenes(gene.list, fields='symbol')
RecLig = read.csv(file.path(RES_DIR, "RecLigPairs.csv"))  # read csv file RecLig of Receptor Ligand pairs
HSC2_receptors <- intersect(RecLig$Receptor, agm1_genes_common$symbol)        # create list of Receptors
write.csv(HSC2_receptors, file.path(RES_DIR2, "HSC2_receptors.csv"))          # save csv file receptors expressed in HSC generated in vitro in colony type 2



###  Combined Analysis of Ligand-Receptor pairs = Fig 5 and Sup Data 5
RecLig = read.csv(file.path(RES_DIR, "RecLigPairs.csv"))                    # read csv file RecLig of Receptor Ligand pairs
AGMEC_ligands = read.csv(file.path(RES_DIR2, "AgmEC_Ligands.csv"))          # load csv file of Ligands expressed by supportive AGM-EC lines
AgmArtEC_ligands = read.csv(file.path(RES_DIR2, "AgmArtEC_ligands.csv"))    # load csv file Ligands expressed by primary AGM arterial EC
EC_ligands <- union(AGMEC_ligands$x, AgmArtEC_ligands$x)                    # combine ligands from AGM-EC and primary arterial EC
AgmHSC_receptors = read.csv(file.path(RES_DIR2, "AgmHSC_receptors.csv"))    # load csv file of receptors expressed by primary AGM HSC precursors

HSC1_receptors = read.csv(file.path(RES_DIR2, "HSC1_receptors.csv"))        # load csv file of receptors expressed by HSC from colony type 1
HSC2_receptors = read.csv(file.path(RES_DIR2, "HSC2_receptors.csv"))        # load csv file of receptors expressed by HSC from colony type 2
HSC_receptors <- union(HSC1_receptors$x, HSC2_receptors$x)                  # combine receptors expressed by HSC from both colony types
write.csv(HSC_receptors, file.path(RES_DIR2, "HSC1_2_receptors.csv"))       # save csv file of receptors expressed by HSC generated in vitro = Sup Data 5 (Tab F)

###Create Receptor Ligand list with combined primary EC and AGM-EC
AgmHSC_RecLig <- RecLig[(RecLig$Ligand %in% EC_ligands$x & RecLig$Receptor %in% AgmHSC_receptors$x),]    # find Rec Ligand Pairs involving AGM HSC precursors
HSC_RecLig <- RecLig[(RecLig$Ligand %in% EC_ligands$x & RecLig$Receptor %in% HSC_receptors),]            # find Rec Ligand Pairs involving HSC generated in vitro
write.csv(AgmHSC_RecLig, file.path(RES_DIR2, "AgmHSC_RecLig.csv"))          # save csv file of receptor ligand pairs involving AGM HSC precursors = Sup Data 5 (Tab A)
write.csv(HSC_RecLig, file.path(RES_DIR2, "HSC_RecLig.csv"))                # save csv file of receptor ligand pairs involving HSC generated in vitro = Sup Data 5 (Tab B)












### scRNAseq analysis of hematopoietic populations generated following AGM-EC vs Engineered nich culture = Fig 6, Sup Fig 11
###LOAD UNPROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample11-12.RDS"))

agm <- updateCDS(agm)     
agm <- estimateSizeFactors(agm)
agm <- estimateDispersions(agm)
disp_table = dispersionTable(agm)
disp_table = disp_table %>% mutate(excess_disp =
  (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 2500)$gene_id)
agm = setOrderingFilter(agm, top_subset_genes)
agm <- preprocessCDS(agm, num_dim = 10)                    
agm <- reduceDimension(agm, reduction_method = 'UMAP')   
agm <- clusterCells(agm,
                        method = 'louvain',
                        res = 1e-3,    ### adjust resolution of clusters
                        verbose = T)


###LOAD PROCESSED CDS###
agm <- readRDS(file.path(RES_DIR, "Sample11-12_P.RDS"))

### UMI and genes plotted = Sup Fig 11d
agm <- detectGenes(agm, min_expr=0.1)
pData(agm)$UMI <- Matrix::colSums(exprs(agm))
summary(pData(agm)$UMI)
summary(pData(agm)$num_genes_expressed)
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
        panel.grid = element_blank(), axis.title = element_blank(), axis.text.y = element_text(size=15))
sample_color <- c("EC" = "red3", "Eng" = "blue4")
ggplot(pData(agm), aes(x=sample, y=UMI, fill=sample)) +geom_boxplot() +  scale_y_continuous(trans='log10', limits=c(100,100000), labels=number) + simple_theme2 +
   scale_fill_manual(values = sample_color) 
ggplot(pData(agm), aes(x=sample, y=num_genes_expressed, fill=sample)) +geom_boxplot()  +
                  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 +   scale_fill_manual(values = sample_color) 


### Plot Psuedotime Trajectory in UMAP = Fig 6e
agm <- partitionCells(agm)
agm <- learnGraph(agm,  RGE_method = 'SimplePPT')
get_correct_root_state <- function(cds, cell_phenotype, root_type){
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
root_node_ids = get_correct_root_state(agm,
                                      cell_phenotype =
                                        'Cluster', "4")
agm <- orderCells(agm, root_pr_nodes = root_node_ids)
plot_cell_trajectory(agm, cell_size = 0.01) +
    scale_color_viridis(option = "viridis") +simple_theme


### Classify cell types
Cdh5_id <- row.names(subset(fData(agm), gene_short_name == "Cdh5"))
Kdr_id <- row.names(subset(fData(agm), gene_short_name == "Kdr"))
Ptprc_id <- row.names(subset(fData(agm), gene_short_name == "Ptprc"))
Pdzk1ip1_id <- row.names(subset(fData(agm), gene_short_name == "Pdzk1ip1"))
Hlf_id <- row.names(subset(fData(agm), gene_short_name == "Hlf"))
Procr_id <- row.names(subset(fData(agm), gene_short_name == "Procr"))
Kit_id <- row.names(subset(fData(agm), gene_short_name == "Kit"))
Cd34_id <- row.names(subset(fData(agm), gene_short_name == "Cd34"))
Flt3_id <- row.names(subset(fData(agm), gene_short_name == "Flt3"))
Vwf_id <- row.names(subset(fData(agm), gene_short_name == "Vwf"))
Mpl_id <- row.names(subset(fData(agm), gene_short_name == "Mpl"))
Gypa_id <- row.names(subset(fData(agm), gene_short_name == "Gypa"))
Klf1_id <- row.names(subset(fData(agm), gene_short_name == "Klf1"))
Il7r_id <- row.names(subset(fData(agm), gene_short_name == "Il7r"))
Ptcra_id <- row.names(subset(fData(agm), gene_short_name == "Ptcra"))
Ebf1_id <- row.names(subset(fData(agm), gene_short_name == "Ebf1"))
Adgre1_id <- row.names(subset(fData(agm), gene_short_name == "Adgre1"))
C1qb_id <- row.names(subset(fData(agm), gene_short_name == "C1qb"))
Fcnb_id <- row.names(subset(fData(agm), gene_short_name == "Fcnb"))
Elane_id <- row.names(subset(fData(agm), gene_short_name == "Elane"))
Siglech_id <- row.names(subset(fData(agm), gene_short_name == "Siglech"))
Prss34_id <- row.names(subset(fData(agm), gene_short_name == "Prss34"))
Cd200r3_id <- row.names(subset(fData(agm), gene_short_name == "Cd200r3"))
Cd7_id <- row.names(subset(fData(agm), gene_short_name == "Cd7"))
Fcgr3_id <- row.names(subset(fData(agm), gene_short_name == "Fcgr3"))
Cd14_id <- row.names(subset(fData(agm), gene_short_name == "Cd14"))
Cd36_id <- row.names(subset(fData(agm), gene_short_name == "Cd36"))
Cebpa_id <- row.names(subset(fData(agm), gene_short_name == "Cebpa"))
Gata1_id <- row.names(subset(fData(agm), gene_short_name == "Gata1"))
Gfi1b_id <- row.names(subset(fData(agm), gene_short_name == "Gfi1b"))
Irf8_id <- row.names(subset(fData(agm), gene_short_name == "Irf8"))
Zeb2_id <- row.names(subset(fData(agm), gene_short_name == "Zeb2"))
Gfi1_id <- row.names(subset(fData(agm), gene_short_name == "Gfi1"))

cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "HSC", classify_func=function(x) {x[Pdzk1ip1_id,] > 0 & (x[Hlf_id,] > 0 | x[Procr_id,] > 0) & x[Kdr_id,] == 0})
cth <- addCellType(cth, "EC", classify_func=function(x) {x[Cdh5_id,] > 0 & x[Kdr_id,] > 1 & x[Ptprc_id,] == 0})
cth <- addCellType(cth, "MPP", classify_func=function(x) {x[Pdzk1ip1_id,] == 0 & x[Kit_id,] > 0 & x[Cd34_id,] > 0 & x[Flt3_id,] > 0})
cth <- addCellType(cth, "GMP", classify_func=function(x) {(x[Hlf_id,] == 0 & x[Cd34_id,] > 0 & x[Fcgr3_id,] > 0 & x[Cebpa_id,] > 0) | (x[Irf8_id,] > 0 & x[Zeb2_id,] > 0 & x[C1qb_id,] == 0 & x[Fcnb_id,] == 0 & x[Cd14_id,] == 0 & x[Gfi1_id,] == 0 )})
cth <- addCellType(cth, "MEP", classify_func=function(x) {x[Gypa_id,] == 0 & x[Cd34_id,] == 0 & x[Gata1_id,] > 0 & x[Vwf_id,] == 0 & x[Prss34_id,] == 0 & x[Cd200r3_id,] == 0 & x[Gfi1b_id,] > 0})
cth <- addCellType(cth, "Meg", classify_func=function(x) {x[Vwf_id,] > 1 & x[Mpl_id,] > 0 & x[Hlf_id,] ==0})
cth <- addCellType(cth, "Ery", classify_func=function(x) {x[Gypa_id,] > 1 & x[Klf1_id,] > 1})
cth <- addCellType(cth, "T", classify_func=function(x) {x[Ptcra_id,] >1 & x[Il7r_id,] > 0})
cth <- addCellType(cth, "B", classify_func=function(x) {x[Ebf1_id,] >0 & x[Hlf_id,] ==0 & x[Pdzk1ip1_id,] ==0 & x[Il7r_id,] > 0})
cth <- addCellType(cth, "Mac", classify_func=function(x) {x[C1qb_id,] > 0 & x[Cd36_id,] >0 })
cth <- addCellType(cth, "DC", classify_func=function(x) {x[Siglech_id,] > 0 & x[Cd7_id,] > 0 })
cth <- addCellType(cth, "Neu", classify_func=function(x) {x[Fcnb_id,] > 1 & x[Elane_id,] > 1})
cth <- addCellType(cth, "Mast", classify_func=function(x) {x[Prss34_id,] > 0 & x[Cd200r3_id,] > 0})
cth <- addCellType(cth, "Mono", classify_func=function(x) {x[Cd14_id,] > 0 & x[Adgre1_id,] > 0 & x[C1qb_id,] == 0 & x[Fcnb_id,] == 0 })
agm <- classifyCells(agm, cth, method="markers-only")

### Plot UMAP of Cell Types = Fig 6f
tgray=alpha("gray", 0.05)
celltype_color <- c("1_Ambiguous" = tgray, "2_B" = "royalblue1", "3_DC" = "turquoise",
                    "4_EC" = "azure4", "5_Ery" = "coral4", "6_GMP" = "olivedrab3", "7_HSC" = "orangered3", 
                    "8_Mac" = "khaki4", "9_Mast" = "yellow", "10_Meg" = "coral", "11_MEP" = "gold2", "12_Mono" = "khaki4",
                    "13_MPP" = "pink", "14_Neu" = "darkgreen", "15_T" = "navy", "16_Unknown" = tgray)
plot_cell_clusters(agm,
                   color_by = 'CellType',
                   cell_size = 0.8,
                   show_group_id = F)  + scale_color_manual(values = celltype_color) +simple_theme

### Plot UMAP by sample (EC = AGM-EC co-cultured vs Eng = Engineered conditions) = Fig 6g
tgray=alpha("gray", 0)
sample_color <- c("1_EC" = "red3", "2_Eng" = tgray)
sample_color2 <- c("1_EC" = tgray, "2_Eng" = "blue4")  
plot_cell_clusters(agm,
                   color_by = 'sample',
                   cell_size = 0.8,
                   show_group_id = F) +
                   scale_color_manual(values = sample_color) +simple_theme  
plot_cell_clusters(agm,
                   color_by = 'sample',
                   cell_size = 0.8,
                   show_group_id = F) +
                   scale_color_manual(values = sample_color2) +simple_theme

### create new field combined sample and celltype = sample_type
pData(agm)$sample_type <- as.character(paste(pData(agm)$sample, pData(agm)$CellType))

### Plot UMAP of Cell Types for AGM-EC co-culture sample only = Fig 6h (left panel)
tgray2=alpha("gray", 0.05)
tgray=alpha("gray", 0)
celltype_color2 <- c("1_EC Ambiguous" = tgray2, "2_EC B" = "royalblue1", "3_EC DC" = "turquoise",
                    "4_EC EC" = "azure4", "5_EC Ery" = "coral4", "6_EC GMP" = "olivedrab3", "7_EC HSC" = "orangered3", 
                    "8_EC Mac" = "khaki4", "9_EC Mast" = "yellow", "10_EC Meg" = "coral", "11_EC MEP" = "gold2", "12_EC Mono" = "khaki4",
                    "13_EC MPP" = "pink", "14_EC Neu" = "darkgreen", "15_EC Unknown" = tgray2,
                    "16_Eng Ambiguous" = tgray, "17_Eng B" = tgray,
                    "18_Eng Ery" = tgray, "19_Eng GMP" = tgray, "20_Eng HSC" = tgray, 
                    "21_Eng Mac" = tgray, "22_Eng Mast" = tgray, "23_Eng Meg" = tgray, "24_Eng MEP" = tgray, "25_Eng Mono" = tgray,
                    "26_Eng MPP" = tgray, "27_Eng Neu" = tgray, "28_Eng T" = tgray, "29_Eng Unknown" = tgray)
plot_cell_clusters(agm,
                   color_by = 'sample_type',
                   cell_size = 0.8,
                   show_group_id = F)  + scale_color_manual(values = celltype_color2) +simple_theme

### Plot UMAP of Cell Types for Engingeered culture sample only = Fig 6h (right panel)
celltype_color3 <- c("1_EC Ambiguous" = tgray, "2_EC B" = tgray, "3_EC DC" = tgray,
                    "4_EC EC" = tgray, "5_EC Ery" = tgray, "6_EC GMP" = tgray, "7_EC HSC" = tgray, 
                    "8_EC Mac" = tgray, "9_EC Mast" = tgray, "10_EC Meg" = tgray, "11_EC MEP" = tgray, "12_EC Mono" = tgray,
                    "13_EC MPP" = tgray, "14_EC Neu" = tgray, "15_EC Unknown" = tgray,
                    "16_Eng Ambiguous" = tgray2, "17_Eng B" = "royalblue1",
                    "18_Eng Ery" = "coral4", "19_Eng GMP" = "olivedrab3", "20_Eng HSC" = "orangered3", 
                    "21_Eng Mac" = "khaki4", "22_Eng Mast" = "yellow", "23_Eng Meg" = "coral", "24_Eng MEP" = "gold2", "25_Eng Mono" = "khaki4",
                    "26_Eng MPP" = "pink", "27_Eng Neu" = "darkgreen", "28_Eng T" = "navy", "29_Eng Unknown" = tgray2)
plot_cell_clusters(agm,
                   color_by = 'sample_type',
                   cell_size = 0.8,
                   show_group_id = F)  + scale_color_manual(values = celltype_color3) +simple_theme












