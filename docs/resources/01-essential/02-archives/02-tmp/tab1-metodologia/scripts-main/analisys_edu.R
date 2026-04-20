source("~/Bioinformatica/bin/scRNASeq.R")
#set working directoty
setwd("/mnt/beegfs/Bioinformatica/INT/05-INT-24-scRNASeq_Timo_Mdelgado/Analisis/")
 #Leemos los ficheros rds que se crean despues de cellranger
readscRDS(FullPath="/mnt/beegfs/Bioinformatica/INT/05-INT-24-scRNASeq_Timo_Mdelgado/Resultados/")
#By default, merge() will combine the Seurat objects based on the raw count matrices,erasing any previously normalized and scaled data matrices.
mergeSeurat<-mergescRDS(pattern = "SCS003_24", project = "MarioD")
table(mergeSeurat@meta.data$orig.ident)

#Check initial quality metrics (mt, n_features, highly expressed genes ...)
mergeSeurat<-quality(mergeSeurat)
#Filter
mergeSeurat.filter<-filter(mergeSeurat, nFeature_thresh = c(200,6500), mt_thresh = 15, ndims = 12)

#Cluster markers in cluster
find_Markers(mergeSeurat.filter,fc=1, number=2)

#Cell annottation
mergeSeurat.filter<-cell_anottation(mergeSeurat.filter)

#Cluster markers in cell anottation
find_Markers(mergeSeurat.filter, idents = "GeneralCell_HPCA")
find_Markers(mergeSeurat.filter, idents = "GeneralCell_Mouse")
find_Markers(mergeSeurat.filter, idents = "GeneralCell_Monaco")

#Tabula Muris
mergeSeurat.filter <- tabula_muris(mergeSeurat.filter)
find_Markers(mergeSeurat.filter, idents = "cell_anot")

#DEG####
DEG_by_column(mergeSeurat.filter,column="orig.ident",label="WT_vs_KO")

DEG_by_column(mergeSeurat.filter,column="cell_anot",label="Tabula_muris")
DEG_by_column_anot(mergeSeurat.filter,column="orig.ident",anot="cell_anot",label="TabulaMuris")

DEG_by_column(mergeSeurat.filter,column="GeneralCell_Monaco",label="Monaco")
DEG_by_column_anot(mergeSeurat.filter,column="orig.ident",anot="GeneralCell_Monaco",label="Monaco")

DEG_by_column(mergeSeurat.filter,column="GeneralCell_HPCA",label="HPCA")
DEG_by_column_anot(mergeSeurat.filter,column="orig.ident",anot="GeneralCell_HPCA",label="HPCA")

DEG_by_column(mergeSeurat.filter,column="GeneralCell_Mouse",label="RNAseqMouse")
DEG_by_column_anot(mergeSeurat.filter,column="orig.ident",anot="GeneralCell_Mouse",label="RNASeqMouse")


#Cort####
pdf(file = "Cort.pdf", paper = "a4r")
expr_values <- FetchData(mergeSeurat.filter, vars = "Cort")
celdas_expresoras <- rownames(expr_values)[expr_values$Cort > 0]

FeaturePlot(mergeSeurat.filter, features = "Cort", cells = celdas_expresoras, cols = c("lightgrey", "red"))

DotPlot(mergeSeurat, features = "Cort") + scale_size(range = c(0.5, 6))

FeaturePlot(mergeSeurat.filter, features = "Cort", cols = c("lightgrey", "red"), pt.size = 1.5, min.cutoff = "q5")

VlnPlot(mergeSeurat,features = "Cort")

dev.off()
#Other####
pdf(file = "More_plots.pdf", width = 11.69, height = 8.27)

library(dittoSeq)
dittoBarPlot(mergeSeurat.filter, 
             var = "cell_anot", 
             group.by = "orig.ident" ,scale = "count")
dittoBarPlot(mergeSeurat.filter, 
             var = "GeneralCell_Monaco", 
             group.by = "orig.ident" ,scale = "count") 
dittoBarPlot(mergeSeurat.filter, 
             var = "GeneralCell_HPCA", 
             group.by = "orig.ident" ,scale = "count") 
dittoBarPlot(mergeSeurat.filter, 
             var = "GeneralCell_Mouse", 
             group.by = "orig.ident" ,scale = "count") 
dittoBarPlot(mergeSeurat.filter, 
             var = "cell_anot", 
             group.by = "orig.ident" )
dittoBarPlot(mergeSeurat.filter, 
             var = "GeneralCell_Monaco", 
             group.by = "orig.ident" ) 
dittoBarPlot(mergeSeurat.filter, 
             var = "GeneralCell_HPCA", 
             group.by = "orig.ident" ) 
dittoBarPlot(mergeSeurat.filter, 
             var = "GeneralCell_Mouse", 
             group.by = "orig.ident" ) 

#+ scale_fill_manual(values = metadata(mergeSeurat.filter)$color_vectors$celltype)
#+ 

dittoFreqPlot(mergeSeurat.filter, "cell_anot",sample.by = "orig.ident", group.by = "orig.ident")
dittoFreqPlot(mergeSeurat.filter, "GeneralCell_Mouse",sample.by = "orig.ident", group.by = "orig.ident")
dittoFreqPlot(mergeSeurat.filter, "GeneralCell_HPCA",sample.by = "orig.ident", group.by = "orig.ident")
dittoFreqPlot(mergeSeurat.filter, "GeneralCell_Monaco",sample.by = "orig.ident", group.by = "orig.ident")

dev.off()


library(viridis)
dittoHeatmap(mergeSeurat.filter, c(u,d),
             annot.by = c("cell_anot", "orig.ident"),
             scaled.to.max =  TRUE,
             #cluster_cols = TRUE,
             heatmap.colors = viridis(100),
             cluster_cols = TRUE, 
             scale = "none")






markers<- FindMarkers(mergeSeurat.filter, ident.1 = 0)
giveMeMarkers<-function(file,number=3,ident=""){
  x<-read.delim(file, header = T, row.names = 1)
  x<-x[x$p_val_adj<=0.05 & x$avg_log2FC>0,]
  return(rownames(head(x[order(x$avg_log2FC,decreasing = T),],number)))
}




# Visualizar las anotaciones transferidas
DimPlot(mergeSeurat.filter.tabula, group.by = "predicted.id", label = TRUE, repel = TRUE, pt.size = 0.5) + 
  ggtitle("Annotated Thymus Cells - Experimental Data")

DimPlot(mergeSeurat.filter,label = T, group.by = "seurat_clusters", repel = T)

DimPlot(mergeSeurat.filter,label = T, group.by = "GeneralCell_HPCA", repel = T)
DimPlot(mergeSeurat.filter,label = T, group.by = "GeneralCell_HPCA", repel = T, split.by = "orig.ident")

DimPlot(mergeSeurat.filter,label = T, group.by = "PreciseCell_HPCA", repel = T)
DimPlot(mergeSeurat.filter,label = T, group.by = "PreciseCell_HPCA", repel = T, split.by = "orig.ident")

DimPlot(mergeSeurat.filter,label = T, group.by = "GeneralCell_Monaco", repel = T)
DimPlot(mergeSeurat.filter,label = T, group.by = "GeneralCell_Monaco", repel = T, split.by = "orig.ident")

DimPlot(mergeSeurat.filter,label = T, group.by = "PreciseCell_Monaco", repel = T)
DimPlot(mergeSeurat.filter,label = T, group.by = "PreciseCell_Monaco", repel = T, split.by = "orig.ident")

DimPlot(mergeSeurat.filter,label = T, group.by = "GeneralCell_Mouse", repel = T)
DimPlot(mergeSeurat.filter,label = T, group.by = "GeneralCell_Mouse", repel = T, split.by = "orig.ident")

DimPlot(mergeSeurat.filter,label = T, group.by = "PreciseCell_Mouse", repel = T)
DimPlot(mergeSeurat.filter,label = T, group.by = "PreciseCell_Mouse", repel = T, split.by = "orig.ident")


### ensure ident to plot is active
Idents(seurat_obj) <- "sample"

### captures count of cells for the ident with the fewest
maxcells  <- min(table(Idents(seurat_obj)))
                 
### nested object subsetting with downsampling
DoHeatmap(subset(seurat_obj, downsample = maxcells), features = features2plot, etc.)
                 
                 


library("atlasapprox")
mouse_organs <- GetOrgans(organism = 'm_musculus')
print(mouse_organs)

library(ggplot2)

# Convertir tabla a data.frame
tabla_df <- as.data.frame(table(mergeSeurat.filter$GeneralCell_HPCA))

# Crear gráfico
ggplot(tabla_df, aes(x = mergeSeurat.filter$GeneralCell_HPCA, y = Freq)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Frecuencia por Grupo", x = "Grupo", y = "Frecuencia") +
  theme_minimal()


##2 Quality####


#Filter
C7_sc.filter <- subset(C7_sc, subset = nFeature_RNA > 200 & nFeature_RNA < 6500 & percent.mt < 15)

VlnPlot(C7_sc.filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001)

plot1 <- FeatureScatter(C7_sc.filter, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position="none")
plot2 <- FeatureScatter(C7_sc.filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
plot1 + plot2


#High expressed genes
par(mar = c(4, 8, 2, 1))
C <- C7_sc.filter@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
#highly expressed "normal" genes
boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main="High expressed genes after filtering")



#3 Reunion 020425####
a<-DimPlot(mergeSeurat.filter,group.by = "seurat_clusters", label = T) 
manual_anot<-read.table("../manual_anottation.txt", header = T, sep = "\t")
manual_anot_unique <- manual_anot[!duplicated(manual_anot$gen) & !duplicated(manual_anot$gen, fromLast = TRUE), ]

for (type in unique(manual_anot$etapa)){
  pdf(file=paste("Cell_type_", type, "_manual_etapa.pdf", sep = ""), width = 12, height = 8)
  genes = manual_anot[manual_anot$etapa == type, 1]
  b = FeaturePlot(mergeSeurat.filter, features = genes)
  plot(a  + ggtitle(paste("Cell type",type)) + b)
  dev.off()
}

for (type in unique(manual_anot$big_etapa)){
  pdf(file=paste("Cell_type_", type, "_manual_Big_etapa.pdf", sep = ""), width = 24, height = 8)
  cat(paste("Big etapa:", type, "\n"))
  genes = manual_anot[manual_anot$big_etapa == type, 1]
  b = FeaturePlot(mergeSeurat.filter, features = genes)
  plot(a  + ggtitle(paste("Cell type",type)) + b)
  dev.off()
}

a<-DimPlot(mergeSeurat.filter,group.by = "seurat_clusters", label = T)

b<-FeaturePlot(mergeSeurat.filter, features = c("H2-Eb1"))
#  ,"Klrb1c","Klrk1"))
b<-FeaturePlot(mergeSeurat.filter, features = c("Klrc2","Klrb1c","Klrk1"))
c<-FeaturePlot(mergeSeurat.filter, features = c("Zfp990"))
a + b
a + c

#4. Anotacion Manual####

# Cargar librerías
library(Seurat)
library(dplyr)

# 1. Leer el archivo CSV

anotaciones <- read.csv("anotaciones.tsv", stringsAsFactors = FALSE, sep = "\t")

# 2. Convertir cluster a carácter (para evitar problemas de indexación)
anotaciones$cluster <- as.character(anotaciones$cluster)

# 3. Crear un data frame con las anotaciones para cada célula
clusters_obj <- Idents(mergeSeurat.filter)  # seurat_obj es tu objeto Seurat
clusters_obj <- as.character(clusters_obj)

# 4. Mapear anotaciones según el cluster
etapa <- anotaciones$etapa[match(clusters_obj, anotaciones$cluster)]
big_etapa <- anotaciones$big_etapa[match(clusters_obj, anotaciones$cluster)]

# 5. Añadir al objeto Seurat en meta.data
mergeSeurat.filter$manual_stage<-NA
mergeSeurat.filter$manual_stage <- etapa
mergeSeurat.filter$manual_big_stage<-NA
mergeSeurat.filter$manual_big_stage <- big_etapa

# Comprobar
head(mergeSeurat.filter@meta.data)

a <- DimPlot(mergeSeurat.filter, group.by = "manual_stage", label = T)
x <- DimPlot(mergeSeurat.filter, group.by = "manual_big_stage", label = T)
a + x
a <- DimPlot(mergeSeurat.filter, group.by = "manual_stage", label = T, split.by = "orig.ident")
b <- FeaturePlot(mergeSeurat.filter, features = "Zfp990",split.by = "orig.ident")
a+b

DotPlot(mergeSeurat.filter, features = c("Zfp990", "Kif1b","Rab6b"), group.by = "manual_stage", split.by = "orig.ident")
