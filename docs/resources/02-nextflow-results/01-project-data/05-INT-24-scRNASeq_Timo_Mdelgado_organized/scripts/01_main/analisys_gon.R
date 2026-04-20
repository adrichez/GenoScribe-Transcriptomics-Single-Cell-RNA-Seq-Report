setwd("/mnt/beegfs/gariza/Projects/scRNASeq/")


#Si hay problemas con matrixStats
rm(list = ls())
.libPaths("/mnt/beegfs/gariza/Rlibs")
options(repos = c(CRAN="https://cloud.r-project.org"))
"matrixStats" %in% loadedNamespaces()
packageVersion("matrixStats")
remotes::install_version("matrixStats", version = "1.1.0", lib="~/Rlibs", force=TRUE)
library(matrixStats, lib.loc="~/Rlibs")
remotes::install_github("mojaveazure/seurat-disk", dependencies = FALSE)

library(Seurat)
library(Matrix)
library(dplyr)
library(BiocParallel)
library(scater)
library(ggplot2)
library(celldex)
library(SeuratData)
library(patchwork)


my.date<-function(){
  return(format(Sys.time(), "%d-%m-%Y %H:%M:%S"))
}
DEG_by_cluster_WT_vs_KO <- function(Obj, cluster_column="seurat_clusters", 
                                    condition_column="condition", label="Cluster") {
  
  # Asegurarse de usar RNA
  DefaultAssay(Obj) <- "RNA"
  
  for (cl in unique(Obj@meta.data[[cluster_column]])) {
    cat(paste0("[", Sys.time(), "] Analizando cluster ", cl, "\n"))
    
    # Hacer subset del cluster
    sub_obj <- tryCatch({
      cells_in_cluster <- rownames(Obj@meta.data)[Obj@meta.data[[cluster_column]] == cl]
      if (length(cells_in_cluster) == 0) stop("No se encontraron células en este cluster")
      subset(Obj, cells = cells_in_cluster)
    }, error = function(e) {
      message("❌ Error al hacer subset del cluster ", cl, ": ", e$message)
      return(NULL)
    })
    if (is.null(sub_obj)) next
    if (ncol(sub_obj) <= 10) {
      cat("Cluster ", cl, " tiene menos de 10 células. Saltando...\n")
      next
    }
    
    # Limpiar niveles duplicados en condition
    sub_obj@meta.data[[condition_column]] <- droplevels(as.factor(sub_obj@meta.data[[condition_column]]))
    conds <- unique(sub_obj@meta.data[[condition_column]])
    if (length(conds) < 2) {
      cat("Cluster ", cl, " no tiene ambas condiciones. Saltando...\n")
      next
    }
    
    # Definir identidades
    Idents(sub_obj) <- condition_column
    
    # Calcular DEGs
    markers <- tryCatch({
      FindMarkers(sub_obj, ident.1="WT", ident.2="KO", min.pct = 0, logfc.threshold = 0)
    }, error = function(e) {
      message("❌ Error en FindMarkers cluster ", cl, ": ", e$message)
      return(NULL)
    })
    if (is.null(markers)) next
    
    # Guardar tabla completa
    write.table(markers[order(-abs(markers$avg_log2FC)), ],
                file = paste0("DEG_cluster_", cl, "_WT_vs_KO_", label, ".tsv"),
                sep="\t", quote=F, col.names=NA)
    
    # Filtrar DEGs solo por log2FC distinto de 0
    markers_sig <- markers[abs(markers$avg_log2FC) > 0, ]
    if (nrow(markers_sig) == 0) {
      message("Cluster ", cl, " no tiene DEGs (log2FC=0)")
    }
  }
  
  cat("¡DEGs guardados con éxito!\n")
}

DEG_by_cluster_WT_vs_KO <- function(Obj, cluster_column="seurat_clusters", 
                                    condition_column="condition", label="Cluster") {
  
  DefaultAssay(Obj) <- "RNA"
  
  for (cl in unique(Obj@meta.data[[cluster_column]])) {
    cat(paste0("[", Sys.time(), "] Analizando cluster ", cl, "\n"))
    
    # Subset del cluster
    sub_obj <- subset(Obj, subset = !!sym(cluster_column) == cl)
    
    if (ncol(sub_obj) <= 10) next
    
    conds <- unique(sub_obj@meta.data[[condition_column]])
    if (length(conds) < 2) {
      cat("Cluster ", cl, " no tiene ambas condiciones. Saltando...\n")
      next
    }
    
    Idents(sub_obj) <- condition_column
    
    # DEGs WT vs KO
    markers <- FindMarkers(sub_obj, ident.1="WT", ident.2="KO")
    
    # Guardar tabla
    write.table(markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)),],
                file=paste0("DEG_cluster_", cl, "_WT_vs_KO_", label, ".tsv"),
                sep="\t", quote=F, col.names=NA)
    
    # Filtrar genes significativos
    markers_sig <- markers[markers$p_val_adj <= 0.05 & abs(markers$avg_log2FC) > 0,]
    if (nrow(markers_sig) > 0) {
      top_up <- rownames(head(markers_sig[order(markers_sig$avg_log2FC, decreasing = TRUE),], 5))
      top_down <- rownames(head(markers_sig[order(markers_sig$avg_log2FC, decreasing = FALSE),], 5))
      
      pdf(file = paste0("DEG_cluster_", cl, "_WT_vs_KO_", label, "_plots.pdf"), height=11, width=8.5)
      
      # DimPlot
      print(DimPlot(sub_obj, group.by=condition_column, label=TRUE, repel=TRUE) + ggtitle(paste("Cluster", cl)))
      
      # DotPlot
      print(DotPlot(sub_obj, features = c(top_up, top_down), group.by=condition_column) +
              coord_flip() + ggtitle(paste("Top genes Cluster", cl)))
      
      # FeaturePlot: 2 genes de cada lado
      genes_to_plot <- c(head(top_up,2), head(top_down,2))
      print(FeaturePlot(sub_obj, features=genes_to_plot, split.by=condition_column))
      
      dev.off()
    }
  }
  cat("¡DEGs y plots guardados con éxito!\n")
}

Markers_specific <- function(Obj, group="seurat_clusters") {
  
  Idents(Obj) <- group
  
  # 1. Encontrar DEGs para todos los clusters vs todos
  all <- FindAllMarkers(
    Obj,
    only.pos = TRUE,
    logfc.threshold = 0,        # No bloquees logFC pequeños
    min.pct = 0,
    test.use = "wilcox"
  )
  
  # 2. Aplicar tus filtros
  all <- all[all$p_val_adj < 0.05 & all$avg_log2FC > 0, ]
  
  # 3. QUEDARSE SOLO LOS GENES QUE SALEN EN UN ÚNICO CLUSTER
  unique_markers <- all[!duplicated(all$gene) & !duplicated(all$gene, fromLast=TRUE), ]
  
  return(unique_markers)
}


# ruta al directorio con matrix.mtx, barcodes.tsv, features.tsv
data_dir_wt <- "/mnt/beegfs/gariza/Projects/scRNASeq/WT/"
data_dir_ko <- "/mnt/beegfs/gariza/Projects/scRNASeq/KO/"

# Leer datos
ko_counts <- Read10X(data.dir = data_dir_ko)
wt_counts <- Read10X(data.dir = data_dir_wt)

# Crear objetos
ko <- CreateSeuratObject(counts = ko_counts, project = "KO")
wt <- CreateSeuratObject(counts = wt_counts, project = "WT")

# Etiqueta de condición
ko$condition <- "KO"
wt$condition <- "WT"

# Merge / combinar
combined <- merge(ko, y = wt, add.cell.ids = c("KO","WT"), project = "MouseSc")

# QC: porcentaje mitocondrial (para ratón, patrón "^mt-", para humano "^MT-")
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

# Visualizar QC
VlnPlot(combined, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
plot1 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filtrar células según el QC (umbrales sugeridos; ajusta si ves que pierdes demasiadas células)
combined <- subset(combined, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 20)

# Normalización + corrección variable técnica (SCTransform)
combined <- SCTransform(combined, verbose = FALSE)

# Muestra los 10 genes más variables calculados por SCTransform
top10 <- head(VariableFeatures(combined), 10)
plotVariable<-VariableFeaturePlot(combined, assay = "SCT")
top10
# FeaturePlot para esos genes
FeaturePlot(combined, features = top10)


# PCA
combined <- RunPCA(combined, verbose = FALSE)
ElbowPlot(combined)

# Vecindario, clustering
combined <- FindNeighbors(combined, dims = 1:18) #La dimension se determina segun el aplanamiento del ElboePlot
combined <- FindClusters(combined, resolution = c(0.3, 0.5, 0.7, 1))
DimPlot(combined, group.by = "SCT_snn_res.0.5", label = TRUE)
combined <- FindClusters(combined, resolution = 0.5)

DimHeatmap(combined, dims = 1:10, nfeatures = 15, balanced = TRUE)


# Visualizar UMAP

combined <- RunUMAP(combined, dims = 1:18)
DimPlot(combined, reduction = "umap", group.by = "condition")
DimPlot(combined, reduction = "umap", label = TRUE)
DimPlot(combined, reduction = "umap", split.by = "condition")


# Encontrar marcadores por cluster
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
write.csv(markers, "markers_by_cluster.csv")

# Top genes por cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Heatmap de genes marcador
DoHeatmap(combined, features = top_markers$gene) + NoLegend()

# Comparación KO vs WT globalmente o por cluster
Idents(combined) <- "condition"
ko_vs_wt <- FindMarkers(combined, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0.25)
head(ko_vs_wt)
write.csv(ko_vs_wt, "markers_KO_vs_WT.csv")

# Visualización de genes diferenciales
top_ko_vs_wt <- head(ko_vs_wt[order(-ko_vs_wt$avg_log2FC), ], 10)
FeaturePlot(combined, features = rownames(top_ko_vs_wt))
DotPlot(combined, features = unique(top_markers$gene)) + RotatedAxis()

#Reference data
ref <- celldex::MouseRNAseqData()

#SingleR & celldex libraries
suppressPackageStartupMessages(library("SingleR", lib.loc = "/mnt/beegfs/Bioinformatica/bin/R/4.1.2/lib/R/library/"))
suppressPackageStartupMessages(library("celldex"))
ncores = 4

##4 cell anotattion
sobj_counts <- suppressWarnings(as.SingleCellExperiment(combined))
sobj_counts <- suppressWarnings(logNormCounts(sobj_counts))

suppressMessages(ref <- celldex::MouseRNAseqData())
cat(paste("[ " , my.date() , " ]"," Anotatting using MouseRNASeq\n",sep=""))
pred <- SingleR(test=sobj_counts, ref=ref, labels=ref$label.main, BPPARAM=MulticoreParam(ncores))
combined[["GeneralCell_Mouse"]] <- pred$labels

DimPlot(combined, group.by = "GeneralCell_Mouse", label = TRUE, repel = TRUE, split.by = "condition") +
  ggtitle("Cell Type Annotation (SingleR + MouseRNAseqData)")



#Otra normalizacion
combined <- NormalizeData(combined)
str(combined)

combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined)

plotV <- VariableFeaturePlot(combined)
LabelPoints(plot = plotV, points = top10, repel = TRUE)


#Downsampling de KO
ncol(wt)

n_cells <- ncol(wt)

set.seed(123)  # para reproducibilidad
ko_ds <- subset(ko, cells = sample(colnames(ko), n_cells))

# Luego combinas normalmente
combined_ds <- merge(ko_ds, y = wt, add.cell.ids = c("KO", "WT"), project = "MouseSc")

setwd("/mnt/beegfs/gariza/Projects/scRNASeq/Analisis_downsampling/Sin_filtrar")

# QC: porcentaje mitocondrial (para ratón, patrón "^mt-", para humano "^MT-")
combined_ds[["percent.mt"]] <- PercentageFeatureSet(combined_ds, pattern = "^mt-")

# Visualizar QC
VlnPlot(combined_ds, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)
plot1 <- FeatureScatter(combined_ds, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined_ds, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filtrar células según el QC (umbrales sugeridos; ajusta si ves que pierdes demasiadas células)
combined_ds <- subset(combined_ds, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 & percent.mt < 20)

# Normalización + corrección variable técnica (SCTransform)
combined_ds <- SCTransform(combined_ds, verbose = FALSE)

# Muestra los 10 genes más variables, en este caso quitar la normalizacion de SCTransform y usar NormalizeData en su lugar para evitar incompatibilidades
combined_ds <- NormalizeData(combined_ds)
combined_ds <- FindVariableFeatures(combined_ds, selection.method = "vst", nfeatures = 2000)
combined_ds <- ScaleData(combined_ds)

plotVariable <- VariableFeaturePlot(combined_ds)
top10 <- head(VariableFeatures(combined_ds), 10)
LabelPoints(plot = plotVariable, points = top10, repel = TRUE)

# PCA
combined_ds <- RunPCA(combined_ds, verbose = FALSE)
ElbowPlot(combined_ds)

# Vecindario, clustering
combined_ds <- FindNeighbors(combined_ds, dims = 1:18) #La dimensión se determina según el aplanamiento del ElbowPlot
combined_ds <- FindClusters(combined_ds, resolution = c(0.3, 0.5, 0.7, 1))
DimPlot(combined_ds, group.by = "SCT_snn_res.0.5", label = TRUE)
combined_ds <- FindClusters(combined_ds, resolution = 0.5)

DimHeatmap(combined_ds, dims = 1:10, nfeatures = 15, balanced = TRUE)

# Visualizar UMAP
combined_ds <- RunUMAP(combined_ds, dims = 1:18)
DimPlot(combined_ds, reduction = "umap", group.by = "condition")
DimPlot(combined_ds, reduction = "umap", label = TRUE)
DimPlot(combined_ds, reduction = "umap", split.by = "condition")

# Encontrar marcadores por cluster
markers <- FindAllMarkers(combined_ds, min.pct = 0.1, logfc.threshold = 0.2)
head(markers)
write.csv(markers, "markers_by_cluster.csv")

# Top genes por cluster
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Heatmap de genes marcador
DoHeatmap(combined_ds, features = top_markers$gene) + NoLegend()

# Comparación KO vs WT globalmente o por cluster
Idents(combined_ds) <- "condition"
ko_vs_wt_ds <- FindMarkers(combined_ds, ident.1 = "KO", ident.2 = "WT", min.pct = 0.05, logfc.threshold = 0)
ko_vs_wt_ds
write.csv(ko_vs_wt_ds, "markers_KO_vs_WT_ds_all.csv")
FeaturePlot(combined_ds, features = "Zfp990")

# Visualización de genes diferenciales
top_ko_vs_wt_ds <- head(ko_vs_wt_ds[order(-ko_vs_wt_ds$avg_log2FC), ], 10)
FeaturePlot(combined_ds, features = rownames(top_ko_vs_wt_ds))
DotPlot(combined_ds, features = unique(top_markers$gene)) + RotatedAxis()


#Integracion de datos
combined.list <- SplitObject(combined_ds, split.by = "condition")
combined.list <- lapply(combined.list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = combined.list)
combined.anchors <- FindIntegrationAnchors(object.list = combined.list, 
                                           anchor.features = features)

combined_integrated <- IntegrateData(anchorset = combined.anchors)
DefaultAssay(combined_integrated) <- "integrated"

# Escalado de los datos integrados
combined_integrated <- ScaleData(combined_integrated)

# PCA
combined_integrated <- RunPCA(combined_integrated, verbose = FALSE)

# UMAP
combined_integrated <- RunUMAP(combined_integrated, dims = 1:18)

# Vecinos y clustering
combined_integrated <- FindNeighbors(combined_integrated, dims = 1:18)
combined_integrated <- FindClusters(combined_integrated, resolution = 0.5)

DimPlot(combined_integrated, group.by = c("condition","seurat_clusters"), label = T)
DimPlot(combined_integrated, split.by = "condition", label = T)

#dittoSeq
combined_integrated$GeneralCell_Mouse <- combined_ds$GeneralCell_Mouse
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dittoSeq")
library(dittoSeq)

combined_ds$cell_anot <- timo_tabula_seurat$cell_anot[match(colnames(combined_ds),
                                                            colnames(timo_tabula_seurat))]
combined_ds$GeneralCell_Mouse <- combined$GeneralCell_Mouse[match(colnames(combined_ds),
                                                                  colnames(combined))]


vars <- c("seurat_clusters")

# Barplots con escala 'count'
for(v in vars){
  print(dittoBarPlot(combined_integrated, var = v, group.by = "condition", scale = "count"))
}

# Barplots sin escala
for(v in vars){
  print(dittoBarPlot(combined_integrated, var = v, group.by = "condition"))
}

# Frequency plots
for(v in vars){
  print(dittoFreqPlot(combined_integrated, var = v, sample.by = "condition", group.by = "condition"))
}

dev.off()

meta <- combined_integrated@meta.data
tab_counts <- table(meta$seurat_clusters, meta$condition)

results <- data.frame(
  cluster   = rownames(tab_counts),
  WT_count  = NA,
  KO_count  = NA,
  log2FC    = NA,
  p_value   = NA,
  padj      = NA
)

for (cl in rownames(tab_counts)) {
  
  mat <- matrix(c(
    tab_counts[cl, "WT"],
    tab_counts[cl, "KO"],
    sum(tab_counts[-which(rownames(tab_counts)==cl), "WT"]),
    sum(tab_counts[-which(rownames(tab_counts)==cl), "KO"])
  ), nrow = 2, byrow = TRUE)
  
  p <- fisher.test(mat)$p.value
  
  WT <- tab_counts[cl, "WT"]
  KO <- tab_counts[cl, "KO"]
  
  results[results$cluster == cl, "WT_count"] <- WT
  results[results$cluster == cl, "KO_count"] <- KO
  results[results$cluster == cl, "p_value"]  <- p
}

# Ajuste por FDR
results$padj <- p.adjust(results$p_value, method = "BH")

# Cálculo del logFC basado en los conteos
results$log2FC <- log2((results$KO_count) / (results$WT_count))
results

write.csv(results, file = "cluster_differences.csv", row.names = FALSE)

#Markers por cluster
specific <- Markers_specific(combined_integrated)
head(specific)
specific
write.csv(specific, "markers_by_cluster.csv")

#DEGs por cluster
DEG_by_cluster_WT_vs_KO(
  Obj = combined_integrated,
  cluster_column = "seurat_clusters",
  condition_column = "condition",
  label = "Integrated"
)



# Reference data para anotacion
ref <- celldex::MouseRNAseqData()

# SingleR & celldex libraries
suppressPackageStartupMessages(library("SingleR", lib.loc = "/mnt/beegfs/Bioinformatica/bin/R/4.1.2/lib/R/library/"))
suppressPackageStartupMessages(library("celldex"))
ncores = 4

## 4 cell annotation
sobj_counts <- suppressWarnings(as.SingleCellExperiment(combined_ds))
sobj_counts <- suppressWarnings(logNormCounts(sobj_counts))

suppressMessages(ref <- celldex::MouseRNAseqData())
cat(paste("[ " , my.date() , " ]"," Annotating using MouseRNASeq\n",sep=""))
options(matrixStats.useNames.NA = "deprecated")
pred <- SingleR(test = sobj_counts, ref = ref, labels = ref$label.main, BPPARAM = MulticoreParam(ncores))
combined_ds[["GeneralCell_Mouse"]] <- pred$labels

DimPlot(combined_ds, group.by = "GeneralCell_Mouse", label = TRUE, repel = TRUE, split.by = "condition") +
  ggtitle("Cell Type Annotation (SingleR + MouseRNAseqData)")

#Comparar T-cells
# Define las identidades como tipo celular
Idents(combined_ds) <- "GeneralCell_Mouse"

# Crea un subobjeto con solo las T cells
t_cells <- subset(combined_ds, idents = "T cells")

# Dentro de ese subobjeto, define la condición (KO/WT) como identidad activa
Idents(t_cells) <- "condition"

# Comprueba cuántas T cells hay por condición
table(t_cells$condition)

tcell_DEGs <- FindMarkers(
  t_cells,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.1
)
write.csv(tcell_DEGs, "Tcells_KO_vs_WT_all.csv")

#Tabula muris
tabula_muris<-function(Obj){
  suppressWarnings({
    
    # Cargar las librerías necesarias
    #library(Seurat)       # Para trabajar con datos single-cell
    #library(SeuratDisk)   # Para conversión si se necesita en un futuro
    library(dplyr)        # Para manipulación de datos
    
    # 1. Cargar los datos de Tabula Muris en formato .mtx, barcodes.tsv y features.tsv
    # Define la ruta donde se encuentran los archivos
    ruta_tabula <- "/mnt/beegfs/Bioinformatica/INT/05-INT-24-scRNASeq_Timo_Mdelgado/Analisis/Anotacion/Thymus-10X_P7_11"
    
    # Leer los datos en formato estándar de 10x (matriz de expresión)
    timo_tabula <- Read10X(data.dir = ruta_tabula)
    
    # Crear un objeto Seurat a partir de los datos cargados
    timo_tabula_seurat <- CreateSeuratObject(counts = timo_tabula, 
                                             project = "TabulaMuris_Timo", 
                                             min.cells = 3, 
                                             min.features = 200)
    
    # 2. Preprocesar los datos de Tabula Muris
    # Filtrar células y genes con pocas cuentas (ruido)
    timo_tabula_seurat <- subset(timo_tabula_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
    
    # Normalizar los datos (log-normalización para homogeneizar las cuentas)
    timo_tabula_seurat <- NormalizeData(timo_tabula_seurat, verbose = FALSE)
    
    # Identificar genes con alta variabilidad
    timo_tabula_seurat <- FindVariableFeatures(timo_tabula_seurat, selection.method = "vst", nfeatures = 2000)
    
    # Escalar los datos para eliminar el efecto de la profundidad de secuenciación
    timo_tabula_seurat <- ScaleData(timo_tabula_seurat, verbose = FALSE)
    
    # Reducción de dimensionalidad con PCA (Análisis de Componentes Principales)
    timo_tabula_seurat <- RunPCA(timo_tabula_seurat, npcs = 30, verbose = FALSE)
    
    # Generar un mapa UMAP para visualización
    timo_tabula_seurat <- RunUMAP(timo_tabula_seurat, reduction = "pca", dims = 1:18)
    
    # Realizar clustering para identificar grupos de células
    timo_tabula_seurat <- FindNeighbors(timo_tabula_seurat, dims = 1:18)
    timo_tabula_seurat <- FindClusters(timo_tabula_seurat, resolution = 0.5)
    
    # Visualizar los clusters en UMAP
    DimPlot(timo_tabula_seurat, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5) + 
      ggtitle("Tabula Muris - Thymus Data")
    
    #3 Anotar Tabula segun la web
    df <- read.csv("/mnt/beegfs/Bioinformatica/INT/05-INT-24-scRNASeq_Timo_Mdelgado/Analisis/Anotacion/annotations_droplet_Thymus.csv", header = F)
    df$V1<-paste0(df$V1, "-1")
    timo_tabula_seurat[["cell_anot"]] <- NA
    timo_tabula_seurat$cell_anot <- df$V3[match(colnames(timo_tabula_seurat), df$V1)]
    
    # 4. Transferir anotaciones de Tabula Muris a tus datos
    # Encontrar "anclajes" entre los datos de referencia y tus datos
    anchors <- FindTransferAnchors(reference = timo_tabula_seurat, 
                                   query = Obj, 
                                   dims = 1:18)
    
    # Transferir las etiquetas de anotación desde Tabula Muris a tus datos
    predictions <- TransferData(anchorset = anchors, 
                                refdata = timo_tabula_seurat$cell_anot, 
                                dims = 1:18)
    
    # Añadir las anotaciones transferidas como metadatos a tu objeto Seurat
    Obj <- AddMetaData(Obj, metadata = predictions)
    
    # 5. Visualizar las anotaciones en tus datos
    # Generar un UMAP para ver las células anotadas
    Obj <- RunUMAP(Obj, reduction = "pca", dims = 1:18)
    Obj$cell_anot<-Obj$predicted.id
    Obj[["predicted.id"]] <- NULL
    
    return(Obj)
  })
}

combined_ds_tabula <- tabula_muris(combined_ds)
DimPlot(combined_ds_tabula, group.by = "cell_anot", label = TRUE, repel = TRUE)

#Expresion diferencial
t_cells <- subset(combined_ds_tabula, subset = cell_anot == "T cell")

markers_WT_vs_KO <- FindMarkers(
  t_cells,
  ident.1 = "KO",
  ident.2 = "WT",
  group.by = "condition",
  logfc.threshold = 0.25)  

write.csv(markers_WT_vs_KO, "Tcells_KO_vs_WT_tabula_muris.csv")

FeaturePlot(combined_ds, features = "Cenps")

WT <- subset(combined_integrated, subset = condition == "WT")
KO <- subset(combined_integrated, subset = condition == "KO")

markers_WT <- FindAllMarkers(WT, only.pos=TRUE)
markers_KO <- FindAllMarkers(KO, only.pos=TRUE)

####Reunión marina_mario 13/01/25 -> agrupar clusters 0,1 y 2,6, ademas de renombrarlos####
table(combined_ds$seurat_clusters)

# Crear copia de la columna de clusters
combined_ds$cell_annot <- as.character(combined_ds$seurat_clusters)

# Unir clusters 0 y 1 como "Cluster_0_1"
combined_ds$cell_annot[combined_ds$seurat_clusters %in% c(0, 1)] <- "DPre"

# Unir clusters 2 y 6 como "Cluster_2_6"
combined_ds$cell_annot[combined_ds$seurat_clusters %in% c(2, 6)] <- "Cluster 2-6"
combined_ds$cell_annot[combined_ds$seurat_clusters == 3] <- "DPsel"
combined_ds$cell_annot[combined_ds$seurat_clusters == 4] <- "SP"
combined_ds$cell_annot[combined_ds$seurat_clusters == 5] <- "DPblast"
combined_ds$cell_annot[combined_ds$seurat_clusters == 7] <- "preDP (ISP)"
combined_ds$cell_annot[combined_ds$seurat_clusters == 8] <- "CD8 exhausted, naïve"
combined_ds$cell_annot[combined_ds$seurat_clusters == 9] <- "Cluster 9"
combined_ds$cell_annot[combined_ds$seurat_clusters == 10] <- "Monocyte, Macrophages, Dendritic"
combined_ds$cell_annot[combined_ds$seurat_clusters == 11] <- "DN, Precursors"
combined_ds$cell_annot[combined_ds$seurat_clusters == 12] <- "Stroma"
combined_ds$cell_annot[combined_ds$seurat_clusters == 13] <- "Cluster 13"
combined_ds$cell_annot[combined_ds$seurat_clusters == 14] <- "Cluster 14"

DimPlot(combined_ds, group.by = "cell_annot", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP con clusters combinados")

# Comparación KO vs WT globalmente o por cluster
combined_ds$cell_annot <- factor(combined_ds$cell_annot)
Idents(combined_ds) <- combined_ds$cell_annot
cluster_ids <- levels(combined_ds$cell_annot)
all_markers <- list()

for (clus in cluster_ids) {
  # Subset del cluster
  cluster_cells <- subset(combined_ds, idents = clus)
  
  # Comparación KO vs WT dentro del cluster
  markers <- FindMarkers(cluster_cells,
                         ident.1 = "KO",
                         ident.2 = "WT",
                         group.by = "condition")
  
  all_markers[[clus]] <- markers
}

#Otra forma
DEG_by_cluster_KO_vs_WT <- function(
    Obj,
    cluster_column = "cell_annot",
    condition_column = "condition",
    label = "cell_annot"
) {
  
  DefaultAssay(Obj) <- "RNA"
  
  for (cl in unique(Obj@meta.data[[cluster_column]])) {
    cat(paste0("[", Sys.time(), "] Analizando cluster ", cl, "\n"))
    
    # Subset del cluster (robusto)
    cells_use <- rownames(Obj@meta.data[Obj@meta.data[[cluster_column]] == cl, ])
    sub_obj <- subset(Obj, cells = cells_use)
    
    if (ncol(sub_obj) <= 10) next
    
    conds <- unique(sub_obj@meta.data[[condition_column]])
    if (length(conds) < 2) {
      cat("Cluster ", cl, " no tiene ambas condiciones. Saltando...\n")
      next
    }
    
    # Condición como identidad
    Idents(sub_obj) <- sub_obj@meta.data[[condition_column]]
    
    # DEGs KO vs WT
    markers <- FindMarkers(sub_obj, ident.1 = "KO", ident.2 = "WT")
    
    # Guardar tabla completa
    write.table(
      markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ],
      file = paste0("DEG_", cl, "_KO_vs_WT_", label, ".tsv"),
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
    
    # Filtrar significativos
    markers_sig <- markers[
      markers$p_val_adj <= 0.05 & abs(markers$avg_log2FC) > 0,
    ]
    
    if (nrow(markers_sig) > 0) {
      
      top_up   <- rownames(head(markers_sig[order(markers_sig$avg_log2FC, decreasing = TRUE), ], 5))
      top_down <- rownames(head(markers_sig[order(markers_sig$avg_log2FC, decreasing = FALSE), ], 5))
      
      pdf(file = paste0("DEG_", cl, "_KO_vs_WT_", label, "_plots.pdf"),
          height = 11, width = 8.5)
      
      print(
        DimPlot(sub_obj, group.by = condition_column, label = TRUE, repel = TRUE) +
          ggtitle(paste("Cluster:", cl))
      )
      
      print(
        DotPlot(sub_obj, features = c(top_up, top_down), group.by = condition_column) +
          coord_flip() +
          ggtitle(paste("Top genes –", cl))
      )
      
      genes_to_plot <- c(head(top_up, 2), head(top_down, 2))
      print(
        FeaturePlot(sub_obj, features = genes_to_plot, split.by = condition_column)
      )
      
      dev.off()
    }
  }
  
  cat("✅ DEGs y plots guardados con éxito\n")
}

