#Loading usual packages

suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("celldex"))
suppressPackageStartupMessages(library("SingleR"))
suppressPackageStartupMessages(library("scuttle"))
suppressPackageStartupMessages(library("clusterProfiler"))
suppressMessages(require(BiocParallel))
suppressPackageStartupMessages(library(future))
#suppressPackageStartupMessages(library(cli))
#suppressPackageStartupMessages(library(purrr))
#suppressPackageStartupMessages(library(magrittr))
#modifying parameters for scRNASeq

options(ggrepel.max.overlaps = Inf)
options(future.globals.maxSize = 8000 * 1024^2)
my.date<-function(){
  return(format(Sys.time(), "%d-%m-%Y %H:%M:%S"))
}
readscRDS<-function(FullPath, pattern=".rds", filter=".rds", csv=F){
  cat(paste("[ " , my.date() , " ]"," ##########READING INPUT FILES###########\n",sep=""))
  
  for (file in list.files(path=FullPath, pattern = pattern, full.names = T)){
    cat(paste("[ " , my.date() , " ]"," Reading ....", file, "\n", sep=""))
    var=gsub(filter,"",basename(file))
    obj = readRDS(file)
    if (csv) {
      myfile=gsub(filter,"",file)
      Ids = read.csv(paste(myfile, ".csv", sep = ""), header = T, stringsAsFactors = FALSE)
      #print(head(Ids))
      Ids <- Ids[!grepl("doublet|unassigned", Ids$code, ignore.case = TRUE), ]
      #print(head(Ids))
      obj@meta.data$Paciente <- Ids$code[match(rownames(obj@meta.data), Ids$barcode)]
      obj <- obj[, complete.cases(obj@meta.data)]
      
    }
    assign(var,obj, envir = .GlobalEnv)
  }
}
mergescRDS<-function(pattern, project="MyscProject"){
  # Obtener todos los objetos en el entorno
  all_objects <- ls(envir = .GlobalEnv, pattern = pattern)
  cat(paste("[ " , my.date() , " ]"," ###########MERGING DATASETS############\n",sep=""))
  
  # Filtrar los objetos cuyo nombre contenga la palabra "pool" y que sean de clase Seurat
  int_objects <- lapply(all_objects, function(x) {
    if (inherits(get(x), "Seurat") && grepl(pattern, x, ignore.case = TRUE)) {
      return(get(x))  # Retornar el objeto Seurat completo
    }
  })
  
  # Eliminar elementos NULL (si hay objetos que no son Seurat o no cumplen las condiciones)
  int_objects <- Filter(Negate(is.null), int_objects)
  
  if(length(int_objects)>1){
    Obj_merge<-merge(x = int_objects[[1]], y = int_objects[-1], add.cell.ids = ls(pattern = pattern, envir = .GlobalEnv), project=project)
    return(Obj_merge)
  }else{
    cat("You need at lest to objetcs to merge them")
    exit(1)
  }
  
}
quality<-function(Obj, type="Pre", ndims=0){
  outfile = paste(type,"_Quality_", deparse(substitute(Obj)), ".pdf",sep="")
  pdf(file = outfile, width = 11.69, height = 8.27)
  
  suppressWarnings({
    cat(paste("[ " , my.date() , " ]"," ################QUALITY#################\n",sep=""))
    
    cat(paste("[ " , my.date() , " ]"," Calculating percent of mt genes\n",sep=""))
    
    Obj_tmp<-RenameGenesSeurat(obj = Obj, newnames = toupper(rownames(Obj)))
    Obj[["percent.mt"]] <- PercentageFeatureSet(Obj_tmp, pattern = "^MT-")
    
    z<-VlnPlot(Obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.0001, group.by = "orig.ident")
    plot(z)
    
    # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
    # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
    cat(paste("[ " , my.date() , " ]"," Calculating relevant features to check\n",sep=""))
    plot1 <- FeatureScatter(Obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
      theme(legend.position="none")
    plot2 <- FeatureScatter(Obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
      theme(legend.position="none")
    plot(plot1 + plot2)
    
    #High expressed genes
    cat(paste("[ " , my.date() , " ]"," Calculating high expressed genes\n",sep=""))
    C <- Obj@assays$RNA@counts
    C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
    most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
    #highly expressed "normal" genes
    cat(paste("[ " , my.date() , " ]"," plot\n",sep=""))
    par(mar = c(4, 8, 2, 1))
    suppressWarnings(boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
                             col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main="High expressed genes"))
    
    if (type == "Pre"){
      cat(paste("[ " , my.date() , " ]"," Normalizing and Scaling data\n",sep=""))
      Obj <- NormalizeData(Obj, verbose=F)
      Obj <- FindVariableFeatures(Obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
      
      # Identify the 10 most highly variable genes
      top10 <- head(VariableFeatures(Obj), 10)
      
      # plot variable features with and without labels
      plot1 <- VariableFeaturePlot(Obj) + 
        theme(legend.position="top")
      plot2 <- suppressMessages(LabelPoints(plot = plot1, points = top10, repel = TRUE)) + 
        theme(legend.position="none")
      suppressWarnings(plot(plot1 + plot2))
      
      #Scaling
      all.genes <- rownames(Obj)
      Obj <- ScaleData(Obj, features = all.genes, verbose = F)
      
      Obj <- suppressWarnings(RunPCA(Obj,verbose = FALSE))
      plot(VizDimLoadings(Obj, dims = 1:2, reduction = "pca"))
      plot(DimPlot(Obj, reduction = "pca") + NoLegend())
      DimHeatmap(Obj, dims = 1, cells = 500, balanced = TRUE)
      plot(ElbowPlot(Obj))
    }
    
    if (type == "Post") {
      DimHeatmap(Obj, dims = 1:ndims, cells = 500, balanced = TRUE)
    }
    
    dev.off()
    cat(paste("[ " , my.date() , " ]"," Please check", outfile, " for select the best filter thresholds\n",sep=""))
    return(Obj)
  })
}
filter<-function(mergeSeurat,nFeature_thresh,mt_thresh, ndims){
  mergeSeurat <- subset(mergeSeurat, subset = nFeature_RNA > min(nFeature_thresh) & nFeature_RNA < max(nFeature_thresh) & percent.mt < mt_thresh)
  #UMAP
  mergeSeurat <- FindNeighbors(mergeSeurat, dims = 1:ndims, verbose = F)
  mergeSeurat <- FindClusters(mergeSeurat, resolution = 0.5, verbose = F)
  mergeSeurat <- RunUMAP(mergeSeurat, dims = 1:ndims, verbose = F)
  
  #DimPlot(mergeSeurat, reduction = "umap")
  #DimPlot(mergeSeurat, group.by = "Sample_Name")
  
  #print(mergeSeurat)
  #print(mergeSeurat.filtered)
  return(quality(mergeSeurat, type="Post", ndims=ndims))
}
cell_anottation<-function(dataset){
  cat(paste("[ " , my.date() , " ]"," Starting cell annottation\n",sep=""))
  ##4 cell anotattion
  sobj_counts <- suppressWarnings(as.SingleCellExperiment(dataset))
  sobj_counts <- suppressWarnings(logNormCounts(sobj_counts))
  ##HPCA
  suppressMessages(ref <- celldex::HumanPrimaryCellAtlasData())
  cat(paste("[ " , my.date() , " ]"," Anotatting using HPCA...\n",sep=""))
  pred <- SingleR(test=sobj_counts, ref=ref, labels=ref$label.main, BPPARAM=MulticoreParam(8))
  dataset[["GeneralCell_HPCA"]] <- pred$labels
  
  pred <- SingleR(test=sobj_counts, ref=ref, labels=ref$label.fine, BPPARAM=MulticoreParam(8))
  dataset[["PreciseCell_HPCA"]] <- pred$labels
  
  ##Monaco
  cat(paste("[ " , my.date() , " ]"," Anotatting using Monaco...\n",sep=""))
  suppressMessages(ref <- celldex::MonacoImmuneData())
  pred <- SingleR(test=sobj_counts, ref=ref, labels=ref$label.main, BPPARAM=MulticoreParam(8))
  dataset[["GeneralCell_Monaco"]] <- pred$labels
  
  pred <- SingleR(test=sobj_counts, ref=ref, labels=ref$label.fine, BPPARAM=MulticoreParam(8))
  dataset[["PreciseCell_Monaco"]] <- pred$labels
  
  #mouseRNAseq
  suppressMessages(ref <- celldex::MouseRNAseqData())
  cat(paste("[ " , my.date() , " ]"," Anotatting using MouseRNASeq\n",sep=""))
  pred <- SingleR(test=sobj_counts, ref=ref, labels=ref$label.main, BPPARAM=MulticoreParam(8))
  dataset[["GeneralCell_Mouse"]] <- pred$labels
  
  pred <- SingleR(test=sobj_counts, ref=ref, labels=ref$label.fine, BPPARAM=MulticoreParam(8))
  dataset[["PreciseCell_Mouse"]] <- pred$labels
  
  return(dataset)
}
tabula_muris<-function(Obj){
  suppressWarnings({
    
    # Cargar las librerías necesarias
    library(Seurat)       # Para trabajar con datos single-cell
    library(SeuratDisk)   # Para conversión si se necesita en un futuro
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
    timo_tabula_seurat <- RunUMAP(timo_tabula_seurat, reduction = "pca", dims = 1:20)
    
    # Realizar clustering para identificar grupos de células
    timo_tabula_seurat <- FindNeighbors(timo_tabula_seurat, dims = 1:20)
    timo_tabula_seurat <- FindClusters(timo_tabula_seurat, resolution = 0.5)
    
    # Visualizar los clusters en UMAP
    DimPlot(timo_tabula_seurat, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5) + 
      ggtitle("Tabula Muris - Thymus Data")
    
    #3 Anotar Tabula segun la web
    df <- read.csv("Anotacion/annotations_droplet_Thymus.csv", header = F)
    df$V1<-paste0(df$V1, "-1")
    timo_tabula_seurat[["cell_anot"]] <- NA
    timo_tabula_seurat$cell_anot <- df$V3[match(colnames(timo_tabula_seurat), df$V1)]
    
    # 4. Transferir anotaciones de Tabula Muris a tus datos
    # Encontrar "anclajes" entre los datos de referencia y tus datos
    anchors <- FindTransferAnchors(reference = timo_tabula_seurat, 
                                   query = Obj, 
                                   dims = 1:20)
    
    # Transferir las etiquetas de anotación desde Tabula Muris a tus datos
    predictions <- TransferData(anchorset = anchors, 
                                refdata = timo_tabula_seurat$cell_anot, 
                                dims = 1:20)
    
    # Añadir las anotaciones transferidas como metadatos a tu objeto Seurat
    Obj <- AddMetaData(Obj, metadata = predictions)
    
    # 5. Visualizar las anotaciones en tus datos
    # Generar un UMAP para ver las células anotadas
    Obj <- RunUMAP(Obj, reduction = "pca", dims = 1:20)
    Obj$cell_anot<-Obj$predicted.id
    Obj[["predicted.id"]] <- NULL
    
    return(Obj)
  })
}
find_Markers<-function(Obj,fc=0,pval=0.05, number=4, idents="seurat_clusters"){
  outfile = paste(idents, "_", deparse(substitute(Obj)), ".pdf",sep="")
  Idents(Obj)<-idents
  plan("multisession", workers = 12)
  a<-DimPlot(Obj,label = T, group.by = idents, repel = T)
    
  pdf(file = outfile, width = 14.69, height = 8.27)
  
  clusteres<-sort(unique(Obj[[idents]][,1]))
  #cat(clusteres)
  for (cluster in clusteres){
    cat(paste("[ " , my.date() , " ]"," Buscando markers in: ", idents , " [",cluster, "]\n",sep=""))
    if (sum(Obj[[idents]][,1] == cluster) > 50){
      markers<- FindMarkers(Obj, ident.1 = cluster, logfc.threshold=fc, only.pos = T, group.by = idents)
      markers<-markers[markers$p_val_adj<=pval,]
      write.table(markers,file=paste(idents, "_", cluster, "_vs_allother.tsv",sep=""),sep="\t", col.names = NA)
      genes<-rownames(markers)[0:number]
      k<-FeaturePlot(Obj,genes)
      plot(a + ggtitle(paste("Cluster", cluster, "markers")) + k)
    }else{
      cat(paste("[ " , my.date() , " ]"," Not enough cells for: ", cluster, " (only: ",sum(Obj[[idents]][,1] == cluster)," cells)\n",sep=""))
    }
  }
  dev.off()
}
RenameGenesSeurat <- function(obj, newnames) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA
  
  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    #if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
DEG_by_column<-function(Obj,column="orig.ident",label="SeuratCluster"){
  Idents(Obj)<-column
  plan("multisession", workers = 12)
  for (cell in unique(Obj[[column]][,1])) {
    cat(paste("[ " , my.date() , " ]"," Buscando markers in: ", cell , "\n",sep=""))
    if (sum(Obj[[column]][,1] == cell ) > 10){
      
      MarkersCell <- FindMarkers(Obj, ident.1 = cell)
      write.table(MarkersCell[order(MarkersCell$p_val_adj,-abs(MarkersCell$avg_log2FC)),], file = paste("DEG_", cell, "_vs_All_",label,".tsv", sep=""), sep="\t", col.names = NA, quote = F)
      
      MarkersCell <-MarkersCell[MarkersCell$p_val_adj<=0.05 & abs(MarkersCell$avg_log2FC)>0,]
      aa<-rownames(head(MarkersCell[order(MarkersCell$avg_log2FC,decreasing = T),],5))
      bb<-rownames(head(MarkersCell[order(MarkersCell$avg_log2FC,decreasing = F),],5))
      
      pdf(file = paste("DEG_", cell ,"_vs_All_",label,".pdf",sep=""), height = 11.69, width = 8.27)
      xx<-DimPlot(mergeSeurat.filter,label = T, repel = T, group.by = column)
      plot(xx)
      x<-DotPlot(Obj,features = c(aa,bb), group.by=column,col.min=0, col.max=2) + coord_flip() + ggtitle(label)
      plot(x)
      p<-FeaturePlot(Obj,c(aa[1:2],bb[1:2]), split.by = "orig.ident") + labs(caption = paste( cell, " KO/WT", sep="")) + theme(plot.caption = element_text(hjust = -0.5, face = "bold"))
      plot(p & theme(legend.position = "right"))
      dev.off()
    }
  }
}
DEG_by_column_anot<-function(Obj,column="orig.ident",anot="GeneralCell_HPCA",label="HPCA"){
  Idents(Obj)<-column
  plan("multisession", workers = 12)
  for (cell in unique(Obj[[anot]][,1])) {
    cat(paste("[ " , my.date() , " ]"," Buscando markers in: ", cell , " in [",paste0(unique(Obj$orig.ident),collapse= "-"),"]\n",sep=""))
    
    if (sum(Obj[[anot]][,1] == cell ) > 10){
      
      subset_seurat <- Obj %>%
        subset(cells = rownames(Obj@meta.data)[Obj@meta.data[[anot]] == cell])
      
      MarkersCell <- FindMarkers(subset_seurat, ident.1 = unique(Obj$orig.ident)[1], ident.2 = unique(Obj$orig.ident)[2])
      write.table(MarkersCell[order(MarkersCell$p_val_adj,-abs(MarkersCell$avg_log2FC)),], file = paste("DEG_", cell, "_vs_All_",label,".tsv", sep=""), sep="\t", col.names = NA, quote = F)
      
      MarkersCell <-MarkersCell[MarkersCell$p_val_adj<=0.05 & abs(MarkersCell$avg_log2FC)>0,]
      aa<-rownames(head(MarkersCell[order(MarkersCell$avg_log2FC,decreasing = T),],5))
      bb<-rownames(head(MarkersCell[order(MarkersCell$avg_log2FC,decreasing = F),],5))
      
      pdf(file = paste("DEG_", cell ,"_vs_",cell, paste0(unique(Obj$orig.ident),collapse= "-"), "_", label,".pdf",sep=""), height = 11.69, width = 8.27)
      #xx<-DimPlot(mergeSeurat.filter,label = T, repel = T, group.by = column)
      #plot(xx)
      #x<-DotPlot(Obj,features = c(aa,bb), group.by=column,col.min=0, col.max=2) + coord_flip() + ggtitle(label)
      #plot(x)
      p<-FeaturePlot(Obj,c(aa[1:2],bb[1:2]), split.by = "orig.ident") + labs(caption = paste( cell, " KO/WT", sep="")) + theme(plot.caption = element_text(hjust = -0.5, face = "bold"))
      plot(p & theme(legend.position = "right"))
      dev.off()
    }
  }
}
