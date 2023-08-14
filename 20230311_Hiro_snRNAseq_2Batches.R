library(Seurat)
library(xlsx)
library(ggplot2)
library(DoubletFinder)
####ImportData####
##load CellBender output Filtered.H5 files
LH001.data <- Read10X_h5("./data/LH001/20230310_LH001_CellBender_output_filtered.h5")
LH001_seurat<- CreateSeuratObject(counts = LH001.data, project = "LH001")
LH002.data <- Read10X_h5("./data/LH002/20230310_LH002_CellBender_output_filtered.h5")
LH002_seurat <- CreateSeuratObject(counts = LH002.data, project = "LH002")
###QC and remove doublet####
##LH001
LH001_seurat$log10GenesPerUMI <- log10(LH001_seurat$nFeature_RNA)/log10(LH001_seurat$nCount_RNA)
mitogenes<-grep('^mt-',LH001_seurat@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)
LH001_seurat$mitoRatio <- PercentageFeatureSet(object = LH001_seurat, pattern = "^mt-")
LH001_seurat$mitoRatio <- LH001_seurat@meta.data$mitoRatio/100
Ribogenes<-grep('^Rp[sl][[:digit:]]',LH001_seurat@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)
LH001_seurat$RiboRatio <- PercentageFeatureSet(object = LH001_seurat, pattern = "^Rp[sl][[:digit:]]")
LH001_seurat$RiboRatio <- LH001_seurat@meta.data$RiboRatio/100
metadata <- LH001_seurat@meta.data
metadata<-dplyr::rename (metadata,seq_folder=orig.ident,
                         nUMI = nCount_RNA,
                         nGene = nFeature_RNA)
LH001_seurat@meta.data <- metadata
LH001_seurat@meta.data$sample<-LH001_seurat@meta.data$seq_folder
filtered_LH001_seurat <-subset(x = LH001_seurat, 
                               subset= ((nUMI >= 500) & 
                                          (nGene >= 400)&
                                          (log10GenesPerUMI>0.9)&
                                          (mitoRatio < 0.15)&
                                          (RiboRatio<0.05)))
filtered_LH001_seurat <- NormalizeData(filtered_LH001_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_LH001_seurat<- CellCycleScoring(filtered_LH001_seurat, 
                                 g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)), 
                                 s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),
                                 set.ident=T)
filtered_LH001_seurat <- FindVariableFeatures(filtered_LH001_seurat, selection.method = "vst", nfeatures = 2000)
filtered_LH001_seurat <- ScaleData(filtered_LH001_seurat, vars.to.regress = c("mitoRatio",
                                                              "S.Score", 
                                                              "G2M.Score",
                                                              "RiboRatio"))
filtered_LH001_seurat <- RunPCA(filtered_LH001_seurat, 
                                features = VariableFeatures(object = filtered_LH001_seurat))
DimPlot(filtered_LH001_seurat, reduction = "pca")
ElbowPlot(filtered_LH001_seurat)
filtered_LH001_seurat<-RunUMAP(filtered_LH001_seurat,dims = 1:20,reduction = "pca")
filtered_LH001_seurat <- FindNeighbors(object = filtered_LH001_seurat,dims = 1:20)
filtered_LH001_seurat <- FindClusters(object = filtered_LH001_seurat,resolution = c(0.1))
Idents(filtered_LH001_seurat)<-"RNA_snn_res.0.1"
sweep.res.list.LH001 <- paramSweep_v3(filtered_LH001_seurat, PCs = 1:50, sct = FALSE)
sweep.stats.LH001 <- summarizeSweep(sweep.res.list.LH001, GT = FALSE)
bcmvn.LH001 <- find.pK(sweep.stats.LH001)
ggplot(bcmvn.LH001,mapping = aes(pK,BCmetric,group=1))+
  geom_point()+
  geom_line()
homotypic.prop.LH001 <- modelHomotypic(filtered_LH001_seurat@meta.data$RNA_snn_res.0.1)           
nExp_poi.LH001 <- round(0.04*nrow(filtered_LH001_seurat@meta.data))  
nExp_poi.LH001.adj <- round(nExp_poi.LH001*(1-homotypic.prop.LH001))
filtered_LH001_seurat <- doubletFinder_v3(filtered_LH001_seurat, PCs = 1:50, pN = 0.25, pK = 0.005, 
                                   nExp = nExp_poi.LH001.adj, reuse.pANN = FALSE, sct = FALSE)
DimPlot(object = filtered_LH001_seurat,group.by = "DF.classifications_0.25_0.005_164")
VlnPlot(filtered_LH001_seurat, features = "nFeature_RNA", 
        group.by = "DF.classifications_0.25_0.005_164", pt.size = 0.1)
filtered_LH001_seurat_singlet<-filtered_LH001_seurat[, filtered_LH001_seurat@meta.data[, "DF.classifications_0.25_0.005_164"] == "Singlet"]
##LH002
LH002_seurat$log10GenesPerUMI <- log10(LH002_seurat$nFeature_RNA)/log10(LH002_seurat$nCount_RNA)
mitogenes<-grep('^mt-',LH002_seurat@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)
LH002_seurat$mitoRatio <- PercentageFeatureSet(object = LH002_seurat, pattern = "^mt-")
LH002_seurat$mitoRatio <- LH002_seurat@meta.data$mitoRatio/100
Ribogenes<-grep('^Rp[sl][[:digit:]]',LH002_seurat@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)
LH002_seurat$RiboRatio <- PercentageFeatureSet(object = LH002_seurat, pattern = "^Rp[sl][[:digit:]]")
LH002_seurat$RiboRatio <- LH002_seurat@meta.data$RiboRatio/100
metadata <- LH002_seurat@meta.data
metadata<-dplyr::rename (metadata,seq_folder=orig.ident,
                         nUMI = nCount_RNA,
                         nGene = nFeature_RNA)
LH002_seurat@meta.data <- metadata
LH002_seurat@meta.data$sample<-LH002_seurat@meta.data$seq_folder
filtered_LH002_seurat <-subset(x = LH002_seurat, 
                               subset= ((nUMI >= 500) & 
                                          (nGene >= 400)&
                                          (log10GenesPerUMI>0.9)&
                                          (mitoRatio < 0.15)&
                                          (RiboRatio<0.05)))
filtered_LH002_seurat <- NormalizeData(filtered_LH002_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
filtered_LH002_seurat<- CellCycleScoring(filtered_LH002_seurat, 
                                         g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)), 
                                         s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),
                                         set.ident=T)
filtered_LH002_seurat <- FindVariableFeatures(filtered_LH002_seurat, selection.method = "vst", nfeatures = 2000)
filtered_LH002_seurat <- ScaleData(filtered_LH002_seurat, vars.to.regress = c("mitoRatio",
                                                                              "S.Score", 
                                                                              "G2M.Score",
                                                                              "RiboRatio"))
filtered_LH002_seurat <- RunPCA(filtered_LH002_seurat, 
                                features = VariableFeatures(object = filtered_LH002_seurat))
DimPlot(filtered_LH002_seurat, reduction = "pca")
ElbowPlot(filtered_LH002_seurat)
filtered_LH002_seurat<-RunUMAP(filtered_LH002_seurat,dims = 1:20,reduction = "pca")
filtered_LH002_seurat <- FindNeighbors(object = filtered_LH002_seurat,dims = 1:20)
filtered_LH002_seurat <- FindClusters(object = filtered_LH002_seurat,resolution = c(0.1))
Idents(filtered_LH002_seurat)<-"RNA_snn_res.0.1"
sweep.res.list.LH002 <- paramSweep_v3(filtered_LH002_seurat, PCs = 1:50, sct = FALSE)
sweep.stats.LH002 <- summarizeSweep(sweep.res.list.LH002, GT = FALSE)
bcmvn.LH002 <- find.pK(sweep.stats.LH002)
ggplot(bcmvn.LH002,mapping = aes(pK,BCmetric,group=1))+
  geom_point()+
  geom_line()
homotypic.prop.LH002 <- modelHomotypic(filtered_LH002_seurat@meta.data$RNA_snn_res.0.1)     
nExp_poi.LH002 <- round(0.04*nrow(filtered_LH002_seurat@meta.data)) 
nExp_poi.LH002.adj <- round(nExp_poi.LH002*(1-homotypic.prop.LH002))
filtered_LH002_seurat <- doubletFinder_v3(filtered_LH002_seurat, PCs = 1:50, pN = 0.25, pK = 0.14, 
                                          nExp = nExp_poi.LH002.adj, reuse.pANN = FALSE, sct = FALSE)
DimPlot(object = filtered_LH002_seurat,group.by = "DF.classifications_0.25_0.14_147")
VlnPlot(filtered_LH002_seurat, features = "nFeature_RNA", 
        group.by = "DF.classifications_0.25_0.14_147", pt.size = 0.1)
filtered_LH002_seurat_singlet<-filtered_LH002_seurat[, filtered_LH002_seurat@meta.data[, "DF.classifications_0.25_0.14_147"] == "Singlet"]
###Integrate two batches####
combined_seurat<-merge(filtered_LH001_seurat_singlet,filtered_LH002_seurat_singlet)
split_seurat <- SplitObject(combined_seurat, split.by = "sample")
split_seurat <- split_seurat[c("LH001", "LH002")]
i<-1
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], 
                                        g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)), 
                                        s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),
                                        set.ident=T)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   vars.to.regress = c("mitoRatio",
                                                       "S.Score", 
                                                       "G2M.Score",
                                                       "RiboRatio"),
                                   variable.features.n=3000)
}
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features) 
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
integ_anchors1<-integ_anchors
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
seurat_integrated <- RunPCA(object = seurat_integrated,
                            ndims.print = 1:50,nfeatures.print = 5)
PCAPlot(seurat_integrated) 
ElbowPlot(seurat_integrated,ndims = 50)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:50,
                             reduction = "pca")
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:50)
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.1))
DefaultAssay(seurat_integrated)<-"RNA"
Idents(seurat_integrated)<-"integrated_snn_res.0.1"
seurat_integrated <- RenameIdents(object = seurat_integrated,
                        "0" = "LepR+",
                        "1" = "Hematopoietic",
                        "2" = "Hematopoietic",
                        "3" = "Endothelial cell",
                        "4" = "Adipocyte",
                        "5" = "Hematopoietic",
                        "6" = "Hematopoietic",
                        "7" = "Pericyte")
###Extract and integrate adipocyte and LepR cells from two batches####
stromal_seurat<-subset(seurat_integrated, ident = c("LepR+","Adipocyte"))
stromal_seurat@meta.data$batch<-"0"
stromal_seurat@meta.data[stromal_seurat@meta.data$sample=="LH001",]$batch<-"batch1"
stromal_seurat@meta.data[stromal_seurat@meta.data$sample=="LH002",]$batch<-"batch2"
split_seurat0<- SplitObject(stromal_seurat, split.by = "batch")
split_seurat0<- split_seurat0[c("batch1", "batch2")]
i<-1
for (i in 1:length(split_seurat0)) {
  split_seurat0[[i]] <-NormalizeData(object =split_seurat0[[i]] )
  split_seurat0[[i]] <- CellCycleScoring(split_seurat0[[i]], 
                                         g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)), 
                                         s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),
                                         set.ident=T)
  split_seurat0[[i]] <- SCTransform(split_seurat0[[i]], 
                                    vars.to.regress = c("S.Score", 
                                                        "G2M.Score"),
                                    variable.features.n=3000)
}
integ_features0<- SelectIntegrationFeatures(object.list = split_seurat0, 
                                            nfeatures = 3000) 
split_seurat0<- PrepSCTIntegration(object.list = split_seurat0, 
                                   anchor.features = integ_features0) 
integ_anchors0 <- FindIntegrationAnchors(object.list = split_seurat0, 
                                         normalization.method = "SCT", 
                                         anchor.features = integ_features0)
seurat_integrated0<- IntegrateData(anchorset = integ_anchors0, 
                                   normalization.method = "SCT")
seurat_integrated0<- RunPCA(object = seurat_integrated0,
                            ndims.print = 1:50,nfeatures.print = 5)
PCAPlot(seurat_integrated0) 
ElbowPlot(seurat_integrated0,ndims = 50)
seurat_integrated0<- RunUMAP(seurat_integrated0, 
                             dims = 1:50,
                             reduction = "pca")

seurat_integrated0<- FindNeighbors(object = seurat_integrated0, 
                                   dims = 1:50)
seurat_integrated0<- FindClusters(object = seurat_integrated0,
                                  resolution = c(0.2))
ii<-1
DefaultAssay(seurat_integrated0)<-"RNA"
Idents(seurat_integrated0)<-"integrated_snn_res.0.2"
seurat_integrated0<-RenameIdents(seurat_integrated0,
                                 "0"="LepR+_Adipo",
                                 "1"="LepR+_Adipo",
                                 "2"="Adipocyte",
                                 "3"="LepR+_Osteo")
seurat_integrated0@meta.data$celltype<-Idents(seurat_integrated0)
seurat_integrated0@meta.data$celltype<-factor(seurat_integrated0@meta.data$celltype,
                                              levels=c("LepR+_Osteo","LepR+_Adipo","Adipocyte"))
UmapAllCells<-DimPlot(object = seurat_integrated0,
                      pt.size = 1,label=T)
ggsave(filename = "./result/20230311_LeprAdipocyte_UmapAllCells.pdf",
       plot =UmapAllCells,
       device="pdf",width = 9,height=8)
DefaultAssay(seurat_integrated0) <- "RNA"
Cxcl12AllCells<-FeaturePlot(object = seurat_integrated0,features = "Cxcl12",
                            label = T,min.cutoff = "q10",
                            pt.size = 1,
                            cols = c("grey","red"))
ggsave(filename = "./result/20230311_LeprAdipocyte_Cxcl12AllCells.pdf",
       plot =Cxcl12AllCells,
       device="pdf",width = 8,height=8)
Fabp4AllCells<-FeaturePlot(object = seurat_integrated0,features = "Fabp4",
                           label = T,min.cutoff = "q10",
                           max.cutoff = "q80",pt.size = 1,
                           cols = c("grey","red"))
ggsave(filename = "./result/20230311_LeprAdipocyte_Fabp4AllCells.pdf",
       plot =Fabp4AllCells,
       device="pdf",width = 8,height=8)
Idents(seurat_integrated0)<-factor(x = Idents(seurat_integrated0),
                                   levels = c("LepR+_Osteo","LepR+_Adipo","Adipocyte"))
Vln0<-VlnPlot(object = seurat_integrated0,
              features = c("Lepr","Cxcl12","Kitl","Fabp4","Adipoq",
                                                 "Lpl","Sp7","Bglap","Alpl","Postn","Spp1",
                                                 "Plin1","Plin2","Plin3","Plin4","Plin5",
                                                 "Pnpla2", "Lipe", "Mgll","Tnc",
                                                 "Adrb1", "Adrb2", "Adrb3",
                                                 "Il7r","Il17a","Il6ra"),
              assay = "RNA",cols = c("#619cff","#f8766d","#00ba38"))
ggsave(filename = "./result/20230311_LeprAdipocyte_Vlnplot_AllCells.pdf",plot = Vln0,
       device="pdf",width = 10,height=20)
Dotplot0_1<-DotPlot(seurat_integrated0,features = c("Lepr", "Cxcl12",
                                                  "Kitl", "Il7", "Fabp4", 
                                                  "Lpl", "Adipoq"),
                  assay = "RNA") + scale_size(range = c(1,10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "./result/20230311_LeprAdipocyte_Dotplot_AllCells_1.pdf",plot = Dotplot0_1,
       device="pdf",width = 6,height=3)
Dotplot0_2<-DotPlot(seurat_integrated0,features = c("Alpl", "Spp1", "Sp7", "Postn", "Bglap"),
                    assay = "RNA") + scale_size(range = c(1,10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "./result/20230311_LeprAdipocyte_Dotplot_AllCells_2.pdf",plot = Dotplot0_2,
       device="pdf",width = 6,height=3)
Dotplot0_3<-DotPlot(seurat_integrated0,features = c("Plin1","Plin2","Plin3","Plin4","Plin5",
                                                    "Pnpla2", "Lipe", "Mgll"),
                    assay = "RNA") + scale_size(range = c(1,10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "./result/20230311_LeprAdipocyte_Dotplot_AllCells_3.pdf",plot = Dotplot0_3,
       device="pdf",width = 6,height=3)
Dotplot0_4<-DotPlot(seurat_integrated0,features = c("Adrb1", "Adrb2", "Adrb3",
                                                    "Il6ra"),
                    assay = "RNA") + scale_size(range = c(1,10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(filename = "./result/20230311_LeprAdipocyte_Dotplot_AllCells_4.pdf",plot = Dotplot0_4,
       device="pdf",width = 6,height=3)
pdf(file ="./result/20230311_LeprAdipocyte_pheatmap_BYScillus.pdf",
    width =10,height = 10)
Scillus::plot_heatmap(dataset = ScaleData(seurat_integrated0,use_raster=F),
                      n = 12,
                      markers =FindAllMarkers(seurat_integrated0, logfc.threshold = 0.1, min.pct = 0, only.pos = T),
                      sort_var = c("celltype"),
                      anno_var = c("celltype"),
                      anno_colors = list(c("#619cff","#f8766d","#00ba38")),
                      hm_limit = c(-1,0,1),
                      hm_colors = c("purple","black","yellow")
)
dev.off()
###Extract tdtomato expression####
GENES<-rownames(seurat_integrated0)
CELLS<-colnames(seurat_integrated0)
indMyGene <- which(GENES %in% c("tdWPRE","Lepr"))
CountMyGeneMyCell<-GetAssayData(object = seurat_integrated0, slot = 'counts')[indMyGene, ]
CountMyGeneMyCell==0
tdWPRE_zero_Cells<-colnames(CountMyGeneMyCell)[CountMyGeneMyCell[2,]==0]
tdWPRE_nonzero_Cells<-colnames(CountMyGeneMyCell)[CountMyGeneMyCell[2,]!=0]
seurat_integrated0@meta.data$tdt<-"NA"
seurat_integrated0@meta.data[tdWPRE_zero_Cells,]$tdt<-"zero"
seurat_integrated0@meta.data[tdWPRE_nonzero_Cells,]$tdt<-"nonzero"
####Cell type proportion####
library("aod")
Data1<-seurat_integrated0@meta.data
for (iii in c("LepR+_Osteo","LepR+_Adipo","Adipocyte")){
  Data1$CluterorNOT<-0
  Data1[Data1$celltype==iii,]$CluterorNOT<-1
  tdt_zero<-binom::binom.confint(x = nrow(Data1[Data1$tdt=="zero"&Data1$CluterorNOT=="1",]),
                           n =nrow(Data1[Data1$tdt=="zero",]),methods = "logit")
  tdt_nonzero<-binom::binom.confint(x = nrow(Data1[Data1$tdt=="nonzero"&Data1$CluterorNOT=="1",]),
                            n =nrow(Data1[Data1$tdt=="nonzero",]),methods = "logit")
  write.csv(x = tdt_zero,paste("./result/20230202_Binom_CI_",iii,"_tdt_zero.csv",sep=""))
  write.csv(x = tdt_nonzero,paste("./result/20230202_Binom_CI_",iii,"_tdt_nonzero.csv",sep=""))
}
#https://github.com/rpolicastro/scProportionTest
library("scProportionTest")
seurat_data<-seurat_integrated0
prop_test<-sc_utils(seurat_data)
prop_test@meta_data$celltype<-factor(prop_test@meta_data$celltype,
                                     levels=c("LepR+_Osteo","LepR+_Adipo","Adipocyte"))
prop_test <- permutation_test(
  prop_test, cluster_identity = "celltype",
  sample_1 = "zero", sample_2 = "nonzero",
  sample_identity = "tdt"
)
ClusterDiff<-permutation_plot(prop_test)
write.csv(x = prop_test@results$permutation,file = "./result/20230311_LeprAdipocyte_ClusterDiff.csv")
##Divide into tdt- and tdt+ cells####
GENES0<-rownames(seurat_integrated0)
CELLS0<-colnames(seurat_integrated0)
names(seurat_integrated0)
dim(x = seurat_integrated0)
indMyGene <- which(GENES0 %in% c("tdWPRE","Lepr"))
seurat_integrated0_tdWPRE_zero<-subset(seurat_integrated0, subset=(tdt=="zero"))
seurat_integrated0_tdWPRE_nonzero<-subset(seurat_integrated0, subset=(tdt=="nonzero"))
###tdt+ cells###
DefaultAssay(seurat_integrated0_tdWPRE_nonzero)<-"RNA"
seurat_integrated0_tdWPRE_nonzero$celltype<-factor(seurat_integrated0_tdWPRE_nonzero$celltype,
                                                   levels=c("LepR+_Osteo","LepR+_Adipo","Adipocyte"))
UmapAllCells_tdWPRE_nonzero<-DimPlot(object = seurat_integrated0_tdWPRE_nonzero,group.by = "celltype",
                                     pt.size = 1,label=T,cols=c("#619cff","#f8766d","#00ba38"))
ggsave(filename = "./result/20230311_LeprAdipocyte_tdWPRE_nonzero_UmapAllCells.pdf",
       plot =UmapAllCells_tdWPRE_nonzero,
       device="pdf",width = 9,height=8)
DefaultAssay(seurat_integrated0_tdWPRE_nonzero) <- "RNA"
Cxcl12AllCells_tdWPRE_nonzero<-FeaturePlot(object = seurat_integrated0_tdWPRE_nonzero,features = "Cxcl12",
                                           label = T,min.cutoff = "q10",
                                           pt.size = 1,
                                           cols = c("grey","red"))
ggsave(filename = "./result/20230311_LeprAdipocyte_tdWPRE_nonzero_Cxcl12AllCells.pdf",
       plot =Cxcl12AllCells_tdWPRE_nonzero,
       device="pdf",width = 8,height=8)
Fabp4AllCells_tdWPRE_nonzero<-FeaturePlot(object = seurat_integrated0_tdWPRE_nonzero,features = "Fabp4",
                                          label = T,min.cutoff = "q10",
                                          max.cutoff = "q80",pt.size = 1,
                                          cols = c("grey","red"))
ggsave(filename = "./result/20230311_LeprAdipocyte_tdWPRE_nonzero_Fabp4AllCells.pdf",
       plot =Fabp4AllCells_tdWPRE_nonzero,
       device="pdf",width = 8,height=8)

