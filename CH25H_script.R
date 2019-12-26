library(Seurat)
library(ggplot2)
library(cowplot)

#データの読み込み
file<-"GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"

#１行目と２行目が変なので３行目から読み込む
ex<-data.table::fread(file, stringsAsFactors=FALSE, sep="\t", data.table = FALSE,skip=2)
rownames(ex)<-ex[,1]
e<-ex[,-c(1,16293)]

#１行目を読み込む
L1<-readLines(con = file, n = 1)
L1<-unlist(strsplit(L1,"\t"))

#２行目を読み込む
L2<-readLines(con = file, n = 2)
L2<-unlist(strsplit(L2,"\t"))

#発現値データに列名、行名をつける
colnames(e)<-L1[-1]

#サンプル情報の読み込み
df<-read.csv("GSE120575_patient_ID_single_cells.txt",skip=19,sep="\t",as.is=T)
df<-df[,!is.na(df[1,])]

#データの整理
df<-df[1:16291,]
rownames(df)<-df[,2]
common<-intersect(colnames(e),rownames(df))
df<-df[common,]
e<-e[,common]
df[,1]<-sub(" ","",df[,1])
colnames(e)<-df[,1]
rownames(df)<-df[,1]


#CD4陽性細胞を取り出す 
ex_CD4 <- e["CD4",]
ex_cd4_posi <- e[,which(ex_CD4 > 0)]

E <- ex_cd4_posi
common<-intersect(colnames(E),rownames(df))
df<-df[common,]
e<-E[,common]
df[,1]<-sub(" ","",df[,1])
colnames(e)<-df[,1]
rownames(df)<-df[,1]


################################################
#TAM（CD4+T細胞）だけを抽出したRAWデータ。
write.csv(E,"ex_CD4.csv")
################################################


e<-read.csv("ex_CD4.csv",as.is=T,row=1)

#パッケージの読み込み
require(Seurat)

#スーラオブジェクトの作成
EX <- CreateSeuratObject(counts = e, project = "GSE120575", min.cells = 0)


#ミトコンドリア遺伝子の含有率を調べる
EX[["percent.mt"]] <- PercentageFeatureSet(EX, pattern = "^MT-")

#遺伝子数、カウント数、ミトコンドリア遺伝子含有率を可視化
VlnPlot(EX, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

#正規化ここは飛ばしていいはず
#Supplementary_files_format_and_content: tab-delimited text file containing log2(TPM+1) values (transcript-per-million reads) for 55,737 transcripts (rows) and 16,291 cells (columns)
#EX <- NormalizeData(EX, normalization.method = "LogNormalize", scale.factor = 10000)

#分散を調べる
EX <- FindVariableFeatures(EX, selection.method = "vst", nfeatures = 2000)

#スケーリングする
all.genes <- rownames(EX)
EX <- ScaleData(EX, features = all.genes)

#主成分スコアの計算
EX <- RunPCA(EX, features = VariableFeatures(object = EX))

#クラスターを探す
EX <- FindNeighbors(EX, dims = 1:20)
EX <- FindClusters(EX, resolution = 0.8)

#UMAPの計算
EX <- RunUMAP(EX, dims = 1:20)

#saveRDS(EX,file="EX.rds")






################################################
#TAM（CD4+T細胞）のクラスタリング
png("result/TAM_pos_UMAP.png")
DimPlot(EX, reduction = "umap")
dev.off()
################################################

################################################
#陽性細胞の発現値の可視化
png("result/TAM_pos_UMAP_CD4.png")
FeaturePlot(EX, features="CD4")
dev.off()

png("result/TAM_pos_UMAP_CH25H.png")
FeaturePlot(EX, features="CH25H")
dev.off()
################################################


#カウントデータ
e <- as.matrix(EX@assays$RNA@counts)

#CH25Hの情報
mygroup <- factor(ifelse(e["CH25H",] > 0, "Positive","Negative"))

EX[["mygroup"]] <- mygroup
deg <- FindMarkers(EX, ident.1="Positive" ,ident.2="Negative",group.by = "mygroup")
deg$avg_logFC <- 2^deg$avg_logFC
colnames(deg) [2] <- "FC"
write.csv(deg, "CH25H_DEG.csv")

