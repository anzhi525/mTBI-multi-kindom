
##2.1 GMM与菌---------------------

#输入分组
metadata = read.table (file = "metadataFBT.txt", row.names = 1, header = T, sep = "\t")


GMM_abundance = read.table (file = "diff_GMM_GBM.txt", sep = "\t", row.names = 1, header = T)
GMM_abundance <- GMM_abundance[which(rowSums(GMM_abundance) > 0),]
GMM_abundance <- GMM_abundance[,rownames(metadata)]
GMM_abundance <-as.data.frame( t(GMM_abundance))



#输入菌表
metabolomic = read.table (file = "diff_species.txt", sep = "\t", row.names = 1, header = T)
metabolomic <- metabolomic[which(rowSums(metabolomic) > 0),]
metabolomic <- as.data.frame(t(metabolomic))
metabolomic <- metabolomic[rownames(GMM_abundance),]


res <-  corr.test(metabolomic,GMM_abundance,use='pairwise',method='spearman',adjust='BH',alpha=0.05)
names(res)#查看列表每个项目

res$stars

#
# 保存需要的文件
library(reshape2)
phenotyper <- res$r#提取r数据
phenotyper <- na.omit(phenotyper)

phenotypep <- res$p#提取r数据
phenotypep <- na.omit(phenotypep)

phenotypepadj <- res$p.adj#提取r数据
phenotypepadj <- na.omit(phenotypepadj)

statistic <-res$stars
statistic <- na.omit(statistic)
#class(res)
phenotyperlong <- melt(phenotyper, id.vars = colnames(phenotyper))
colnames(phenotyperlong) <- c("species","GMM","Correlation")
# View(phenotyperlong)

phenotypeplong <- melt(phenotypep, id.vars = colnames(phenotypep))
colnames(phenotypeplong) <- c("species","GMM","pvalue")
# View(phenotypeplong)


phenotypepadjlong <- melt(phenotypepadj, id.vars = colnames(phenotypepadj))
colnames(phenotypepadjlong) <- c("species","GMM","padjust")


phenotypestatisticlong <- melt(statistic, id.vars = colnames(statistic))
colnames(phenotypestatisticlong) <- c("species","GMM","statistic")


hebingphenotyperp <- merge(phenotyperlong, phenotypeplong , by = c("species","GMM"))
hebingphenotyperpj <- merge(hebingphenotyperp,phenotypepadjlong, by = c("species","GMM"))
hebingphenotyperpj <- merge(hebingphenotyperpj,phenotypestatisticlong, by = c("species","GMM"))
# View(hebingphenotyperp)


library(openxlsx)
write.xlsx(hebingphenotyperpj, paste0("GMM/species_GMM_corr_un.all.xlsx"),
           append = FALSE,
           colNames = TRUE,
           rowNames = FALSE)
write.table(hebingphenotyperpj,file = "GMM/species_GMM_corr_un.all.txt",sep = "\t",quote = F,row.names = F)

data_merge1 <- subset(hebingphenotyperpj, hebingphenotyperpj$pvalue <= 0.05)
data_merge_filtered <- subset(data_merge1, abs(data_merge1$Correlation) >= 0.5)


#保存文件（筛选的）
write.xlsx(data_merge_filtered, paste0("GMM/species_GMM_corr_filter.all.xlsx"),
           append = FALSE,
           colNames = TRUE,
           rowNames = FALSE)

write.table(data_merge_filtered,file = "GMM/species_GMM_corr_filter.all.txt",sep = "\t",quote = F,row.names = F)

#画图


library(readxl)


write.table(phenotyper, file=paste0("GMM/heatmapGMM.cor.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)
write.table(phenotypepadj, file=paste0("GMM/heatmapGMM.padj.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)


# write.table(plotmat, file=paste0("heatmap.cor.14.26.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)
# write.table(stars, file=paste0("heatmap.padj.14.26.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)




plotmat = read.table (file = "GMM/heatmapGMM.cor.txt", row.names = 1, header = T, sep = "\t")
stars = read.table (file = "GMM/heatmapGMM.padj.txt", row.names = 1, header = T, sep = "\t")

rows_to_delete <- c("Streptococcus_sp_263_SSPC",
                    "Streptococcus_cristatus",
                    "Rothia_dentocariosa",
                    "Isoptericola_variabilis",
                    "Klebsiella_aerogenes",
                    "Lawsonibacter_sp_NSJ_51",
                    "Cronobacter_sakazakii",
                    "Pseudocitrobacter_faecalis",
                    "Intestinibacter_bartlettii",
                    "Saccharomyces_cerevisiae",
                    "Parascardovia_denticolens",
                    "Lachnoanaerobaculum_sp_ICM7",
                    "Parvimonas_micra",
                    "Christensenella_minuta",
                    "Eubacterium_sp_AM28_29",
                    "Streptococcus_mutans",
                    "Eubacterium_sp_NSJ_61",
                    "Candidatus_Avimonas_narfia",
                    "Clostridium_sp_Marseille_P3244",
                    "Peptoniphilus_sp_Marseille_P3761",
                    "Actinomyces_oricola",
                    "Candidatus_Pararuminococcus_gallinarum",
                    "Campylobacter_hominis",
                    "Roseburia_hominis",
                    "Barnesiella_intestinihominis")


plotmat <- plotmat[!(rownames(plotmat) %in% rows_to_delete), ]
stars <- stars[!(rownames(stars) %in% rows_to_delete), ]



if (!is.null(stars)){
  ssmt <- stars< 0.01
  stars[ssmt] <-'**'
  smt <- stars >0.01& stars <0.05
  stars[smt] <- '*'
  stars[!ssmt&!smt]<- ''
} else {
  stars <- F
}


#去r的最大值与最小值，并且保留4位数然后*100，且将负值取绝对值
maxr <- round(max(plotmat),4)*100
minr <- abs(round(min(plotmat),4)*100)

library("stringi")
#设置颜色
#EE0000FF", high = "#3B4992FF"
#设置颜色，-1到-0.3是蓝到白，-0.3-0.3是白，0.3-1是从白到红
mycol <- c(colorRampPalette(c('#3B4992FF','white'))(minr),
           colorRampPalette(c('white','white'))(3),
           colorRampPalette(c('white','#EE0000FF'))(maxr))
#max.rownames0 <- max(stri_length(rownames(plotmat)))
#max.colnames0 <- max(stri_length(colnames(plotmat)))

library(pheatmap)
library(ComplexHeatmap)
#画所有的图（不筛选r的），并且显示所有p小于0.05的*


annotation = as.data.frame(read.table(paste0("pheatmap.GMM_GBM.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))

#替换矩阵的名字
plotmat <- plotmat[,rownames(annotation)]
colnames(plotmat) <- annotation$Level1

stars <- stars[,rownames(annotation)]
colnames(stars) <- annotation$Level1

annotation2 = as.data.frame(read.table(paste0("pheatmap.species.ann.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation2 <- annotation2[rownames(plotmat),]


library(RColorBrewer)


#去r的最大值与最小值，并且保留4位数然后*100，且将负值取绝对值
maxr <- round(max(plotmat),4)
minr <- round(min(plotmat),4)
#?colorRamp2()
library(circlize)
# col_main <- colorRamp2(
#   c(minr, 0, maxr), 
#   c(  "#3B4992FF","white","#EE0000FF")
# )


annotation_col <- annotation[c("Modules","Enrichment")]
annotation_raw <- annotation2[c("Family","Enrichment")]


#填充颜色设置
# library(RColorBrewer)
# getPalette0 = colorRampPalette(brewer.pal(14, "Set3"))
# getPalette0(14)
# level1_col  <- getPalette0(14)
# level1_col  <-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")
# names(level1_col) <- unique(annotation$Level1)



Modules_col <-c("#FED9A6","#FB9A99")
names(Modules_col) <- unique(annotation$Modules)


Enrichment_col<- c("#1B9E77", "#D95F02")
names(Enrichment_col) = c("FBT" , "Control" )


#设置pathway的注释
pathway_ann <- HeatmapAnnotation(
  # Level1 = annotation$Level1,
  Modules = annotation$Modules,
  Enrichment =annotation$Enrichment,
  col = list(
    # Level1 = level1_col ,
    Modules = Modules_col,
    Enrichment = Enrichment_col #分类
  )
)

# kindom_col <-c("#35B779","#EA0F6B")
# names(kindom_col) <- unique(annotation2$Kingdom)


#设置phylum颜色
# getPalette_phylum = colorRampPalette(brewer.pal(5, "Set3"))
# Phylum_col <- getPalette_phylum(5)
Phylum_col <- c( "#44598e", "#64a0c0", "#7ec4b7", "#deebcd")
names(Phylum_col) <- unique(annotation2$Phylum)


#设置family颜色
getPalette_family = colorRampPalette(brewer.pal(19, "Set2"))
Family_col <- getPalette_family(19)#Family_col <- getPalette_family(10)
names(Family_col) <- unique(annotation2$Family)



#设置菌的注释
species_ann <- rowAnnotation(
  # Kingdom= annotation2$Kingdom,
  Phylum = annotation2$Phylum,
  Family = annotation2$Family,
  Enrichment =annotation2$Enrichment,
  col = list( 
    # Kingdom=kindom_col,
    Phylum = Phylum_col , 
    Family = Family_col,
    Enrichment = Enrichment_col #分类
  )
)



spe <- rownames(plotmat)
PA <- colnames(plotmat)
column_order <- PA
row_order <- spe

cell_fun <-function(j, i, x, y, width, height, fill) {
  grid.text(stars[i,j], 
            x = x, 
            y =  y-height*0.3,#把y坐标向下移一点
            gp = gpar(fontsize = 15,col="black"))
}

mycol <- c("#93CBE1","#C8EAE0","#FDE6EC","#FEBCCC","#45A8CB")
p=Heatmap(plotmat,cluster_rows = T,cluster_columns = T,
          #设置方格长宽
          #heatmap_width和heatmap_height控制整个热图的大小（包括图例），width和height只控制热图主体的大小
          width = unit(3, "cm"),
          height = unit(17, "cm"),
          #width = 12,height = 11,
          #设置行列name
          row_names_gp = gpar(fontsize = 11,col ="black"),
          column_names_gp = gpar(fontsize = 11,col ="black"),
          row_names_max_width = unit(100,"mm"),##行名与图例距离
          heatmap_legend_param = list(
            at = c(-0.4,-0.2, 0, 0.2,0.4)
          ),
          cell_fun = cell_fun,
          #设置格子分隔
          rect_gp = gpar(col = "grey", lwd = 1),
          #图例name
          name = "SCC bg.adj",
          top_annotation = pathway_ann ,
          #bottom_annotation = column_ha, #对应的注释
          left_annotation = species_ann,
          
          row_order=row_order,
          column_order=column_order,
          #设置颜色
          col = mycol)#,
p

row_order(p)
column_order(p)


# ggsave(file ="pathway_jun.pdf", p, width = 250, height = 400, units="mm")

pdf("GMM/GMM_species2.pdf",  width=10, height=15)
p
dev.off()


png("GMM/GMM_species2.png",  width=700, height=1000,units = "px")
p
dev.off()

write.table(names(column_order(p)),file = "GMM/GMM_order.txt",quote = F,row.names = F)
 write.table(names(row_order(p)),file = "GMM/species_order.txt",quote = F,row.names = F)


##2.2 species与短酸---------------------
metadata = read.table (file = "metadataFBT.txt", row.names = 1, header = T, sep = "\t")

SCFA = read.table (file = "SCFA.txt", sep = "\t", row.names = 1, header = T, comment.char="",check.names=F,quot="", stringsAsFactors = F)
SCFA <- SCFA[rownames(metadata),]


#输入菌表

species = read.table(paste0("diff_species.txt"), header=T, row.names=1, sep="\t", comment.char="")
species <- species[,rownames(SCFA)]
species_filtered1 <- t(species)


res <-  corr.test(species_filtered1,SCFA,use='pairwise',method='spearman',adjust='BH',alpha=0.05)
names(res)#查看列表每个项目

res$stars

#
# 保存需要的文件
library(reshape2)
phenotyper <- res$r#提取r数据
phenotyper <- na.omit(phenotyper)

phenotypep <- res$p#提取r数据
phenotypep <- na.omit(phenotypep)

phenotypepadj <- res$p.adj#提取r数据
phenotypepadj <- na.omit(phenotypepadj)

statistic <-res$stars
statistic <- na.omit(statistic)
#class(res)
phenotyperlong <- melt(phenotyper, id.vars = colnames(phenotyper))
colnames(phenotyperlong) <- c("species","SCFA","Correlation")
# View(phenotyperlong)

phenotypeplong <- melt(phenotypep, id.vars = colnames(phenotypep))
colnames(phenotypeplong) <- c("species","SCFA","pvalue")
# View(phenotypeplong)


phenotypepadjlong <- melt(phenotypepadj, id.vars = colnames(phenotypepadj))
colnames(phenotypepadjlong) <- c("species","SCFA","padjust")


phenotypestatisticlong <- melt(statistic, id.vars = colnames(statistic))
colnames(phenotypestatisticlong) <- c("species","SCFA","statistic")


hebingphenotyperp <- merge(phenotyperlong, phenotypeplong , by = c("species","SCFA"))
hebingphenotyperpj <- merge(hebingphenotyperp,phenotypepadjlong, by = c("species","SCFA"))
hebingphenotyperpj <- merge(hebingphenotyperpj,phenotypestatisticlong, by = c("species","SCFA"))
# View(hebingphenotyperp)


library(openxlsx)
write.xlsx(hebingphenotyperpj, paste0("GMM/species_SCFA_corr_un.all.xlsx"),
           append = FALSE)
write.table(hebingphenotyperpj,file = "GMM/species_SCFA_corr_un.all.txt",sep = "\t",quote = F,row.names = F)

data_merge1 <- subset(hebingphenotyperpj, hebingphenotyperpj$pvalue <= 0.05)
data_merge_filtered <- subset(data_merge1, abs(data_merge1$Correlation) >= 0.5)


#保存文件（筛选的）
write.xlsx(data_merge_filtered, paste0("GMM/species_SCFA_corr_filter.all.xlsx"),
           append = FALSE)

write.table(data_merge_filtered,file = "GMM/species_SCFA_corr_filter.all.txt",sep = "\t",quote = F,row.names = F)

#画图


library(readxl)
library(circlize)

species_SCFA = read.table (file = "GMM/species_SCFA_corr_un.all.txt",  header = T, sep = "\t")
library("tidyr")
library(reshape2)
#纵轴(画图为横轴)  ~ 横轴 , value.var = "pvalues"，中间填充的值，如果缺失则为NA
#acast为将“纵轴设置为行名，而dcast行名为1,2,3,4
species_SCFA_p   <- acast(species_SCFA, species~SCFA , value.var = "pvalue")
species_SCFA_cor <- acast(species_SCFA, species~SCFA , value.var = "Correlation")


##删除没有相关性的细菌

rows_to_delete <- c("Streptococcus_sp_263_SSPC",
                    "Streptococcus_cristatus",
                    "Rothia_dentocariosa",
                    "Isoptericola_variabilis",
                    "Klebsiella_aerogenes",
                    "Lawsonibacter_sp_NSJ_51",
                    "Cronobacter_sakazakii",
                    "Pseudocitrobacter_faecalis",
                    "Intestinibacter_bartlettii",
                    "Saccharomyces_cerevisiae",
                    "Parascardovia_denticolens",
                    "Lachnoanaerobaculum_sp_ICM7",
                    "Parvimonas_micra",
                    "Christensenella_minuta",
                    "Eubacterium_sp_AM28_29",
                    "Streptococcus_mutans",
                    "Eubacterium_sp_NSJ_61",
                    "Candidatus_Avimonas_narfia",
                    "Clostridium_sp_Marseille_P3244",
                    "Peptoniphilus_sp_Marseille_P3761",
                    "Actinomyces_oricola",
                    "Candidatus_Pararuminococcus_gallinarum",
                    "Campylobacter_hominis",
                    "Roseburia_hominis",
                    "Barnesiella_intestinihominis")


species_SCFA_p <- species_SCFA_p[!(rownames(species_SCFA_p) %in% rows_to_delete), ]
species_SCFA_cor <- species_SCFA_cor[!(rownames(species_SCFA_cor) %in% rows_to_delete), ]


maxr1 <- round(max(species_SCFA_cor),4)
minr1 <- round(min(species_SCFA_cor),4)
col_linchuang_jun <- colorRamp2(
  c(minr1, 0, maxr1), 
  c(  "#FFFFB3","white","#BEBADA")
)


#rm (plotmat_p)
if (!is.null(species_SCFA_p)){
  ssmt <- species_SCFA_p< 0.01
  species_SCFA_p[ssmt] <-'**'
  smt <- species_SCFA_p >0.01& species_SCFA_p <0.05
  species_SCFA_p[smt] <- '*'
  species_SCFA_p[!ssmt&!smt]<- ''
} else {
  v <- F
}



cell_fun1 <-function(j, i, x, y, width, height, fill) {
  grid.text(species_SCFA_p[i,j], 
            x = x, y = y-height*0.3,
            gp = gpar(fontsize = 15,col="black"))}

# a <- rownames(linchuang_jun_cor)#菌
sc <- read.table("Pathway/SCFA_order.txt", header = T)
spe <- read.table("GMM/species_order.txt", header = T)
row_order<- spe$x
column_order<- sc$x

p1=Heatmap(species_SCFA_cor,
           name = "species_clinical cor",
           #设置单元格内的星号
           cell_fun = cell_fun1,
           #聚类
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           #设置行列name
           row_names_gp = gpar(fontsize = 11,col ="black"),
           column_names_gp = gpar(fontsize = 11,col ="black"),
           row_names_max_width = unit(100,"mm"),##行名与图例距离
           #设置顺序
           row_order=row_order,
           column_order =column_order, 
           #column_order=row_order(p),
           #设置方格长宽
           #heatmap_width和heatmap_height控制整个热图的大小（包括图例），width和height只控制热图主体的大小
           width = unit(2, "cm"),
           height = unit(17, "cm"),
           #width = 12,height = 11,
           #width = unit(2.28, "cm"),
           #设置格子分隔
           rect_gp = gpar(col = "grey", lwd = 1),
           #设置行不聚类
           #cluster_cols = T)#,
           #设置颜色
           col = col_linchuang_jun)#,
p1
row_order(p1)
column_order(p1)

pdf("GMM/specie_SCFA.2.pdf",  width=5, height=10)
p1
dev.off()


png("GMM/specie_SCFA.2.png",  width=500, height=700,units = "px")
p1
dev.off()

write.table(names(row_order(p1)),file = "GMM/SCFA_order.txt",quote = F,row.names = F)




##2.3 GMM与短酸---------------------


# .计算相关性-
library(psych)

metadata = read.table (file = "metadataFBT.txt", row.names = 1, header = T, sep = "\t")

SCFA = read.table (file = "SCFA.txt", sep = "\t", row.names = 1, header = T, comment.char="",check.names=F,quot="", stringsAsFactors = F)
SCFA <- SCFA[rownames(metadata),]

#输入pathway
pathway = as.data.frame(t(read.table("diff_GMM_GBM.txt",header = T,row.names = 1, sep = "\t", comment.char="",check.names=F,quot="", stringsAsFactors = F)))

pathway <- pathway[,which(colSums(pathway) > 0)]

pathway <- pathway[rownames(SCFA),]


res <-  corr.test(pathway,SCFA,use='pairwise',method='spearman',adjust='BH',alpha=0.05)
names(res)#查看列表每个项目

res$stars

#
#保存需要的文件
library(reshape2)
phenotyper <- res$r#提取r数据
phenotyper <- na.omit(phenotyper)

phenotypep <- res$p#提取r数据
phenotypep <- na.omit(phenotypep)

phenotypepadj <- res$p.adj#提取r数据
phenotypepadj <- na.omit(phenotypepadj)

statistic <-res$stars
statistic <- na.omit(statistic)
#class(res)
phenotyperlong <- melt(phenotyper, id.vars = colnames(phenotyper))
colnames(phenotyperlong) <- c("GMM","SCFA","Correlation")
# View(phenotyperlong)

phenotypeplong <- melt(phenotypep, id.vars = colnames(phenotypep))
colnames(phenotypeplong) <- c("GMM","SCFA","pvalue")
# View(phenotypeplong)


phenotypepadjlong <- melt(phenotypepadj, id.vars = colnames(phenotypepadj))
colnames(phenotypepadjlong) <- c("GMM","SCFA","padjust")


phenotypestatisticlong <- melt(statistic, id.vars = colnames(statistic))
colnames(phenotypestatisticlong) <- c("GMM","SCFA","statistic")


hebingphenotyperp <- merge(phenotyperlong, phenotypeplong , by = c("GMM","SCFA"))
hebingphenotyperpj <- merge(hebingphenotyperp,phenotypepadjlong, by = c("GMM","SCFA"))
hebingphenotyperpj <- merge(hebingphenotyperpj,phenotypestatisticlong, by = c("GMM","SCFA"))
# View(hebingphenotyperp)


library(openxlsx)
write.xlsx(hebingphenotyperpj, paste0("GMM/GMM_SCFA_corr_un.all.xlsx"))
write.table(hebingphenotyperpj,file = "GMM/GMM_SCFA_corr_un.all.txt",sep = "\t",quote = F,row.names = F)

data_merge1 <- subset(hebingphenotyperpj, hebingphenotyperpj$pvalue <= 0.05)
data_merge_filtered <- subset(data_merge1, abs(data_merge1$Correlation) >= 0.5)


#保存文件（筛选的）
write.xlsx(data_merge_filtered, paste0("GMM/GMM_SCFA_corr_filter.all.xlsx"))

write.table(data_merge_filtered,file = "GMM/GMM_SCFA_corr_filter.all.txt",sep = "\t",quote = F,row.names = F)


##########################################画图

pathway_SCFA = read.table (file = "GMM/GMM_SCFA_corr_un.all.txt",  header = T, sep = "\t")
library("tidyr")
library(reshape2)
#纵轴(画图为横轴)  ~ 横轴 , value.var = "pvalues"，中间填充的值，如果缺失则为NA
#acast为将“纵轴设置为行名，而dcast行名为1,2,3,4
# linchuang_pathway_p   <- acast(linchuang_pathway, pathway~linchuang , value.var = "p.value")
# linchuang_pathway_cor <- acast(linchuang_pathway, pathway~linchuang , value.var = "correlation")

pathway_SCFA_p   <- acast(pathway_SCFA, GMM~SCFA , value.var = "pvalue")
pathway_SCFA_cor <- acast(pathway_SCFA, GMM~SCFA , value.var = "Correlation")



annotation = as.data.frame(read.table(paste0("pheatmap.GMM_GBM.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))

#替换矩阵的名字

annotation <- annotation[rownames(pathway_SCFA_p),]
# plotmat <- plotmat[rownames(annotation),]

rownames(pathway_SCFA_p) <- annotation$Level1
pathway_SCFA_p <- as.data.frame(t(pathway_SCFA_p))


pathway_SCFA_cor <- pathway_SCFA_cor[rownames(annotation),]
rownames(pathway_SCFA_cor) <- annotation$Level1
pathway_SCFA_cor <- as.data.frame(t(pathway_SCFA_cor))


maxr1 <- round(max(pathway_SCFA_cor),4)
minr1 <- round(min(pathway_SCFA_cor),4)
col_pathway_SCFA <- colorRamp2(
  c(minr1, 0, maxr1), 
  # c(  "#f1f1b8","white", "#C795FF")
  c(  "#FFFFB3","white","#8DD3C7")
)


#rm (plotmat_p)
if (!is.null(pathway_SCFA_p)){
  ssmt <- pathway_SCFA_p< 0.01
  pathway_SCFA_p[ssmt] <-'**'
  smt <- pathway_SCFA_p >0.01& pathway_SCFA_p <0.05
  pathway_SCFA_p[smt] <- '*'
  pathway_SCFA_p[!ssmt&!smt]<- ''
} else {
  pathway_SCFA_p <- F
}



cell_fun2 <-function(j, i, x, y, width, height, fill) {
  grid.text(pathway_SCFA_p[i,j], 
            x = x, y = y-height*0.3,
            gp = gpar(fontsize = 15,col="black"))}


PA <- read.table("GMM/GMM_order.txt", header = T,sep = "\t")
SCFA <- read.table("GMM/SCFA_order.txt", header = T,sep = "\t")
row_order<- SCFA$x
column_order<- PA$x

p2=Heatmap(pathway_SCFA_cor,
           name = "pathway_SCFA cor",
           #设置单元格内的星号
           cell_fun = cell_fun2,
           #设置方格长宽
           #heatmap_width和heatmap_height控制整个热图的大小（包括图例），width和height只控制热图主体的大小
           width = unit(3, "cm"),
           height = unit(2, "cm"),
           #设置行列name
           row_names_gp = gpar(fontsize = 11,col ="black"),
           column_names_gp = gpar(fontsize = 11,col ="black"),
           #聚类
           row_names_side = "left",
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           #width = 12,height = 11,
           #width = unit(2.28, "cm"),
           #设置顺序
           #row_order=c("Clostridium_bolteae","Blautia_sp_CAG_257","Clostridium_clostridioforme","Clostridium_symbiosum","Hungatella_hathewayi","Clostridium_citroniae","Anaerotruncus_colihominis","Ruminococcus_gnavus","Erysipelatoclostridium_ramosum","Sellimonas_intestinalis","Tyzzerella_nexilis","Clostridium_scindens","Dialister_sp_CAG_357","Barnesiella_intestinihominis","Oscillibacter_sp_57_20","Coprococcus_comes","Roseburia_sp_CAG_471","Faecalibacterium_prausnitzii","Bifidobacterium_adolescentis","Firmicutes_bacterium_CAG_110","Oscillibacter_sp_CAG_241","Eubacterium_ramulus","Gemmiger_formicilis"),
           column_order=column_order,
           row_order=row_order,
           #设置格子分隔
           rect_gp = gpar(col = "grey", lwd = 1),
           #设置行不聚类
           #cluster_cols = T)#,
           #设置颜色
           col = col_pathway_SCFA)#,
p2


pdf("GMM/GMM_SCFA.2.pdf",  width=5, height=5)
p2
dev.off()


png("GMM/GMM_SCFA.2.png",  width=300, height=400,units = "px")
p2
dev.off()

# write.table(row_order,file = "Pathway/SCFA_order.txt",quote = F,row.names = F)




