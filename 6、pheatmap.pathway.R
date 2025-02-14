setwd("F:/0A其他工作/14、MTBI项目0812/8、个性化绘图/4、功能与菌相关性热图")


##2.1 PATHWAY与差异菌---------------------

#输入分组
metadata = read.table (file = "metadata59.new.txt", row.names = 1, header = T, sep = "\t")

species = read.table(paste0("diff_species.txt"), header=T, row.names=1, sep="\t", comment.char="")
species <- species[,rownames(metadata)]
species_filtered1 <- t(species)

Pathway = t(read.table (file = "diff_pathway.txt", sep = "\t", row.names = 1, header = T))
Pathway <- Pathway[which(rowSums(Pathway) > 0),]
data_filtered1 <- Pathway[rownames(species_filtered1),]

library(psych)
res <-  corr.test(data_filtered1,species_filtered1,use='pairwise',method='spearman',adjust='BH',alpha=0.05)
names(res)#查看列表每个项目

res$stars

#
# 保存需要的文件
library(reshape2)
phenotyper <- as.data.frame(res$r)#提取r数据
phenotyper <- na.omit(phenotyper)
phenotyper <- as.matrix(phenotyper)

phenotypep <- as.data.frame(res$p)#提取r数据
phenotypep <- na.omit(phenotypep)
phenotypep <- as.matrix(phenotypep)

phenotypepadj <- as.data.frame(res$p.adj)#提取r数据
phenotypepadj <- na.omit(phenotypepadj)
phenotypepadj <- as.matrix(phenotypepadj)

statistic <-as.data.frame(res$stars)
statistic <- na.omit(statistic)
statistic <- as.matrix(statistic)
#class(res)
phenotyperlong <- melt(phenotyper, id.vars = colnames(phenotyper))
colnames(phenotyperlong) <- c("Pathway","species","Correlation")
# View(phenotyperlong)

phenotypeplong <- melt(phenotypep, id.vars = colnames(phenotypep))
colnames(phenotypeplong) <- c("Pathway","species","pvalue")
# View(phenotypeplong)


phenotypepadjlong <- melt(phenotypepadj, id.vars = colnames(phenotypepadj))
colnames(phenotypepadjlong) <- c("Pathway","species","padjust")


phenotypestatisticlong <- melt(statistic, id.vars = colnames(statistic))
colnames(phenotypestatisticlong) <- c("Pathway","species","statistic")


hebingphenotyperp <- merge(phenotyperlong, phenotypeplong , by = c("Pathway","species"))
hebingphenotyperpj <- merge(hebingphenotyperp,phenotypepadjlong, by = c("Pathway","species"))
hebingphenotyperpj <- merge(hebingphenotyperpj,phenotypestatisticlong, by = c("Pathway","species"))
# View(hebingphenotyperp)


library(openxlsx)
write.xlsx(hebingphenotyperpj, paste0("Pathway/species_Pathway_corr_un.all.xlsx"),
           append = FALSE,
           colNames = TRUE,
           rowNames = FALSE)
write.table(hebingphenotyperpj,file = "Pathway/species_Pathway_corr_un.all.txt",sep = "\t",quote = F,row.names = F)

data_merge1 <- subset(hebingphenotyperpj, hebingphenotyperpj$pvalue <= 0.05)
data_merge_filtered <- subset(data_merge1, abs(data_merge1$Correlation) >= 0.5)


#保存文件（筛选的）
write.xlsx(data_merge_filtered, paste0("Pathway/species_Pathway_corr_filter.all.xlsx"),
           append = FALSE,
           colNames = TRUE,
           rowNames = FALSE)

write.table(data_merge_filtered,file = "Pathway/species_Pathway_corr_filter.all.txt",sep = "\t",quote = F,row.names = F)

#画图


library(readxl)


write.table(phenotyper, file=paste0("Pathway/heatmapPathway.cor.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)
write.table(phenotypepadj, file=paste0("Pathway/heatmapPathway.padj.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)

plotmat <- as.matrix(read.table (file = "Pathway/heatmapPathway.cor.txt", row.names = 1, header = T, sep = "\t", comment.char="",check.names=F,quot="", stringsAsFactors = F))
plotmat[is.na(plotmat)] <- 0.0000
stars <- as.matrix(read.table (file = "Pathway/heatmapPathway.padj.txt", row.names = 1, header = T, sep = "\t", comment.char="",check.names=F,quot="", stringsAsFactors = F))
stars[is.na(stars)] <- 1

library(reshape2)
plotmat1<-melt(
  plotmat,                       #待转换的数据集名称
  #id.vars=c("Fc gamma R_mediated phagocytosis","biosynthesis_globo series"),  #要保留的主字段
  #variable.name="Year", 
  varnames = c("Var1", "Var2"),#转换后的分类字段名称（维度）
  value.name="cor"#转换后的度量值名称
)

stars1<-melt(
  stars,                       #待转换的数据集名称
  #id.vars=c("Fc gamma R_mediated phagocytosis","biosynthesis_globo series"),  #要保留的主字段
  #variable.name="Year",         #转换后的分类字段名称（维度）
  value.name="padj"#转换后的度量值名称
)




data_totol<-merge(plotmat1,stars1,by=c("Var1","Var2"),all=T)
write.table(data_totol, file=paste0("Pathway/pathway.species.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)




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
# mycol <- c("#FEBCCC""#93CBE1","#C8EAE0",)
mycol <- c(colorRampPalette(c('#31688E','white'))(minr),
           colorRampPalette(c('white','white'))(3),
           colorRampPalette(c('white','#FB9A99'))(maxr))


#max.rownames0 <- max(stri_length(rownames(plotmat)))
#max.colnames0 <- max(stri_length(colnames(plotmat)))

library(pheatmap)
library(ComplexHeatmap)
#画所有的图（不筛选r的），并且显示所有p小于0.05的*


annotation = as.data.frame(read.table(paste0("pheatmap.pathway.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))

#替换矩阵的名字

annotation <- annotation[rownames(plotmat),]
# plotmat <- plotmat[rownames(annotation),]

rownames(plotmat) <- annotation$Level3

stars <- stars[rownames(annotation),]
rownames(stars) <- annotation$Level3


annotation2 = as.data.frame(read.table(paste0("pheatmap.species.ann.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation2 <- annotation2[colnames(plotmat),]
# 
library(RColorBrewer)



#去r的最大值与最小值，并且保留4位数然后*100，且将负值取绝对值
maxr <- round(max(plotmat),4)
minr <- round(min(plotmat),4)
#?colorRamp2()
library(circlize)
##这个是原本默认配色的方案
# col_main <- colorRamp2(
#   c(minr, 0, maxr), 
#   c(  "#3B4992FF","white","#EE0000FF")
# )

# 
# annotation_col <- annotation[c("Level1","Enrichment")]
# annotation_raw <- annotation2[c("Family","Enrichment")]
# 

#填充颜色设置
library(RColorBrewer)
getPalette0 = colorRampPalette(brewer.pal(6, "Set3"))
getPalette0(6)
level1_col  <- getPalette0(6)
# level1_col <- c("#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC")
names(level1_col) <- unique(annotation$Level1)


getPalette1 = colorRampPalette(brewer.pal(28, "Set1"))
level2_col <-getPalette1(28)
names(level2_col) <- unique(annotation$Level2)


Enrichment_col<- c("#1B9E77", "#D95F02")
names(Enrichment_col) = c("mTBI" , "Control" )



#设置pathway的注释
pathway_ann <- rowAnnotation(
  Level1 = annotation$Level1,
  Level2 = annotation$Level2,
  Enrichment =annotation$Enrichment,
  col = list(
    Level1 = level1_col ,
    Level2 = level2_col,
    Enrichment = Enrichment_col #分类
  )
)

#设置kindom颜色
getPalette_kindom = colorRampPalette(brewer.pal(4, "Set1"))
kindom_col <- getPalette_kindom(4)#Kingdom
# kindom_col <-c("#35B779","#EA0F6B")
names(kindom_col) <- unique(annotation2$Kingdom)


#设置phylum颜色
getPalette_phylum = colorRampPalette(brewer.pal(9, "Set3"))
Phylum_col <- getPalette_phylum(9)
# Phylum_col <- c("#3b374c", "#44598e", "#64a0c0", "#7ec4b7", "#deebcd")
names(Phylum_col) <- unique(annotation2$Phylum)


#设置family颜色
getPalette_family = colorRampPalette(brewer.pal(19, "Set2"))
Family_col <- getPalette_family(19)#Family_col <- getPalette_family(10)
names(Family_col) <- unique(annotation2$Family)



#设置菌的注释
species_ann <- HeatmapAnnotation(
  Kingdom= annotation2$Kingdom,
  Phylum = annotation2$Phylum,
  Family = annotation2$Family,
  Enrichment =annotation2$Enrichment,
  col = list( 
    Kingdom=kindom_col,
    Phylum = Phylum_col , 
    Family = Family_col,
    Enrichment = Enrichment_col #分类
  )
)

# column_order<-c("Fc gamma R_mediated phagocytosis","Homologous recombination","Pancreatic cancer","Melanoma","Glycerophospholipid metabolism","Natural killer cell mediated cytotoxicity","Chlorocyclohexane and chlorobenzene degradation","Longevity regulating pathway_multiple species","Notch signaling pathway","Lysosome","Nitrotoluene degradation","Amino sugar and nucleotide sugar metabolism","Phototransduction","Glycosaminoglycan degradation","Colorectal cancer","Glycosphingolipid biosynthesis_globo series","Glycosphingolipid biosynthesis_ganglio series")
# row_order<-c("Clostridium_bolteae","Blautia_sp_CAG_257","Clostridium_clostridioforme","Clostridium_symbiosum","Hungatella_hathewayi","Clostridium_citroniae","Ruminococcus_gnavus","Erysipelatoclostridium_ramosum","Sellimonas_intestinalis","Tyzzerella_nexilis","Clostridium_scindens","Barnesiella_intestinihominis","Firmicutes_bacterium_CAG_110","Oscillibacter_sp_CAG_241","Gemmiger_formicilis","Anaerotruncus_colihominis","Dialister_sp_CAG_357","Oscillibacter_sp_57_20","Coprococcus_comes","Roseburia_sp_CAG_471","Faecalibacterium_prausnitzii","Bifidobacterium_adolescentis","Eubacterium_ramulus")

pa <- rownames(plotmat)#菌
spe <- colnames(plotmat)#通路
column_order <- spe
row_order <- pa

cell_fun <-function(j, i, x, y, width, height, fill) {
  grid.text(stars[i,j], 
            x = x, 
            y =  y-height*0.3,#把y坐标向下移一点
            gp = gpar(fontsize = 15,col="black"))
}


#方案1
# mycol <- c("#FEBCCC""#93CBE1","#C8EAE0",)
# 方案2
# mycol <- c("#547298","#8E9DBA","#DADEE7","#F1DDD6","#DA909B","#D24846")
#?Heatmap()
p=Heatmap(as.matrix(plotmat),cluster_columns = T,cluster_rows = T,
          #设置方格长宽
          #heatmap_width和heatmap_height控制整个热图的大小（包括图例），width和height只控制热图主体的大小
          width = unit(20, "cm"),
          height = unit(36, "cm"),
          #width = 12,height = 11,
          #设置行列name
          row_names_gp = gpar(fontsize = 11,col ="black"),
          column_names_gp = gpar(fontsize = 11,col ="black"),
          row_names_max_width = unit(110,"mm"),##行名与图例距离
          heatmap_legend_param = list(
            at = c(-0.4,-0.2, 0, 0.2,0.4)
          ),
          #设置行不聚类
          #cluster_cols = T)#,
          #设置单元格内的星号
          cell_fun = cell_fun,
          #设置格子分隔
          rect_gp = gpar(col = "grey", lwd = 1),
          #图例name
          name = "correlation",
          #聚类树的高度 和 宽度
          #column_dend_height/row_dend_widht
          #设置注释
          top_annotation = species_ann ,
          #bottom_annotation = column_ha, #对应的注释
          left_annotation = pathway_ann,
          row_order=row_order,
          column_order=column_order,
          #设置颜色
          col = mycol)#,
p

row_order(p)
column_order(p)


# ggsave(file ="pathway_jun.pdf", p, width = 400, height = 400, units="mm")

pdf("Pathway/pathway_species.pdf",  width=20, height=20)
p
dev.off()


png("pathway/pathway_species.png",  width=1600, height=1500,units = "px")
p
dev.off()


write.table(names(column_order(p)),file = "Pathway/species_order.txt",quote = F,row.names = F)
write.table(names(row_order(p)),file = "Pathway/pathway_order.txt",quote = F,row.names = F)

