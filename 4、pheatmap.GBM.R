setwd("F:/0A其他工作/14、MTBI项目0812/8、个性化绘图/4、功能与菌相关性热图")
##2.1 GMM与菌---------------------

#输入分组
metadata = read.table (file = "metadata59.new.txt", row.names = 1, header = T, sep = "\t")


GBM_abundance = read.table (file = "diff_GBM.txt", sep = "\t", row.names = 1, header = T)
GBM_abundance <- GBM_abundance[which(rowSums(GBM_abundance) > 0),]
GBM_abundance <- GBM_abundance[,rownames(metadata)]
GBM_abundance <-as.data.frame( t(GBM_abundance))



#输入菌表
metabolomic = read.table (file = "diff_species.txt", sep = "\t", row.names = 1, header = T)
metabolomic <- metabolomic[which(rowSums(metabolomic) > 0),]
metabolomic <- as.data.frame(t(metabolomic))
metabolomic <- metabolomic[rownames(GBM_abundance),]


res <-  corr.test(metabolomic,GBM_abundance,use='pairwise',method='spearman',adjust='BH',alpha=0.05)
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
colnames(phenotyperlong) <- c("species","GBM","Correlation")
# View(phenotyperlong)

phenotypeplong <- melt(phenotypep, id.vars = colnames(phenotypep))
colnames(phenotypeplong) <- c("species","GBM","pvalue")
# View(phenotypeplong)


phenotypepadjlong <- melt(phenotypepadj, id.vars = colnames(phenotypepadj))
colnames(phenotypepadjlong) <- c("species","GBM","padjust")


phenotypestatisticlong <- melt(statistic, id.vars = colnames(statistic))
colnames(phenotypestatisticlong) <- c("species","GBM","statistic")


hebingphenotyperp <- merge(phenotyperlong, phenotypeplong , by = c("species","GBM"))
hebingphenotyperpj <- merge(hebingphenotyperp,phenotypepadjlong, by = c("species","GBM"))
hebingphenotyperpj <- merge(hebingphenotyperpj,phenotypestatisticlong, by = c("species","GBM"))
# View(hebingphenotyperp)


library(openxlsx)
write.xlsx(hebingphenotyperpj, paste0("GBM/species_GBM_corr_un.all.xlsx"),
           append = FALSE,
           colNames = TRUE,
           rowNames = FALSE)
write.table(hebingphenotyperpj,file = "GBM/species_GBM_corr_un.all.txt",sep = "\t",quote = F,row.names = F)

data_merge1 <- subset(hebingphenotyperpj, hebingphenotyperpj$pvalue <= 0.05)
data_merge_filtered <- subset(data_merge1, abs(data_merge1$Correlation) >= 0.5)


#保存文件（筛选的）
write.xlsx(data_merge_filtered, paste0("GBM/species_GBM_corr_filter.all.xlsx"),
           append = FALSE,
           colNames = TRUE,
           rowNames = FALSE)

write.table(data_merge_filtered,file = "GBM/species_GBM_corr_filter.all.txt",sep = "\t",quote = F,row.names = F)

#画图


library(readxl)


write.table(phenotyper, file=paste0("GBM/heatmapGBM.cor.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)
write.table(phenotypepadj, file=paste0("GBM/heatmapGBM.padj.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)


# write.table(plotmat, file=paste0("heatmap.cor.14.26.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)
# write.table(stars, file=paste0("heatmap.padj.14.26.txt"), quote = F, sep="\t", eol = "\n", na = "NA", dec = ".",row.names = T, col.names = T)




plotmat = read.table (file = "GBM/heatmapGBM.cor.txt", row.names = 1, header = T, sep = "\t")
stars = read.table (file = "GBM/heatmapGBM.padj.txt", row.names = 1, header = T, sep = "\t")

# rows_to_delete <- c("Streptococcus_sp_263_SSPC",
#                     "Streptococcus_cristatus",
#                     "Rothia_dentocariosa",
#                     "Isoptericola_variabilis",
#                     "Klebsiella_aerogenes",
#                     "Lawsonibacter_sp_NSJ_51",
#                     "Cronobacter_sakazakii",
#                     "Pseudocitrobacter_faecalis",
#                     "Intestinibacter_bartlettii",
#                     "Saccharomyces_cerevisiae",
#                     "Parascardovia_denticolens",
#                     "Lachnoanaerobaculum_sp_ICM7",
#                     "Parvimonas_micra",
#                     "Christensenella_minuta",
#                     "Eubacterium_sp_AM28_29",
#                     "Streptococcus_mutans",
#                     "Eubacterium_sp_NSJ_61",
#                     "Candidatus_Avimonas_narfia",
#                     "Clostridium_sp_Marseille_P3244",
#                     "Peptoniphilus_sp_Marseille_P3761",
#                     "Actinomyces_oricola",
#                     "Candidatus_Pararuminococcus_gallinarum",
#                     "Campylobacter_hominis",
#                     "Roseburia_hominis",
#                     "Barnesiella_intestinihominis")

# 
# plotmat <- plotmat[!(rownames(plotmat) %in% rows_to_delete), ]
# stars <- stars[!(rownames(stars) %in% rows_to_delete), ]
# 


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
mycol <- c(colorRampPalette(c('#31688E','white'))(minr),
           colorRampPalette(c('white','white'))(3),
           colorRampPalette(c('white','#FB9A99'))(maxr))

#max.rownames0 <- max(stri_length(rownames(plotmat)))
#max.colnames0 <- max(stri_length(colnames(plotmat)))

library(pheatmap)
library(ComplexHeatmap)
#画所有的图（不筛选r的），并且显示所有p小于0.05的*


annotation = as.data.frame(read.table(paste0("pheatmap.GBM.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))

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


 annotation_col <- annotation[c("Enrichment")]
annotation_raw <- annotation2[c("Family","Enrichment")]


#填充颜色设置
# library(RColorBrewer)
# getPalette0 = colorRampPalette(brewer.pal(14, "Set3"))
# getPalette0(14)
# level1_col  <- getPalette0(14)
# level1_col  <-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462")
# names(level1_col) <- unique(annotation$Level1)


# 
# Modules_col <-c("#FED9A6","#FB9A99")
# names(Modules_col) <- unique(annotation$Modules)


Enrichment_col<- c("#1B9E77", "#D95F02")
names(Enrichment_col) = c("mTBI" , "Control" )



#设置pathway的注释
pathway_ann <- HeatmapAnnotation(
  # Level1 = annotation$Level1,
  # Modules = annotation$Modules,
  Enrichment =annotation$Enrichment,
  col = list(
    # Level1 = level1_col ,
    # Modules = Modules_col,
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
species_ann <- rowAnnotation(
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

# mycol <- c("#93CBE1","#C8EAE0","#FDE6EC","#FEBCCC","#45A8CB")
p=Heatmap(plotmat,cluster_rows = T,cluster_columns = T,
          #设置方格长宽
          #heatmap_width和heatmap_height控制整个热图的大小（包括图例），width和height只控制热图主体的大小
          width = unit(5.5, "cm"),
          height = unit(20, "cm"),
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
          name = "Species_GBM corr",
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

pdf("GBM/GBM_species2.pdf",  width=10, height=15)
p
dev.off()


png("GBM/GBM_species2.png",  width=700, height=1000,units = "px")
p
dev.off()

write.table(names(column_order(p)),file = "GBM/GBM_order.txt",quote = F,row.names = F)
 write.table(names(row_order(p)),file = "GBM/species_order.txt",quote = F,row.names = F)

