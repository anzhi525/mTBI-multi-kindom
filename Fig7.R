#Fig7 B-----------------------------

result_df <- read.table(paste0("roc_zong.txt"),sep = "\t",header = T, comment.char="",check.names=F,quot="",stringsAsFactors = F)#读取相应模型的roc表
result_df$Group<- factor(result_df$Group ,levels = c("mTBI","Stroke","MG","AD","PD","SCZ","MDD"),ordered = TRUE)

first_row_df <- result_df %>% 
  group_by(Group) %>% 
  slice(1) 

legend_labels <- paste0(first_row_df$Group, ": ", round(first_row_df$AUC,3))

Family_col <- c("#FF6666", "#FF9966", "#FFCC66","#33A02C","#1F78B4","#7570B3","#f6a5c0" )

roc_pic2 <- ggplot(result_df, aes(x = `1-specificities`, y = sensitivities)) +
  geom_path(aes(color = Group), linewidth = 1.2, alpha = 0.7) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "darkgrey", linetype = 4) +
  scale_color_manual(values = Family_col, labels = legend_labels) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5),
    legend.justification = c(1, 0),
    legend.position = c(0.95, 0.1),
    legend.title = element_blank(),
    legend.background = element_rect(fill = NULL, linewidth = 0.2, linetype = "solid", color = "black")
  )
roc_pic2
ggsave(filename =paste0("roc_zong.pdf"), roc_pic2, width = 110, height = 100, units="mm")
ggsave(filename =paste0("roc_zong.png"), roc_pic2, width = 110, height = 100, units="mm")


#Fig7 D-----------------------------
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library("stringi")
library(circlize)

plotmat = read.table (file = "heatmap.cor.txt", row.names = 1, header = T, sep = "\t")
stars = read.table (file = "heatmap.qval.txt", row.names = 1, header = T, sep = "\t")

if (!is.null(stars)){
  ssmt <- stars< 0.01
  stars[ssmt] <-'**'
  smt <- stars >0.01& stars <0.05
  stars[smt] <- '*'
  stars[!ssmt&!smt]<- ''
} else {
  stars <- F
}

maxr <- round(max(plotmat),4)*100
minr <- abs(round(min(plotmat),4)*100)

annotation2 = as.data.frame(read.table(paste0("heatmap_anno.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation2 <- as.data.frame(annotation2[rownames(plotmat),])
colnames(annotation2) <- "Kingdom"

kindom_col <-c("#B2DF8A",  "#FCCDE5","#B3CDE3","#BEBADA")
names(kindom_col) <- unique(annotation2$Kingdom)

species_ann <- rowAnnotation(
  Kingdom= annotation2$Kingdom,
  col = list( 
    Kingdom=kindom_col
  )
)

maxr <- round(max(plotmat),4)*100
minr <- abs(round(min(plotmat),4)*100)

mycol <- c(colorRampPalette(c('#1F78B4','white'))(minr),
           colorRampPalette(c('white','white'))(3),
           colorRampPalette(c('white','#FB9A99'))(maxr))

spe <- rownames(plotmat)
PA <- c("mTBI", "Stroke", "MG","AD","PD","SCZ","MDD")
column_order <- PA
row_order <- spe

cell_fun <-function(j, i, x, y, width, height, fill) {
  grid.text(stars[i,j], 
            x = x, 
            y =  y-height*0.3,
            gp = gpar(fontsize = 15,col="black"))
}

p=Heatmap(plotmat,cluster_rows = T,
          width = unit(3.5, "cm"),
          height = unit(10, "cm"),
          row_names_gp = gpar(fontsize = 11,col ="black"),
          column_names_gp = gpar(fontsize = 11,col ="black"),
          row_names_max_width = unit(100,"mm"),
          heatmap_legend_param = list(
            at = c(-0.4,-0.2, 0, 0.2,0.4)
          ),
          cell_fun = cell_fun,
          rect_gp = gpar(col = "grey", lwd = 1),
           name = "Coef.",
          row_order=row_order,
          column_order=column_order,
           left_annotation = species_ann,
          col = mycol)
p

pdf("heatmap_zong_mtbi_qval.pdf",  width=6, height=5)
p
dev.off()

png("heatmap_zong_mtbi_qval.png",  width=500, height=450)
p
dev.off()




