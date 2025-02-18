#Fig6 B-----------------------------
hebing <- read.table(file = paste0("bacteria.txt"),header = T)
hebing$Group <- factor(hebing$Group, levels = cl_list)

col = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
        "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FF7F00")  

p1 <- ggplot(hebing,aes(x=Group,y=AUC,fill = Group ))+  
  stat_boxplot(geom = "errorbar",width = 0.15)+  
  geom_boxplot(size = 0.5,fill = "white",outlier.shape = NA)+  
  geom_jitter(aes(fill=Group),width = 0.2,shape=21,size=2.5)+  
  scale_y_continuous(limits = c(0.95, 1), breaks = seq(0.95, 1, by = 0.05))+
  theme_bw()+  
  scale_fill_manual(values = col)+  
  ggtitle("")+   
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(,size=13,angle = 45, hjust = 1),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15),
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank())+
  guides(fill = FALSE)
p1

ggsave(filename =paste0("boxplot_15.3.pdf"), p1, width = 130, height = 110, units="mm")
ggsave(filename =paste0("boxplot_15.3.png"), p1, width = 130, height = 110, units="mm")


#Fig6 C-----------------------------
df_wide_new1 <- as.matrix(read.table(paste0("feature_all.txt"),header = T, row.names = 1,sep = '\t',quot="",check.names=F))

annotation_row1 = as.data.frame(read.table(paste0("anno.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation_row1$Group<- factor(annotation_row1$Group,
                               levels = c("one","two","three","four"),ordered = TRUE)


annotation_row2 = as.data.frame(read.table(paste0("feature_anno.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation_row2$Kindom<- factor(annotation_row2$Kindom,
                                levels = c("Bacteria","Virus","Archaea","Fungi"),ordered = TRUE)

ann_colors = list(
  Group = c( one= "#3E4989", two = "#386CB0", three = "#A6CEE3", four = "#CCEBC5"),
  Kindom = c( Bacteria = "#7570B3" , Virus = "#BC80BD" , Archaea  = "#DECBE4",Fungi  =  "#FDDAEC"),
  Change= c(up = "#DE77AE",down = "#7FBC41"),
  importance =   c(colorRampPalette(c('#3B4992FF','white'))(60),
                           colorRampPalette(c('white','white'))(29),
                           colorRampPalette(c('white','#EE0000FF'))(60)))
mycol <- c(colorRampPalette(c("#31688E","#FED9A6"))(100))

df_wide_new2 <- df_wide_new1[rownames(annotation_row2),]
df_wide_new2[is.na(df_wide_new2)] <- 0
df_wide_new2<- df_wide_new2[,c("A","B","F","V",
                               "AB","AF","AV","BV","FV","BF",
                               "ABV","ABF","AFV","BFV",
                               "ABFV")]
library(ggplot2)
library(pheatmap)
p3=pheatmap(t(df_wide_new2),
            cellwidth = 8,cellheight = 8,
            cluster_cols = FALSE,
            cluster_rows = FALSE,
            color = mycol,
            fontsize_number = 12,
            number_color = "white",
            fontsize_row=11,
            fontsize_col=8,
            annotation_row = annotation_row1,
            annotation_col = annotation_row2,
            annotation_colors = ann_colors
)
p3

ggsave(filename =paste0("heatmap.feature.pdf"), p3, width = 450, height = 130, units="mm",limitsize = FALSE)
ggsave(filename =paste0("heatmap.feature.png"), p3, width = 450, height = 130, units="mm",limitsize = FALSE)

