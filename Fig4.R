#Fig5 B-----------------------------
setwd("09_feature_ggpub")
library("RColorBrewer")
library("ggplot2")
library(grid)
library(cowplot)

breaklist <- seq(-1,1,by=0.001)
purple_green <- rev(brewer.pal(n=11,name="PiYG"))
col_purple_green <- colorRampPalette(purple_green)(length(breaklist))

color <- c("#8dd3c7","#bebada","#80b1d3","#fccde5",
           "#feb24c","#33A02C","#B2DF8A","#FFFF99")

data <- read.table(file = paste0("bacteria.txt"),header = T,sep="\t", comment.char="",quot="",check.names=F)
   
data = as.data.frame(data[(order(data$coef)), ])
data$feature <- factor(data$feature,levels = data$feature)

p_bub <- ggplot(data,
                aes(x=feature,y=coef,color=-log10(qval),size=5))+#,size=
  geom_point(alpha=0.7, shape=19)+

  labs(x="",y="Coef")+
  scale_color_gradientn(colours = col_purple_green,name="-log10(FDR)")+
  scale_size(guide = "none")+
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 18),
    axis.title = element_text(size = 13),
    axis.text.x = element_blank(),  
    axis.text.y = element_text(size = 15),
    legend.title = element_text(size = 18), 
    legend.text = element_text(size = 15) 
  )
p_bub


p_cat <- ggplot(data,aes(y=0.1,x=feature))+
  geom_tile(aes(fill=Phylum),width=1, height=0.2)+ 
  labs(fill="Phylum")+ 
  scale_fill_manual(values=color) +
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 15,angle = 90,vjust = 0.5,hjust = 1,
                                   colour = "black"), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15)
        )
p_cat

p <- plot_grid(
	  p_bub,
	  NULL,
	  p_cat,
	  align = "v",
	  ncol = 1,
	  rel_heights = c(5, -0.35, 4)
	)

ggsave(file=paste0("bacteria.feature.pdf"),p,width = 10, height = 6.8,limitsize = FALSE)
ggsave(file=paste0("bacteria.feature.png"),p,width = 10, height = 6.8,limitsize = FALSE)


#Fig5 C-----------------------------
result_df <- read.csv(file = paste0("roc.csv"),header = T)

result_df$group <- factor(result_df$group ,
                       levels = c("A","B","V","F"),ordered = TRUE)
first_row_df <- result_df %>% group_by(group) %>% slice(1) 
legend_labels <- paste0(first_row_df$group, ": ", round(first_row_df$AUC,3))

getPalette_family = colorRampPalette(brewer.pal(4, "Set1"))
Family_col <- getPalette_family(4)

roc_pic2 <- ggplot(result_df, aes(x = 1-specificities, y = sensitivities)) +
  geom_path(aes(color = group), linewidth = 1.2, alpha = 0.7) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color = "darkgrey", linetype = 4) +
  scale_color_manual(values = Family_col, labels = legend_labels) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(hjust = 0.5),
    legend.justification = c(1, 0),
    legend.position = c(0.95, 0.25),
    legend.title = element_blank(),
    legend.background = element_rect(fill = NULL, linewidth = 0.2, linetype = "solid", color = "black")
  )
roc_pic2
ggsave(path = paste0(outdir), filename =paste0("roc_plot.pdf"), roc_pic2, width = 110, height = 110, units="mm")
ggsave(path = paste0(outdir), filename =paste0("roc_plot.png"), roc_pic2, width = 110, height = 110, units="mm")


#Fig5 D-----------------------------
library(ggplot2)
library(ggsci)
pathway<-read.table('V_top20.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
pathway = as.data.frame(pathway[(order(-pathway$importance)), ])
rownames(pathway) <- pathway$species
pathway$species <- factor(pathway$species,levels = rev(pathway$species))

my_cols<-c("#4477b5","#acd1e1","#f6f6ca","#f9cc81","#da342e")

p <- ggplot(data = pathway, 
            aes(x = importance, y = species,fill = importance))+ 
  geom_bar(stat = "identity",width = 0.7,colour="black") +
  scale_fill_gradientn(colors=my_cols)+
  labs(x = "Importance",y = "",title = "")+ 
  geom_vline(xintercept=0,lty=2,col="grey60",lwd=1.5) +
  theme(plot.title = element_text(size=15,face="bold" ,hjust = 0.5),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "right",
        legend.text = element_text(size = 10),
        axis.text.x = element_text(size = 9,
                                   family = "sans"),
        axis.text.y = element_text(family = "sans",size = 12),
        legend.title = element_blank(),
        text = element_text(family = "sans",size = 10),
        strip.text.x = element_text(size = 0.5), 
        strip.text= element_text(family = "sans",size = 7,face="italic"),
        panel.background = element_rect(fill = 'grey92'),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(
          colour="grey50",
          size=.2,
          linetype=2,
          lineend=1),
        legend.key = element_rect(fill = 'white'))
p
ggsave(paste0("V_top20.pdf"), p, width = 200, height = 180, units="mm")
ggsave(paste0("V_top20.png"), p, width = 200, height = 180, units="mm")

