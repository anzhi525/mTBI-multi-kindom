#Fig3 B-----------------------------
library(ggalluvial)
library(RColorBrewer)

data <- read.table(file ="pathway.txt",header = T,sep="\t", comment.char="",quot="",check.names=F)
sankey_colors<-  colors <- c("#FFFF99", "#CCEBC5","#B2DF8A", "#A6CEE3","#FDDAEC", "#BEBADA")
p.sankey.ED <- ggplot(data = data,                      
                      aes(axis1 = data$Level1,                           
                          axis3 = data$Level3)) +  
  scale_x_discrete(limits = c("", "", "")) +  
  geom_alluvium(aes(fill = data$Level1, alpha =1),curve_type = "arctangent") +  
  scale_color_manual(values = sankey_colors)+  scale_fill_manual(values = sankey_colors)+  
  geom_stratum(alpha = 0,color = adjustcolor( "white", alpha.f = 1),size=1.2, fill = 'white') +  
   geom_text(stat = "stratum",cex=6, aes(label = after_stat(stratum)),hjust = 1) +  

  theme_void()+  
  theme(legend.position="none",
        axis.text = element_text(size = 16),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid=element_blank())
p.sankey.ED

breaklist <- seq(-1,1,by=0.001)
purple_green <- rev(brewer.pal(n=11,name="PiYG"))
col_purple_green <- colorRampPalette(purple_green)(length(breaklist))

color <-c("#FFFF99", "#CCEBC5","#B2DF8A", "#A6CEE3","#FDDAEC", "#BEBADA")
data$Level3 <- factor(data$Level3, levels = rev(unique(data$Level3)))

p_bub <- ggplot(data,
                aes(x=coef,y=Level3,color=-log10(qval),size=5))+#,size=genenum_cat
  geom_point(alpha=0.7, shape=19)+
  scale_x_continuous(breaks = seq(-1.5,1.5, by = 0.5))+
  labs(x="Coef",y="")+
  scale_color_gradientn(colours = col_purple_green,name="-log10(FDR)")+
  scale_size(guide = "none")+
  theme_bw() +
  theme(
    axis.title.y = element_text(size = 15),
    axis.title = element_text(size = 13),
    axis.text.x = element_text(angle = 90, vjust = 0.5,size = 14),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14)
  )

p_bub

p_cat = ggdraw()+ draw_plot(p.sankey.ED)+
        draw_plot(p_bub,scale=0.51,x=0.36,y=-0.52,width=0.85,height=1.98)

ggsave(file=paste0("KEGGpathway.pdf") ,p_cat,width = 260, height = 200, units="mm",limitsize = FALSE)
ggsave(file=paste0("KEGGpathway.png"),p_cat,width = 260, height = 200, units="mm",limitsize = FALSE)

#Fig3 C-----------------------------
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(readxl)

df <- read_xlsx(path = "GBM.xlsx", col_names = T, sheet = 1) %>%
  dplyr::mutate(start = 0, 
                end = N)

pdf(file = "GBM.pdf",
    height =2,
    width = 2)
circos.par("start.degree" = 90,
           "track.margin" = c(0.02, 0.02),
           "cell.padding" = c(0.02, 0.02, 0, 0))

df1 <- df %>% 
  dplyr::select(feature, start, end, Description) %>%
  dplyr::rename(ID = feature) %>%
  dplyr::mutate(Description = str_remove(Description, pattern = "\\s\\(.*"))

feature_color <- c("#D51F26", "#272E6A", "#208A42", "#89288F", "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC",
                 "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", "#3BBCA8", "#6E4B9E", "#0C727C",
                 "#D8A767", "#7DD06F", "#844081", "#688EC1", "#C17E73", "#6CD3A7", "#597873", "#7B6FD0",
                 "#CF4A31", "#D0CD47", "#CBC594", "#D19EC4", "#5A7E36", "#D4477D", "#403552", "#76D73C",
                 "#96CED5", "#CE54D1", "#C48736", "#FFB300", "#803E75", "#FF6800", "#A6BDD7", "#C10020", "#CEA262",
                 "#817066", "#007D34", "#F6768E", "#00538A", "#FF7A5C", "#53377A", "#FF8E00", "#B32851", "#F4C800",
                 "#7F180D", "#93AA00", "#593315", "#F13A13", "#232C16", "#faa818", "#41a30d", "#fbdf72", "#367d7d",
                 "#d33502", "#6ebcbc", "#37526d", "#916848", "#f5b390", "#342739", "#bed678", "#a6d9ee", "#0d74b6",
                 "#60824f", "#725ca5", "#e0598b", "#371377", "#7700FF", "#9E0142", "#FF0080", "#DC494C", "#F88D51",
                 "#FAD510", "#FFFF5F", "#88CFA4", "#238B45", "#02401B", "#0AD7D3", "#046C9A", "#A2A475", "grey35",
                 "#D52126", "#88CCEE", "#FEE52C", "#117733", "#CC61B0", "#99C945", "#2F8AC4", "#332288", "#E68316",
                 "#661101", "#F97B72", "#DDCC77", "#11A579", "#E73F74", "#A6CDE2", "#1E78B4", "#74C476", "#34A047",
                 "#F59899", "#E11E26", "#FCBF6E", "#F47E1F", "#CAB2D6", "#6A3E98", "#FAF39B", "#B15928", "#1a1334",
                 "#01545a", "#017351", "#03c383", "#aad962", "#fbbf45", "#ef6a32", "#ed0345", "#a12a5e", "#710162",
                 "#3B9AB2", "#2a7185", "#a64027", "#9cdff0", "#022336", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00",
                 "#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6")

circos.initializeWithIdeogram(df1, plotType = NULL)
circos.track(
  ylim = c(0, 1), 
  track.height = 0.08,  
  bg.border = NA,  
  bg.col = feature_color, 
  panel.fun = function(x, y) {
    ylim = get.cell.meta.data("ycenter")  
    xlim = get.cell.meta.data("xcenter")
    sector.name = get.cell.meta.data("sector.index")  
    track.index = get.current.track.index()
    ideogram.height=10
  } )

df2 <- df %>% 
  dplyr::select(feature, start, end, pval) %>%
  dplyr::rename(ID = feature) %>%
  dplyr::mutate(`-log10pvalue` = -log10(pval)) %>%
  dplyr::select(1,2,3,5)

summary(df2$`-log10pvalue`)

col_fun1 = colorRamp2(breaks = c(1, 2, 4, 6, 8, 10), colors =c("#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494"))

circos.genomicTrackPlotRegion(
  df2,
  track.height = 0.09, 
  bg.border = NA, 
  stack = TRUE,  
  
  panel.fun = function(region, value, ...) {
    circos.genomicRect(
      region, 
      value, 
      col = col_fun1(value[[1]]), 
      border = NA, ...
      ) 
    ylim = get.cell.meta.data("ycenter")  
    xlim = get.cell.meta.data("xcenter")
    sector.name = get.cell.meta.data("sector.index") 
    circos.text(xlim, ylim + 1, sector.name, cex = 0.3, niceFacing = FALSE)  # 添加 GO Term 
  } )

df3 <- df %>%
  dplyr::select(feature, start, end, coef) %>%
  dplyr::rename(ID = feature)

col_fun2 = colorRamp2(breaks = c(-2, 0, 1), colors =c("#33A02C","#ffffff","#482878"))

summary(df3$coef)

df3_up <- df3 %>% dplyr::filter(coef > 0)
df3_down <- df3 %>% dplyr::filter(coef < 0)

circos.genomicTrack(
  df3_up, 
  ylim = c(0.05, 1), 
  track.height = 0.45, 
  bg.col = "#f0f0f0", 
  bg.border = NA,  
  track.margin = c(0, 0), 
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data("sector.index")  
    circos.genomicRect(region, value, 
                       col = col_fun2(value[[1]]), 
                       border = NA, 
                       ytop.column = 1, ybottom = 0,
                       ...) 
  } )

circos.genomicTrack(
  df3_down, 
  ylim = c(-1.5, 0), 
  track.height = 0.25, 
  bg.col = NA, 
  bg.border = NA, 
  track.margin = c(0, 0), 
  panel.fun = function(region, value, ...) {
    sector.name = get.cell.meta.data("sector.index")  
    circos.genomicRect(region, value, 
                       col = col_fun2(value[[1]]), 
                       border = NA,
                       ytop = 0,
                       ytop.column = 0, 
                       ybottom = -1.5,
                       ybottom.column = -1.5,
                       ...) 
  } )
circos.clear()

pvalue_legend <- Legend(
  title = "",
  labels = rev(c(1, 2, 4, 6, 8, 10)),
  type = "points", pch = NA,
  background = rev(c("#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")),
  labels_gp = gpar(fontsize = 5), grid_height = unit(0.3, "cm"), grid_width = unit(0.3, "cm"))

legend_list <- lgd_list_vertical <- packLegend(pvalue_legend)
pushViewport(viewport(x = 0.9, y = 0.13))
grid.draw(lgd_list_vertical)
upViewport()

dev.off()