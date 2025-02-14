setwd("F:/0A其他工作/14、MTBI项目0812/8、个性化绘图/1、火山图")

###加载包###
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(readxl)
###导入数据###

data <- read_excel("Maaslin2.xlsx", sheet = "筛选0.001")

data

cl_list <- c("Archaea","Bacteria","Fungi","Virus")
###数据处理——根据FDR进行差异基因筛选###
for(cl in cl_list){
  ###数据处理——根据FDR进行差异基因筛选###
  table <- read_excel("Maaslin2.xlsx", sheet = "筛选0.001")
  
  data <- subset(table,kindom==cl)
  
cut_off_FDR =0.05 #设置FDR的阈值
cut_off_log2FC =0 #设置log2FC的阈值
data$Sig <- ifelse(data$qval < cut_off_FDR & data$coef > cut_off_log2FC, "Up",
                   ifelse(data$qval < cut_off_FDR & data$coef < cut_off_log2FC, "Down", "no"))
data = data.frame(data)
table(data$Sig) #查看数据统计情况

###绘图——基础火山图###"#31688E""#D53E4F"
p1 <- ggplot(data, aes(x =coef , y=-log10(qval), colour=Sig)) + #x、y轴取值限制，颜色根据"Sig"
  geom_point(alpha=0.65, size=2) +  #点的透明度、大小
  scale_color_manual(values=c("#1F78B4", "#d2dae2","#D53E4F")) + xlim(c(-3, 3)) +  #调整点的颜色和x轴的取值范围
  geom_vline(xintercept=c(-cut_off_log2FC,cut_off_log2FC),lty=4,col="black",lwd=0.8) + #添加x轴辅助线,lty函数调整线的类型："twodash"、"longdash"、"dotdash"、"dotted"、"dashed"、"solid"、"blank"
  geom_hline(yintercept = -log10(cut_off_FDR), lty=4,col="black",lwd=0.8) +  #添加y轴辅助线
  labs(x="Coef.(by MaAsLin2)", y="-log10FDR") +  #x、y轴标签
  ggtitle(cl) + #标题
  theme_bw() + # 主题，help(theme)查找其他个性化设置
  theme(plot.title = element_text(hjust = 0),
        legend.position="right", 
        legend.title = element_blank()
  ) 
p1 #出图



###添加基因名标记###
p2 <- p1 + geom_text_repel(
  data = subset(data, data$qval < cut_off_FDR & abs(data$coef) >= cut_off_log2FC),# 可以设置跟上面不同的阈值，用数值替换即可
  aes(label = feature), size = 3,
  box.padding = unit(0.5, "lines"),
  point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )
p2 #出图


###添加差异最显著的上调和下调的前几个基因###
#下载与加载包#

library(dplyr) # 用于数据处理
library(gt) # 制作表格
#数据预处理——上下调合并绘制同一类型标签框#
# 从原始数据中筛选出差异显著的上调基因中FDR值最小的前10个，并按log2FC的绝对值降序排列
# 从原始数据中筛选出差异显著的下调基因中FDR值最小的前10个，并按log2FC的绝对值降序排列
# 将上调和下调基因的结果合并成一个包含前20个基因的数据框
top_20 <- bind_rows(
  data %>%
    filter(Sig == 'Up') %>%
    arrange(qval, desc(abs(coef))) %>%
    head(5),
  data %>%
    filter(Sig == 'Down') %>%
    arrange(qval, desc(abs(coef))) %>%
    head(5)
)
top_20 %>% gt()  #将数据制成表
#绘图——添加基因标签框图#
p4 <- p1 +
  geom_text_repel(data = top_20,
                   aes(coef, -log10(qval), label = feature),
                   size = 3 ,
                   box.padding = unit(0.5, "lines"),
                   point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE)
p4

ggsave( filename =paste0(cl,".pdf"), p4, width = 150, height = 150, units="mm",limitsize = FALSE)
ggsave(filename =paste0(cl,".png"), p4, width = 150, height = 150, units="mm",limitsize = FALSE)
}