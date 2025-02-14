setwd("F:/0A其他工作/14、MTBI项目0812/8、个性化绘图/2、单界ROC曲线及多界箱线图/2、多界箱型图")

##处理数据
cl_list <- c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")

hebing <- NULL

for(cl in cl_list){
  file <- read.csv(paste0("mTBI_Discovery_selfcv_repeat_max_",cl,".csv"), row.names=1)
  file$Group <- cl
  
  if(is.null(hebing)){
    hebing <- file
  }else{
    hebing <- rbind(hebing,file)
  }
}

write.csv(hebing,file = "mTBI_Discovery_selfcv_repeat_max_ALL.csv",quote = F,row.names = F)


#####画图###
hebing$Group <- factor(hebing$Group, levels = cl_list)


col = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
        "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C", "#FF7F00")  

col = c("#440154", "#482878", "#3E4989", "#31688E", "#26828E", "#1F9E89", "#35B779", "#6DCD59",
        "#B4DD2B", "#FDE725", "#F9F871", "#F3E881", "#E6C88A", "#D69A87", "#CB6B60")  

col = c("#6A3D9A", "#7570B3", "#BEBADA", "#BC80BD", "#DECBE4", "#D69A87", "#ED7953", "#FB9F3A",
        "#F7D038", "#F0F921", "#D9EE42", "#ABED6E", "#6CE448", "#33A02C", "#A6CEE3")  

p1 <- ggplot(hebing,aes(x=Group,y=AUC,fill = Group ))+  
  stat_boxplot(geom = "errorbar",width = 0.15)+  
  geom_boxplot(size = 0.5,fill = "white",outlier.shape = NA)+  
  geom_jitter(aes(fill=Group),width = 0.2,shape=21,size=2.5)+  
  scale_y_continuous(limits = c(0.85, 1), breaks = seq(0.85, 1, by = 0.05))+
  theme_bw()+  
  scale_fill_manual(values = col)+  
  # scale_color_manual(values=c("black","black","black","black","black","black"))+   
  ggtitle("")+   
  theme_bw()+  
  theme(axis.title.x = element_blank(),
        panel.grid.major = element_blank(),        
        panel.grid.minor = element_blank())+
  guides(fill = FALSE)
p1


ggsave(filename =paste0("boxplot_15.3.pdf"), p1, width = 130, height = 110, units="mm")
ggsave(filename =paste0("boxplot_15.3.png"), p1, width = 130, height = 110, units="mm")

