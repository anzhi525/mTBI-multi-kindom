setwd("F:/0A其他工作/14、MTBI项目0812/3、多界预测模型/mTBI_code_Results_最终")

#先把15个表的20个合为一体
out <- matrix(NA, 15, 3)          #一个9*111的空矩阵
all10data<- NA

cl_list <- c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")

# cl <- "A"
for(i in 1:length(cl_list)){
  library(readxl)
  auc <- read.csv(paste0(cl_list[i],"/mTBI_Discovery_selfcv_repeat_max_",cl_list[i],".csv"), header=T,row.names=1,  comment.char="",check.names=F,quot="",fileEncoding="utf8")
  auc_mean <- mean(auc$`"AUC"`)
  # pvalue <- as.data.frame(read.csv2(paste0("D:/rwork/刘睿娜/抑郁症/结果交付1/预测/",sl,"/",cl,"/MDD_Discovery_model_A3_sig_",cl,".csv"), header=T,row.names=1,sep=",",  comment.char="",check.names=F,quot="",fileEncoding="utf8"))
  # pvalue$`"p value"`[1]
  res <- c(auc$`"AUC"`,auc_mean)
  data.frame(AUC=res,
                Group = rep(cl_list[i], 21))
  all10data <- rbind(all10data,res)
    
}
all10data <- as.data.frame(all10data)
all10data = all10data[-1,] 


rownames(all10data) <- cl_list
colnames(all10data) <- c(1:20,"mean_auc")
save_data <- t(all10data)

write.table(save_data, file=paste0("3、heatmap/mean.auc_all.txt"), append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)







all10data
cl_list1 <- c("A","B","F","V",
              "AB","AF","AV","BV","FV","BF",
              "ABV","ABF","AFV","BFV",
              "ABFV")
PLOT_DATA <- all10data[cl_list1,]

cl_list2 <- c("A","B","F","V",
              "AB","AF","AV","BV","FV","BF",
              "ABV","ABF","AFV","BFV",
              "ABFV")
rownames(PLOT_DATA) <- cl_list2

all10data_2 <- round(all10data, 3)

all10data_2 <- all10data_2[cl_list1,]

all10data_2_text = t(read.table(paste0("3、heatmap/mean.auc_all_text.txt"), header=T, row.names=1, sep="\t", comment.char="",quot="",check.names=F))
all10data_2_text2<- all10data_2_text[cl_list1,]

#画图============================================
library(pheatmap)
library(ggplot2)

#画图需要用到r数据和pvalue数据（res$r和res$p）
#提取r值-------------------------------
# metaboliter <- res$r#提取r数据
# 
# #提取p值----------------------------------
# metabolitep <- res$p#提取r数据
# 
# if (!is.null(metabolitep)){
#   ssmt <- metabolitep< 0.01
#   metabolitep[ssmt] <-'**'
#   smt <- metabolitep >0.01& metabolitep <0.05
#   metabolitep[smt] <- '*'
#   metabolitep[!ssmt&!smt]<- ''
# } else {
#   metabolitep <- F
# }


#去r的最大值与最小值，并且保留4位数然后*100，且将负值取绝对值
# maxr <- round(max(metaboliter),4)*100
# minr <- abs(round(min(metaboliter),4)*100)

#给刻度做一个分割
# bk <- seq(-1,1,by =0.01)
# length(bk)
# bk
library("stringi")
#设置颜色
#EE0000FF", high = "#3B4992FF" 

#设置颜色，-1到-0.3是蓝到白，-0.3-0.3是白，0.3-1是从白到红
mycol <- c(colorRampPalette(c("#31688E","#FED9A6"))(100))
# mycol <- c(colorRampPalette(c('#F67B7B','#430052'))(100))
# mycol <- c(colorRampPalette(c('#20938C','#3B4992'))(100),
#            colorRampPalette(c('#3B4992','#FCD8D8'))(100),
#            colorRampPalette(c('#FCD8D8','#F35252'))(100))
# # colorRampPalette(c('white','white'))(3),
           # colorRampPalette(c('white','#772A26'))(100))
# max.rownames0 <- max(stri_length(rownames(metaboliter)))
# max.colnames0 <- max(stri_length(colnames(metaboliter)))

# mycol2 <- c(colorRampPalette(c('#9FA6CA','#F57373'))(100))


PLOT_DATA2 <- PLOT_DATA
PLOT_DATA2[,] <- "*"

annotation_row1 = as.data.frame(read.table(paste0("3、heatmap/anno.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
#ann_colors = c(P = "#772A26", H = "#1C4386")
ann_colors = list(
  Group = c( one= "#3E4989", two = "#386CB0", three = "#A6CEE3", four = "#CCEBC5")
)
PLOT_DATA
#画所有的图（不筛选r的），并且显示所有p小于0.05的*
PLOT_DATA<- PLOT_DATA[c(c("A","B","F","V",
                          "AB","AF","AV","BV","FV","BF",
                          "ABV","ABF","AFV","BFV",
                          "ABFV")),]
all10data_2 <- all10data_2[c(c("A","B","F","V",
                               "AB","AF","AV","BV","FV","BF",
                               "ABV","ABF","AFV","BFV",
                               "ABFV")),]
p=pheatmap(PLOT_DATA,
           #设置方格长宽
           cellwidth = 30,cellheight = 30,
           #设置行不聚类
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           #设置颜色
           color = mycol,
           #图例设置范围及
           legend_breaks = seq(-1,1,by=0.2),
           display_numbers = all10data_2,
           # display_numbers = all10data_2_text,
           fontsize_number = 9,#设置格子中显示的*大小
           #设置显示文字的颜色
           number_color = "white",
           #图片的标题
           # main = "correlation coefficient",
           #是否显示行列名
           #show_rownames=T,show_colnames=T,
           fontsize_row=14,#y轴名的字体大小
           fontsize_col=14,#x轴名的字体大小
           #angle_col=0#x轴名的角度
           annotation_row = annotation_row1,
           #annotation_row = annotation_row2,
           # annotation_col = annotation_row1,
           #annotation_names_row = TRUE, 
           #annotation_names_col = TRUE,
           annotation_colors = ann_colors
)
#ggraph::scale_edge_color_discrete(labels = function(x) str_wrap(x, 30) )
p

ggsave( filename =paste0("3、heatmap/heatmap.2.pdf"), p, width = 300, height = 200, units="mm",limitsize = FALSE)
ggsave(filename =paste0("3、heatmap/heatmap.2.png"), p, width = 300, height = 200, units="mm",limitsize = FALSE)
#?ggsave()

PLOT_DATA2<- PLOT_DATA2[c(c("A","B","F","V",
                            "AB","AF","AV","BV","FV","BF",
                            "ABV","ABF","AFV","BFV",
                            "ABFV")),]
#画所有的图（不筛选r的），并且显示所有p小于0.05的*
p2=pheatmap(PLOT_DATA,
           #设置方格长宽
           cellwidth = 30,cellheight = 30,
           #设置行不聚类
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           #设置颜色
           color = mycol,
           #图例设置范围及
           legend_breaks = seq(-1,1,by=0.2),
           display_numbers = PLOT_DATA2,
           # display_numbers = all10data_2_text,
           fontsize_number = 12,#设置格子中显示的*大小
           #设置显示文字的颜色
           number_color = "white",
           #图片的标题
           # main = "correlation coefficient",
           #是否显示行列名
           #show_rownames=T,show_colnames=T,
           fontsize_row=14,#y轴名的字体大小
           fontsize_col=14,#x轴名的字体大小
           #angle_col=0#x轴名的角度
           annotation_row = annotation_row1,
           #annotation_row = annotation_row2,
           # annotation_col = annotation_row1,
           #annotation_names_row = TRUE, 
           #annotation_names_col = TRUE,
           annotation_colors = ann_colors
)
#ggraph::scale_edge_color_discrete(labels = function(x) str_wrap(x, 30) )
p2

ggsave( filename =paste0("3、heatmap/heatmap.1.pdf"), p2, width = 300, height = 200, units="mm",limitsize = FALSE)
ggsave(filename =paste0("3、heatmap/heatmap.1.png"),  p2, width = 300, height = 200, units="mm",limitsize = FALSE)
#?ggsave()













#把特征矩阵加到一起分析一下-------------------------------

#先把15个表的20个合为一体
# out <- matrix(NA, 15, 3)          #一个9*111的空矩阵
all10data<-matrix(NA, nrow=1, ncol=3)

colnames(all10data) <- c("feature","importance","Group")

cl_list <- c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")


for(i in 1:length(cl_list)){
  library(readxl)
  feature_data <- read.csv(paste0(cl_list[i],"/feature_com_",cl_list[i],".csv"), header=T,  comment.char="",check.names=F,quot="",fileEncoding="utf8")

  all10data00 <- data.frame(feature=feature_data,
             Group = rep(cl_list[i], nrow(feature_data)))
  colnames(all10data00) <- c("feature","importance","Group")
  all10data <- rbind(all10data,all10data00)
  
}
all10data = all10data[-1,] 
all10data <- as.data.frame(all10data)

# 使用 gsub 函数删除引号
all10data$feature <- gsub('"', '', all10data$feature)
all10data$feature <- gsub("'", "", all10data$feature)

library(tidyr)

# 使用pivot_longer函数将宽表格转换为长表格
# df_long <- pivot_longer(all10data, 
#                         cols = starts_with("变量"), 
#                         names_to = "变量", 
#                         values_to = "值")

df_wide <- as.data.frame(spread(all10data, key = Group, value = importance ))

rownames(df_wide) <- df_wide$feature
df_wide_new <- df_wide[, -c(1)]

cl_list1 <-c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")
df_wide_new1 <- df_wide_new[,cl_list1]

cl_list2 <-c("A","B","F","V",
             "AB","AF","AV","BV","FV","BF",
             "ABV","ABF","AFV","BFV",
             "ABFV")
colnames(df_wide_new1) <- cl_list2


write.table(df_wide_new1, file=paste0("3、heatmap/feature_all.txt"), append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)


# 将所有非NA值转换为1
df_wide_new1[!is.na(df_wide_new1)] <- 1

# 将所有NA值转换为0
df_wide_new1[is.na(df_wide_new1)] <- 0

write.table(df_wide_new1, file=paste0("3、heatmap/feature_all_01.txt"), append = F, quote = F, sep="\t", eol = "\n", na = "NA", dec = ".", row.names = T, col.names = T)




#画图=====================================


df_wide_new1 <- as.matrix(read.table(paste0("3、heatmap/feature_all.txt"),header = T, row.names = 1,sep = '\t',quot="",check.names=F))

annotation_row1 = as.data.frame(read.table(paste0("3、heatmap/anno.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation_row1$Group<- factor(annotation_row1$Group,
                               levels = c("one","two","three","four"),ordered = TRUE)


annotation_row2 = as.data.frame(read.table(paste0("3、heatmap/feature_anno.txt"), header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F,quot=""))
annotation_row2$Kindom<- factor(annotation_row2$Kindom,
                                levels = c("Bacteria","Viruses","Archaea","Fungi"),ordered = TRUE)

#annotation_row2$importance<- factor(annotation_row2$importance,ordered = TRUE)
colnames(annotation_row1)

  
#ann_colors = c(P = "#772A26", H = "#1C4386")
ann_colors = list(
  Group = c( one= "#3E4989", two = "#386CB0", three = "#A6CEE3", four = "#CCEBC5"),
  Kindom = c( Bacteria = "#7570B3" , Viruses = "#BC80BD" , Archaea  = "#DECBE4",Fungi  =  "#FDDAEC"),# 没有fungi
  Change= c(up = "#DE77AE",down = "#7FBC41"),
  importance =   c(colorRampPalette(c('#3B4992FF','white'))(60),
                           colorRampPalette(c('white','white'))(29),
                           colorRampPalette(c('white','#EE0000FF'))(60)))
mycol <- c(colorRampPalette(c("#31688E","#FED9A6"))(100))

df_wide_new2 <- df_wide_new1[rownames(annotation_row2),]
# 将所有NA值转换为0
df_wide_new2[is.na(df_wide_new2)] <- 0
df_wide_new2<- df_wide_new2[,c("A","B","F","V",
                               "AB","AF","AV","BV","FV","BF",
                               "ABV","ABF","AFV","BFV",
                               "ABFV")]

library(pheatmap)
# rowSums(df_wide_new1)
# ?pheatmap()
p3=pheatmap(t(df_wide_new2),
            # scale = "column",
            #设置方格长宽
            cellwidth = 8,cellheight = 10,
            #设置行不聚类
            cluster_cols = FALSE,
            cluster_rows = FALSE,
            #设置颜色
            color = mycol,
            #图例设置范围及
            #legend_breaks = seq(-1,1,by=0.2),
            # display_numbers = PLOT_DATA2,
            # display_numbers = all10data_2_text,
            fontsize_number = 12,#设置格子中显示的*大小
            #设置显示文字的颜色
            number_color = "white",
            #图片的标题
            # main = "correlation coefficient",
            #是否显示行列名
            #show_rownames=T,show_colnames=T,
            fontsize_row=11,#y轴名的字体大小
            fontsize_col=8,#x轴名的字体大小
            #angle_col=0#x轴名的角度
            annotation_row = annotation_row1,
            #annotation_row = annotation_row2,
            annotation_col = annotation_row2,
            # annotation_names_row = TRUE,
            # annotation_names_col = TRUE,
            annotation_colors = ann_colors
)

library(ggplot2)

ggsave( filename =paste0("3、heatmap/heatmap.feature.pdf"), p3, width = 500, height = 150, units="mm",limitsize = FALSE)

ggsave(filename =paste0("3、heatmap/heatmap.feature.png"),  p3, width = 500, height = 150, units="mm",limitsize = FALSE)
#?ggsave()

